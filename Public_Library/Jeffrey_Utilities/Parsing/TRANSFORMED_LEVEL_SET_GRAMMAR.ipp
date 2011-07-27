//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_TRANSFORMED_LEVEL_SET_GRAMMAR_IPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_TRANSFORMED_LEVEL_SET_GRAMMAR_IPP

#include <boost/function.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/operator/arithmetic.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/qi/action/action.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/string/lit.hpp>
#include <boost/spirit/home/support/argument.hpp>
#include <boost/spirit/home/support/context.hpp>

#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/NEGATE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/ROTATE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SCALE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/TRANSLATE_FUNCTION.h>
#include <Jeffrey_Utilities/Parsing/BASIC_LEVEL_SET_GRAMMAR.h>
#include <Jeffrey_Utilities/Parsing/RATIONAL_REAL_PARSER.h>
#include <Jeffrey_Utilities/Parsing/ROTATION_GRAMMAR.h>
#include <Jeffrey_Utilities/Parsing/VECTOR_GRAMMAR.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include <Jeffrey_Utilities/Parsing/TRANSFORMED_LEVEL_SET_GRAMMAR.h>

namespace PhysBAM
{

void
TRANSFORMED_LEVEL_SET_GRAMMAR_BASE::
Visit_Grammar_Elements(boost::function< void ( const char* ) > v)
{
    v("negate(<level-set>)");
    v("rotate[<rot>](<level-set>)");
    v("scale[<c>](<level-set>)");
    v("translate[<x0>](<level-set>)");
}

template< class T, int D, class T_ITERATOR, class T_SKIPPER >
TRANSFORMED_LEVEL_SET_GRAMMAR< T, D, T_ITERATOR, T_SKIPPER >::
TRANSFORMED_LEVEL_SET_GRAMMAR()
    : TRANSFORMED_LEVEL_SET_GRAMMAR::base_type(
          transformed_level_set,
          "PhysBAM::TRANSFORMED_LEVEL_SET_GRAMMAR< T, D, T_ITERATOR, T_SKIPPER >"
      )
{
    namespace phx = boost::phoenix;
    namespace qi  = boost::spirit::qi;

    static const RATIONAL_REAL<T> real;

    static const VECTOR_GRAMMAR< T, D, T_ITERATOR, T_SKIPPER > vector;
    static const ROTATION_GRAMMAR< T, D, T_ITERATOR, T_SKIPPER > rotation;

    static const BASIC_LEVEL_SET_GRAMMAR< T, D, T_ITERATOR, T_SKIPPER > basic_level_set;

    typedef COMPOSE_FUNCTION< NEGATE_FUNCTION, FUNCTION_TYPE > POST_NEGATE_FUNCTION_TYPE;
    RULE_TYPE negate =
        ("negate(" >> transformed_level_set >> ')')
        [ qi::_val = phx::construct< POST_NEGATE_FUNCTION_TYPE >(NEGATE_FUNCTION(), qi::_1) ];

    typedef ROTATION< VECTOR<T,D> > ROTATION_TYPE;
    typedef ROTATE_FUNCTION< ROTATION_TYPE > ROTATE_FUNCTION_TYPE;
    typedef COMPOSE_FUNCTION< FUNCTION_TYPE, ROTATE_FUNCTION_TYPE > PRE_ROTATE_FUNCTION_TYPE;
    RULE_TYPE rotate =
        ("rotate[" >> rotation >> "](" >> transformed_level_set >> ')')
        [
            qi::_val = phx::construct< PRE_ROTATE_FUNCTION_TYPE >(
                qi::_2,
                phx::construct< ROTATE_FUNCTION_TYPE >(
                    phx::bind(&ROTATION_TYPE::Inverse, qi::_1)
                )
            )
        ];

    typedef SCALE_FUNCTION<T> SCALE_BY_SCALAR_FUNCTION_TYPE;
    typedef SCALE_FUNCTION< VECTOR<T,D> > SCALE_BY_VECTOR_FUNCTION_TYPE;
    typedef COMPOSE_FUNCTION< FUNCTION_TYPE, SCALE_BY_SCALAR_FUNCTION_TYPE > PRE_SCALE_BY_SCALAR_FUNCTION_TYPE;
    typedef COMPOSE_FUNCTION< FUNCTION_TYPE, SCALE_BY_VECTOR_FUNCTION_TYPE > PRE_SCALE_BY_VECTOR_FUNCTION_TYPE;
    typedef COMPOSE_FUNCTION< SCALE_FUNCTION<T>, PRE_SCALE_BY_SCALAR_FUNCTION_TYPE > POST_PRE_SCALE_BY_SCALAR_FUNCTION_TYPE;
    RULE_TYPE scale = "scale[" >> (
        (real >> "](" >> transformed_level_set >> ')')
        [
            qi::_val = phx::construct< POST_PRE_SCALE_BY_SCALAR_FUNCTION_TYPE >(
                phx::construct< SCALE_FUNCTION<T> >(qi::_1),
                phx::construct< PRE_SCALE_BY_SCALAR_FUNCTION_TYPE >(
                    qi::_2,
                    phx::construct< SCALE_BY_SCALAR_FUNCTION_TYPE >(static_cast<T>(1) / qi::_1)
                )
            )
        ]
      | (vector >> "](" >> transformed_level_set >> ')')
        [
            qi::_val = phx::construct< PRE_SCALE_BY_VECTOR_FUNCTION_TYPE >(
                qi::_2,
                phx::construct< SCALE_BY_VECTOR_FUNCTION_TYPE >(static_cast<T>(1) / qi::_1)
            )
        ]
    );

    typedef TRANSLATE_FUNCTION<T> TRANSLATE_BY_SCALAR_FUNCTION_TYPE;
    typedef TRANSLATE_FUNCTION< VECTOR<T,D> > TRANSLATE_BY_VECTOR_FUNCTION_TYPE;
    typedef COMPOSE_FUNCTION< FUNCTION_TYPE, TRANSLATE_BY_SCALAR_FUNCTION_TYPE > PRE_TRANSLATE_BY_SCALAR_FUNCTION_TYPE;
    typedef COMPOSE_FUNCTION< FUNCTION_TYPE, TRANSLATE_BY_VECTOR_FUNCTION_TYPE > PRE_TRANSLATE_BY_VECTOR_FUNCTION_TYPE;
    RULE_TYPE translate = "translate[" >> (
        (real >> "](" >> transformed_level_set >> ')')
        [
            qi::_val = phx::construct< PRE_TRANSLATE_BY_SCALAR_FUNCTION_TYPE >(
                qi::_2,
                phx::construct< TRANSLATE_BY_SCALAR_FUNCTION_TYPE >(-qi::_1)
            )
        ]
      | (vector >> "](" >> transformed_level_set >> ')')
        [
            qi::_val = phx::construct< PRE_TRANSLATE_BY_VECTOR_FUNCTION_TYPE >(
                qi::_2,
                phx::construct< TRANSLATE_BY_VECTOR_FUNCTION_TYPE >(-qi::_1)
            )
        ]
    );

    transformed_level_set = negate.copy()
                          | rotate.copy()
                          | scale.copy()
                          | translate.copy()
                          | basic_level_set;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_TRANSFORMED_LEVEL_SET_GRAMMAR_IPP
