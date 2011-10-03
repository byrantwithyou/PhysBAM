#if 0
//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef BOOST_PP_IS_ITERATING

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_VECTOR_GRAMMAR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_VECTOR_GRAMMAR_HPP

#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/control/expr_if.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum_shifted_params.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/qi/action/action.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/nonterminal/grammar.hpp>
#include <boost/spirit/home/qi/nonterminal/rule.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/support/argument.hpp>
#include <boost/spirit/home/support/context.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Parsing/RATIONAL_REAL_PARSER.h>

namespace PhysBAM
{

template< class T, int D, class T_ITERATOR, class T_SKIPPER >
struct VECTOR_GRAMMAR;

#define _if_n_comma_real( z, n, data ) BOOST_PP_EXPR_IF( n, >> ',' >> ) real

#define BOOST_PP_ITERATION_LIMITS ( 1, 4 )
#define BOOST_PP_FILENAME_1 <Jeffrey_Utilities/Parsing/VECTOR_GRAMMAR.h>
#include BOOST_PP_ITERATE()

#undef _if_n_comma_real

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_VECTOR_GRAMMAR_HPP

#else // #ifndef BOOST_PP_IS_ITERATING

#define D BOOST_PP_ITERATION()

template< class T, class T_ITERATOR, class T_SKIPPER >
struct VECTOR_GRAMMAR< T, D, T_ITERATOR, T_SKIPPER >
    : boost::spirit::qi::grammar< T_ITERATOR, VECTOR<T,D> ( ), T_SKIPPER >
{
    VECTOR_GRAMMAR()
        : VECTOR_GRAMMAR::base_type(
              vector,
              "PhysBAM::VECTOR_GRAMMAR< T, " BOOST_PP_STRINGIZE( D ) ", T_ITERATOR, T_SKIPPER >"
          )
    {
        namespace phx = boost::phoenix;
        namespace qi  = boost::spirit::qi;

        static const RATIONAL_REAL<T> real;

        vector =
            ('(' >> BOOST_PP_REPEAT( D, _if_n_comma_real, ~ ) >> ')' )
            [
                qi::_val = phx::construct< VECTOR<T,D> >(
                    BOOST_PP_ENUM_SHIFTED_PARAMS( BOOST_PP_INC( D ), qi::_ )
                )
            ];
    }

    boost::spirit::qi::rule< T_ITERATOR, VECTOR<T,D> ( ), T_SKIPPER > vector;
};

#endif // #ifndef BOOST_PP_IS_ITERATING
#endif
