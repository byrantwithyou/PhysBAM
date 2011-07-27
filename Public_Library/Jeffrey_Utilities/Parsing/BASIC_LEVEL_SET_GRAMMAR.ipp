//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_BASIC_LEVEL_SET_GRAMMAR_IPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_BASIC_LEVEL_SET_GRAMMAR_IPP

#include <boost/function.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/qi/action/action.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/string/lit.hpp>
#include <boost/spirit/home/support/argument.hpp>
#include <boost/spirit/home/support/context.hpp>

#include <Jeffrey_Utilities/Level_Sets/STAR24_LEVEL_SET.h>
#include <Jeffrey_Utilities/Level_Sets/TORUS_LEVEL_SET.h>
#include <Jeffrey_Utilities/Level_Sets/TREFOIL_KNOT_LEVEL_SET.h>
#include <Jeffrey_Utilities/Level_Sets/UNIT_CUBE_LEVEL_SET.h>
#include <Jeffrey_Utilities/Level_Sets/UNIT_SPHERE_LEVEL_SET.h>
#include <Jeffrey_Utilities/Parsing/RATIONAL_REAL_PARSER.h>

#include <Jeffrey_Utilities/Parsing/BASIC_LEVEL_SET_GRAMMAR.h>

namespace PhysBAM
{

void
BASIC_LEVEL_SET_GRAMMAR_BASE<1>::
Visit_Grammar_Elements(boost::function< void ( const char* ) > v)
{
    v("unit_interval | unit_cube | unit_sphere");
}

template< class T, class T_ITERATOR, class T_SKIPPER >
BASIC_LEVEL_SET_GRAMMAR< T, 1, T_ITERATOR, T_SKIPPER >::
BASIC_LEVEL_SET_GRAMMAR()
    : BASIC_LEVEL_SET_GRAMMAR::base_type(
          basic_level_set,
          "PhysBAM::BASIC_LEVEL_SET_GRAMMAR< T, 1, T_ITERATOR >"
      )
{
    namespace phx = boost::phoenix;
    namespace qi = boost::spirit::qi;

#ifdef __GNUC__
    static const RATIONAL_REAL<T> real = RATIONAL_REAL<T>();
#else // #ifdef __GNUC__
    static const RATIONAL_REAL<T> real;
#endif // #ifdef __GNUC__

    RULE_TYPE unit_interval =
        (qi::lit("unit_interval") | "unit_cube" | "unit_sphere")
        [ qi::_val = UNIT_SPHERE_LEVEL_SET() ];

    basic_level_set = unit_interval.copy();
}

void
BASIC_LEVEL_SET_GRAMMAR_BASE<2>::
Visit_Grammar_Elements(boost::function< void ( const char* ) > v)
{
    v("unit_circle | unit_sphere");
    v("unit_square | unit_cube");
}

template< class T, class T_ITERATOR, class T_SKIPPER >
BASIC_LEVEL_SET_GRAMMAR< T, 2, T_ITERATOR, T_SKIPPER >::
BASIC_LEVEL_SET_GRAMMAR()
    : BASIC_LEVEL_SET_GRAMMAR::base_type(
          basic_level_set,
          "PhysBAM::BASIC_LEVEL_SET_GRAMMAR< T, 2, T_ITERATOR >"
      )
{
    namespace phx = boost::phoenix;
    namespace qi = boost::spirit::qi;

#ifdef __GNUC__
    static const RATIONAL_REAL<T> real = RATIONAL_REAL<T>();
#else // #ifdef __GNUC__
    static const RATIONAL_REAL<T> real;
#endif // #ifdef __GNUC__

    RULE_TYPE unit_circle =
        (qi::lit("unit_circle") | "unit_sphere")
        [ qi::_val = UNIT_SPHERE_LEVEL_SET() ];
    RULE_TYPE unit_square =
        (qi::lit("unit_square") | "unit_cube")
        [ qi::_val = UNIT_CUBE_LEVEL_SET() ];

    basic_level_set = unit_circle.copy()
                    | unit_square.copy();
}

void
BASIC_LEVEL_SET_GRAMMAR_BASE<3>::
Visit_Grammar_Elements(boost::function< void ( const char* ) > v)
{
    v("star24(<min-radius>{,<max-radius>})");
    v("torus(<minor-radius>)");
    v("trefoil_knot(<minor-radius>)");
    v("unit_cube");
    v("unit_sphere");
}

template< class T, class T_ITERATOR, class T_SKIPPER >
BASIC_LEVEL_SET_GRAMMAR< T, 3, T_ITERATOR, T_SKIPPER >::
BASIC_LEVEL_SET_GRAMMAR()
    : BASIC_LEVEL_SET_GRAMMAR::base_type(
          basic_level_set,
          "PhysBAM::BASIC_LEVEL_SET_GRAMMAR< T, 3, T_ITERATOR >"
      )
{
    namespace phx = boost::phoenix;
    namespace qi  = boost::spirit::qi;

#ifdef __GNUC__
    static const RATIONAL_REAL<T> real = RATIONAL_REAL<T>();
#else // #ifdef __GNUC__
    static const RATIONAL_REAL<T> real;
#endif // #ifdef __GNUC__

    RULE_TYPE star24 = "star24(" >> (
        (real >> ')')
        [ qi::_val = phx::construct< STAR24_LEVEL_SET<T> >(qi::_1) ]
      | (real >> ',' >> real >> ')')
        [ qi::_val = phx::construct< STAR24_LEVEL_SET<T> >(qi::_1, qi::_2) ]
    );
    RULE_TYPE torus =
        ("torus(" >> real >> ')')
        [ qi::_val = phx::construct< TORUS_LEVEL_SET<T> >(qi::_1) ];
    RULE_TYPE trefoil_knot =
        ("trefoil_knot(" >> real >> ')')
        [ qi::_val = phx::construct< TREFOIL_KNOT_LEVEL_SET<T> >(qi::_1) ];
    RULE_TYPE unit_cube =
        qi::lit("unit_cube")
        [ qi::_val = UNIT_CUBE_LEVEL_SET() ];
    RULE_TYPE unit_sphere =
        qi::lit("unit_sphere")
        [ qi::_val = UNIT_SPHERE_LEVEL_SET() ];

    basic_level_set = star24.copy()
                    | torus.copy()
                    | trefoil_knot.copy()
                    | unit_cube.copy()
                    | unit_sphere.copy();
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_BASIC_LEVEL_SET_GRAMMAR_IPP
