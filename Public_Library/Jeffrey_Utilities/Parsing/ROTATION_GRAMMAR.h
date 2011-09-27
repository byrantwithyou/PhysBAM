//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_ROTATION_GRAMMAR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_ROTATION_GRAMMAR_HPP

#include <boost/spirit/home/phoenix/bind/bind_function.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/operator/arithmetic.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/qi/action/action.hpp>
#include <boost/spirit/home/qi/auxiliary/attr.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/nonterminal/grammar.hpp>
#include <boost/spirit/home/qi/nonterminal/rule.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/operator/optional.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/string/lit.hpp>
#include <boost/spirit/home/support/argument.hpp>
#include <boost/spirit/home/support/context.hpp>

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Parsing/RATIONAL_REAL_PARSER.h>
#include <Jeffrey_Utilities/Parsing/VECTOR_GRAMMAR.h>

namespace PhysBAM
{

template< class T, int D, class T_ITERATOR, class T_SKIPPER >
struct ROTATION_GRAMMAR;

template< class T, class T_ITERATOR, class T_SKIPPER >
struct ROTATION_GRAMMAR< T, 1, T_ITERATOR, T_SKIPPER >
    : boost::spirit::qi::grammar< T_ITERATOR, ROTATION< VECTOR<T,1> > ( ), T_SKIPPER >
{
    typedef ROTATION< VECTOR<T,1> > ROTATION_TYPE;

    ROTATION_GRAMMAR()
        : ROTATION_GRAMMAR::base_type(
              rotation,
              "PhysBAM::ROTATION_GRAMMAR< T, 1, T_ITERATOR, T_SKIPPER >"
          )
    {
        namespace qi = boost::spirit::qi;

        rotation = qi::attr(ROTATION_TYPE());
    }

    boost::spirit::qi::rule< T_ITERATOR, ROTATION_TYPE ( ), T_SKIPPER > rotation;
};

template< class T, class T_ITERATOR, class T_SKIPPER >
struct ROTATION_GRAMMAR< T, 2, T_ITERATOR, T_SKIPPER >
    : boost::spirit::qi::grammar< T_ITERATOR, ROTATION< VECTOR<T,2> > ( ), T_SKIPPER >
{
    typedef ROTATION< VECTOR<T,2> > ROTATION_TYPE;

    ROTATION_GRAMMAR()
        : ROTATION_GRAMMAR::base_type(
              rotation,
              "PhysBAM::ROTATION_GRAMMAR<T,2,T_ITERATOR>"
          )
    {
        namespace phx = boost::phoenix;
        namespace qi  = boost::spirit::qi;

        static const RATIONAL_REAL<T> real;

        static const VECTOR_GRAMMAR< T, 2, T_ITERATOR, T_SKIPPER > vector2;

        imag =
            (real >> -qi::lit('*') >> 'i')
          | ('i' >> -qi::lit('*') >> real);
        pm_imag =
            ('+' >> imag) [ qi::_val = qi::_1 ]
          | ('-' >> imag) [ qi::_val =-qi::_1 ];
        complex =
            (real >> pm_imag) [ qi::_val = phx::construct< COMPLEX<T> >(qi::_1, qi::_2) ]
          | (vector2)         [ qi::_val = phx::construct< COMPLEX<T> >(qi::_1) ];
        rotation =
            (complex) [ qi::_val = phx::bind(&ROTATION_TYPE::From_Complex, qi::_1) ]
          | (real)    [ qi::_val = phx::bind(&ROTATION_TYPE::From_Angle  , qi::_1) ];
    }

    boost::spirit::qi::rule< T_ITERATOR, T ( ), T_SKIPPER > imag;
    boost::spirit::qi::rule< T_ITERATOR, T ( ), T_SKIPPER > pm_imag;
    boost::spirit::qi::rule< T_ITERATOR, COMPLEX<T> ( ), T_SKIPPER > complex;
    boost::spirit::qi::rule< T_ITERATOR, ROTATION_TYPE ( ), T_SKIPPER > rotation;
};

template< class T, class T_ITERATOR, class T_SKIPPER >
struct ROTATION_GRAMMAR< T, 3, T_ITERATOR, T_SKIPPER >
    : boost::spirit::qi::grammar< T_ITERATOR, ROTATION< VECTOR<T,3> > ( ), T_SKIPPER >
{
    typedef ROTATION< VECTOR<T,3> > ROTATION_TYPE;

    ROTATION_GRAMMAR()
        : ROTATION_GRAMMAR::base_type(
              rotation,
              "PhysBAM::ROTATION_GRAMMAR<T,3,T_ITERATOR>"
          )
    {
        namespace phx = boost::phoenix;
        namespace qi  = boost::spirit::qi;

        static const RATIONAL_REAL<T> real;

        static const VECTOR_GRAMMAR< T, 3, T_ITERATOR, T_SKIPPER > vector3;
        static const VECTOR_GRAMMAR< T, 4, T_ITERATOR, T_SKIPPER > vector4;

        angle_direction =
            (real >> ',' >> vector3) [ qi::_val = phx::construct< ROTATION_TYPE >(qi::_1, qi::_2) ];
        direction_angle =
            (vector3 >> ',' >> real) [ qi::_val = phx::construct< ROTATION_TYPE >(qi::_2, qi::_1) ];

        imag =
            (real >> -qi::lit('*') >> qi::lit(qi::_r1))
          | (qi::lit(qi::_r1) >> -qi::lit('*') >> real);
        pm_imag =
            ('+' >> imag(qi::_r1)) [ qi::_val = qi::_1 ]
          | ('-' >> imag(qi::_r1)) [ qi::_val =-qi::_1 ];

        imag_vector3 =
            (vector3 >> -qi::lit('*') >> 'i')
          | ('i' >> -qi::lit('*') >> vector3);
        pm_imag_vector3 =
            ('+' >> imag_vector3) [ qi::_val = qi::_1 ]
          | ('-' >> imag_vector3) [ qi::_val =-qi::_1 ];

        quaternion =
            (real >> pm_imag('i') >> pm_imag('j') >> pm_imag('k'))
            [ qi::_val = phx::construct< QUATERNION<T> >(qi::_1, qi::_2, qi::_3, qi::_4) ]
          | (real >> pm_imag_vector3)
            [ qi::_val = phx::construct< QUATERNION<T> >(qi::_1, qi::_2) ]
          | (vector4)
            [ qi::_val = phx::construct< QUATERNION<T> >(qi::_1) ];

        rotation =
            (angle_direction) [ qi::_val = qi::_1 ]
          | (direction_angle) [ qi::_val = qi::_1 ]
          | (quaternion)      [ qi::_val = phx::bind(&ROTATION_TYPE::From_Quaternion     , qi::_1) ]
          | (vector3)         [ qi::_val = phx::bind(&ROTATION_TYPE::From_Rotation_Vector, qi::_1) ];
    }

    boost::spirit::qi::rule< T_ITERATOR, ROTATION_TYPE ( ), T_SKIPPER > angle_direction;
    boost::spirit::qi::rule< T_ITERATOR, ROTATION_TYPE ( ), T_SKIPPER > direction_angle;
    boost::spirit::qi::rule< T_ITERATOR, T ( char ), T_SKIPPER > imag;
    boost::spirit::qi::rule< T_ITERATOR, T ( char ), T_SKIPPER > pm_imag;
    boost::spirit::qi::rule< T_ITERATOR, VECTOR<T,3> ( ), T_SKIPPER > imag_vector3;
    boost::spirit::qi::rule< T_ITERATOR, VECTOR<T,3> ( ), T_SKIPPER > pm_imag_vector3;
    boost::spirit::qi::rule< T_ITERATOR, QUATERNION<T> ( ), T_SKIPPER > quaternion;
    boost::spirit::qi::rule< T_ITERATOR, ROTATION_TYPE ( ), T_SKIPPER > rotation;
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_ROTATION_GRAMMAR_HPP
