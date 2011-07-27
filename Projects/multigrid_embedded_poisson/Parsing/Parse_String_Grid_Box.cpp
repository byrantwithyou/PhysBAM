//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef BOOST_PP_IS_ITERATING

#include <algorithm>
#include <string>

#include <boost/fusion/container/vector/vector10.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/control/expr_if.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/char/char_class.hpp>
#include <boost/spirit/home/qi/numeric/real.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/parse.hpp>

#include "Parse_String_Grid_Box.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

namespace
{

#define DECLARE_PARSE_STRING_GRID_BOX_DIM_HELPER( D ) \
template< class T > \
inline int Parse_String_Grid_Box_Dim_Helper( \
    const std::string& grid_box_str, \
    T (&min_x)[D], T (&max_x)[D]);
DECLARE_PARSE_STRING_GRID_BOX_DIM_HELPER( 2 )
DECLARE_PARSE_STRING_GRID_BOX_DIM_HELPER( 3 )
#undef DECLARE_PARSE_STRING_GRID_BOX_DIM_HELPER

} // namespace

template< class T, int D >
int Parse_String_Grid_Box(
    const std::string& grid_box_str,
    T (&min_x)[D], T (&max_x)[D])
{
    namespace qi = boost::spirit::qi;

    static const char char_of_int[] = { '0', '1', '2', '3' };
    static const char D_char = char_of_int[D];

#ifdef __GNUC__
    static const qi::real_parser<T> real = qi::real_parser<T>();
#else // #ifdef __GNUC__
    static const qi::real_parser<T> real;
#endif // #ifdef __GNUC__

    boost::fusion::vector2<T,T> interval_x;
    std::string::const_iterator it = grid_box_str.begin();
    const bool b = qi::phrase_parse(
        it, grid_box_str.end(),
        '[' >> real >> ',' >> real >> ']' >> '^' >> ('D' | qi::lit(D_char)),
        qi::standard::blank,
        interval_x
    );
    if(b) {
        if(it != grid_box_str.end())
            return 1;
        std::fill(&min_x[0], &min_x[D], boost::fusion::at_c<0>(interval_x));
        std::fill(&max_x[0], &max_x[D], boost::fusion::at_c<1>(interval_x));
        return 0;
    }

    return Parse_String_Grid_Box_Dim_Helper(grid_box_str, min_x, max_x);
}

namespace
{

#define _T_ref_comma_T_ref( z, n, data ) T&, T&
#define _min_x_n_comma_max_x_n( z, n, data ) min_x[n], max_x[n]
#define _if_n_x_real_comma_real( z, n, data ) \
        BOOST_PP_EXPR_IF( n, >> 'x' >> ) '[' >> real >> ',' >> real >> ']'

#define BOOST_PP_ITERATION_LIMITS ( 2, 3 )
#ifdef __GNUC__
// Likely requires adding -I. to your compile command.
#define BOOST_PP_FILENAME_1 "Parsing/Parse_String_Grid_Box.cpp"
#else // #ifdef __GNUC__
#define BOOST_PP_FILENAME_1 "Parse_String_Grid_Box.cpp"
#endif // #ifdef __GNUC__
#include BOOST_PP_ITERATE()

#undef _T_ref_comma_T_refs
#undef _min_x_n_comma_max_x_n
#undef _if_n_x_real_comma_real

} // namespace

#define EXPLICIT_INSTANTIATION( T, D ) \
template int Parse_String_Grid_Box<T,D>( \
    const std::string& grid_box_str, \
    T (&min_x)[D], T (&max_x)[D]);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#else // #ifndef BOOST_PP_IS_ITERATING

#define D BOOST_PP_ITERATION()

template< class T >
inline int Parse_String_Grid_Box_Dim_Helper(
    const std::string& grid_box_str,
    T (&min_x)[D], T (&max_x)[D])
{
    namespace qi = boost::spirit::qi;

#ifdef __GNUC__
    static const qi::real_parser<T> real = qi::real_parser<T>();
#else // #ifdef __GNUC__
    static const qi::real_parser<T> real;
#endif // #ifdef __GNUC__

    BOOST_PP_CAT( boost::fusion::vector, BOOST_PP_MUL( 2, D ) ) <
        BOOST_PP_ENUM( D, _T_ref_comma_T_ref, ~ )
    > box_x( BOOST_PP_ENUM( D, _min_x_n_comma_max_x_n, ~ ) );

    std::string::const_iterator it = grid_box_str.begin();
    const bool b = qi::phrase_parse(
        it, grid_box_str.end(),
        BOOST_PP_REPEAT( D, _if_n_x_real_comma_real, ~ ),
        qi::standard::blank,
        box_x
    );
    if(b && it == grid_box_str.end())
        return 0;

    return 1;
}

#endif // #ifndef BOOST_PP_IS_ITERATING
