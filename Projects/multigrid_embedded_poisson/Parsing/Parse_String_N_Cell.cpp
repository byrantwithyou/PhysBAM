//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef BOOST_PP_IS_ITERATING

#include <algorithm>
#include <string>

#include <boost/fusion/container/vector/vector10.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/facilities/intercept.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/char/char_class.hpp>
#include <boost/spirit/home/qi/numeric/uint.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/parse.hpp>

#include "Parse_String_N_Cell.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace
{

#define DECLARE_PARSE_STRING_N_CELL_DIM_HELPER( D ) \
inline int Parse_String_N_Cell_Dim_Helper( \
    const std::string& n_cell_str, \
    unsigned int (&n_cell)[D]);
DECLARE_PARSE_STRING_N_CELL_DIM_HELPER( 2 )
DECLARE_PARSE_STRING_N_CELL_DIM_HELPER( 3 )
#undef DECLARE_PARSE_STRING_N_CELL_DIM_HELPER

} // namespace

template< int D >
int Parse_String_N_Cell(
    const std::string& n_cell_str,
    unsigned int (&n_cell)[D])
{
    namespace qi = boost::spirit::qi;

    static const char char_of_int[] = { '0', '1', '2', '3' };
    static const char D_char = char_of_int[D];

    std::string::const_iterator it = n_cell_str.begin();
    const bool b = qi::phrase_parse(
        it, n_cell_str.end(),
        qi::uint_ >> '^' >> ('D' | qi::lit(D_char)),
        qi::standard::blank,
        n_cell[0]
    );
    if(b) {
        if(it != n_cell_str.end())
            return 1;
        std::fill(&n_cell[1], &n_cell[D], n_cell[0]);
        return 0;
    }

    return Parse_String_N_Cell_Dim_Helper(n_cell_str, n_cell);
}

namespace
{

#define _n_cell_n( z, n, data ) n_cell[n]
#define _if_n_x_qi_uint( z, n, data ) BOOST_PP_EXPR_IF( n, >> 'x' >> ) qi::uint_

#define BOOST_PP_ITERATION_LIMITS ( 2, 3 )
#ifdef __GNUC__
// Likely requires adding -I. to your compile command.
#define BOOST_PP_FILENAME_1 "Parsing/Parse_String_N_Cell.cpp"
#else // #ifdef __GNUC__
#define BOOST_PP_FILENAME_1 "Parse_String_N_Cell.cpp"
#endif // #ifdef __GNUC__
#include BOOST_PP_ITERATE()

#undef _n_cell_n
#undef _if_n_x_qi_uint

} // namespace

#define EXPLICIT_INSTANTIATION( D ) \
template int Parse_String_N_Cell<D>( \
    const std::string& n_cell_str, \
    unsigned (&n_cell)[D]);
EXPLICIT_INSTANTIATION( 2 )
EXPLICIT_INSTANTIATION( 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#else // #ifndef BOOST_PP_IS_ITERATING

#define D BOOST_PP_ITERATION()

inline int Parse_String_N_Cell_Dim_Helper(
    const std::string& n_cell_str,
    unsigned int (&n_cell)[D])
{
    namespace qi = boost::spirit::qi;

    BOOST_PP_CAT( boost::fusion::vector, D ) <
        BOOST_PP_ENUM_PARAMS( D, unsigned int& BOOST_PP_INTERCEPT )
    > _n_cell( BOOST_PP_ENUM( D, _n_cell_n, ~ ) );

    std::string::const_iterator it = n_cell_str.begin();
    const bool b = qi::phrase_parse(
        it, n_cell_str.end(),
        BOOST_PP_REPEAT( D, _if_n_x_qi_uint, ~ ),
        qi::standard::blank,
        _n_cell
    );
    if(b && it == n_cell_str.end())
        return 0;

    return 1;
}

#endif // #ifndef BOOST_PP_IS_ITERATING
