//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef BOOST_PP_IS_ITERATING

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_BASIC_LEVEL_SET_GRAMMAR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_BASIC_LEVEL_SET_GRAMMAR_HPP

#include <boost/function.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/spirit/home/qi/nonterminal/grammar.hpp>
#include <boost/spirit/home/qi/nonterminal/rule.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int D >
struct BASIC_LEVEL_SET_GRAMMAR_BASE;

template< class T, int D, class T_ITERATOR, class T_SKIPPER >
struct BASIC_LEVEL_SET_GRAMMAR;

#define BOOST_PP_ITERATION_LIMITS ( 1, 3 )
#define BOOST_PP_FILENAME_1 <Jeffrey_Utilities/Parsing/BASIC_LEVEL_SET_GRAMMAR.h>
#include BOOST_PP_ITERATE()

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_BASIC_LEVEL_SET_GRAMMAR_HPP

#else // #ifndef BOOST_PP_IS_ITERATING

#define D BOOST_PP_ITERATION()

template<>
struct BASIC_LEVEL_SET_GRAMMAR_BASE<D>
{ static void Visit_Grammar_Elements(boost::function< void ( const char* ) > v); };

template< class T, class T_ITERATOR, class T_SKIPPER >
struct BASIC_LEVEL_SET_GRAMMAR< T, D, T_ITERATOR, T_SKIPPER >
    : BASIC_LEVEL_SET_GRAMMAR_BASE<D>,
      boost::spirit::qi::grammar<
          T_ITERATOR,
          boost::function< T ( const VECTOR<T,D>& ) > ( ),
          T_SKIPPER
      >
{
    typedef boost::function< T ( const VECTOR<T,D>& ) > FUNCTION_TYPE;

    BASIC_LEVEL_SET_GRAMMAR();

    typedef boost::spirit::qi::rule< T_ITERATOR, FUNCTION_TYPE ( ), T_SKIPPER > RULE_TYPE;
    RULE_TYPE basic_level_set;
};

#undef D

#endif // #ifndef BOOST_PP_IS_ITERATING
