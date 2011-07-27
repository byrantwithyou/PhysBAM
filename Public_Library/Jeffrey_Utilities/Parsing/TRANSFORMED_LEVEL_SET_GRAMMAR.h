//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_TRANSFORMED_LEVEL_SET_GRAMMAR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_TRANSFORMED_LEVEL_SET_GRAMMAR_HPP

#include <boost/function.hpp>
#include <boost/spirit/home/qi/nonterminal/grammar.hpp>
#include <boost/spirit/home/qi/nonterminal/rule.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

struct TRANSFORMED_LEVEL_SET_GRAMMAR_BASE
{ static void Visit_Grammar_Elements(boost::function< void ( const char* ) > v); };

template< class T, int D, class T_ITERATOR, class T_SKIPPER >
struct TRANSFORMED_LEVEL_SET_GRAMMAR
    : TRANSFORMED_LEVEL_SET_GRAMMAR_BASE,
      boost::spirit::qi::grammar<
          T_ITERATOR,
          boost::function< T ( const VECTOR<T,D>& ) > ( ),
          T_SKIPPER
      >
{
    typedef boost::function< T ( const VECTOR<T,D>& ) > FUNCTION_TYPE;

    TRANSFORMED_LEVEL_SET_GRAMMAR();

    typedef boost::spirit::qi::rule< T_ITERATOR, FUNCTION_TYPE ( ), T_SKIPPER > RULE_TYPE;
    RULE_TYPE transformed_level_set;
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_TRANSFORMED_LEVEL_SET_GRAMMAR_HPP
