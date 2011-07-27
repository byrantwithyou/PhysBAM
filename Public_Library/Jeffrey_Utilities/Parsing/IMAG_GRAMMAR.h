//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_IMAG_GRAMMAR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_IMAG_GRAMMAR_HPP

#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/nonterminal/grammar.hpp>
#include <boost/spirit/home/qi/nonterminal/rule.hpp>
#include <boost/spirit/home/qi/numeric/real.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/operator/optional.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/string/lit.hpp>
#include <boost/spirit/home/support/context.hpp>

namespace PhysBAM
{

template< class T, class T_ITERATOR, class T_SKIPPER >
struct IMAG_GRAMMAR
    : boost::spirit::qi::grammar< T_ITERATOR, T ( char ), T_SKIPPER >
{
    IMAG_GRAMMAR()
        : IMAG_GRAMMAR::base_type(
              imag,
              "PhysBAm::IMAG_GRAMMAR< T, T_ITERATOR, T_SKIPPER >"
          )
    {
        namespace phx = boost::phoenix;
        namespace qi  = boost::spirit::qi;

        static const qi::real_parser<T> real;

        imag = (real >> -qi::lit('*') >> qi::lit(qi::_r1))
             | (qi::lit(qi::_r1) >> -qi::lit('*') >> real);
    }

    boost::spirit::qi::rule< T_ITERATOR, T ( char ), T_SKIPPER > imag;
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_IMAG_GRAMMAR_HPP
