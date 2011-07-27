//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_IMAG_PARSER_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_IMAG_PARSER_HPP

#include <boost/mpl/assert.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/domain.hpp>
#include <boost/spirit/home/qi/operator/optional.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/parse.hpp>
#include <boost/spirit/home/qi/parser.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/spirit/home/support/info.hpp>
#include <boost/spirit/home/support/meta_compiler.hpp>

namespace PhysBAM
{

template< class Subject >
struct IMAG_PARSER
    : boost::spirit::qi::unary_parser< IMAG_PARSER< Subject > >
{
    BOOST_MPL_ASSERT((boost::spirit::traits::is_parser< Subject >));

    typedef Subject subject_type;

    template< class Context, class Iterator >
    struct attribute
        : boost::spirit::traits::attribute_of< Subject, Context, Iterator >
    { };

    IMAG_PARSER(const Subject& subject_, const char i_)
        : subject(subject_), i(i_)
    { }

    template< class Iterator, class Context, class Skipper, class Attribute >
    bool parse(
        Iterator& first, const Iterator& last,
        Context& context,
        const Skipper& skipper,
        Attribute& attribute) const
    {
        namespace qi = boost::spirit::qi;

        Iterator it = first;
        if(
            (subject.parse(it, last, context, skipper, attribute)
          && qi::phrase_parse(it, last, -qi::lit('*') >> qi::lit(i), skipper, qi::skip_flag::dont_postskip))
         || (qi::phrase_parse(it = first, last, qi::lit(i) >> -qi::lit('*'), skipper, qi::skip_flag::dont_postskip)
          && subject.parse(it, last, context, skipper, attribute))
        ) {
            first = it;
            return true;
        }
        return false;
    }

    template< class Context >
    boost::spirit::info what(const Context& context) const
    { return boost::spirit::info("IMAG_PARSER", subject.what(context)); }

    Subject subject;
    char i;
};

namespace Result_Of
{

template< class Subject >
struct MAKE_IMAG_PARSER
{
    typedef IMAG_PARSER<
        typename boost::spirit::result_of::compile<
            boost::spirit::qi::domain,
            Subject
        >::type
    > type;
};

} // namespace Result_Of

template< class Subject >
inline typename Result_Of::MAKE_IMAG_PARSER< Subject >::type
Make_Imag_Parser(const Subject& subject, const char i)
{
    return typename Result_Of::MAKE_IMAG_PARSER< Subject >::type(
        boost::spirit::compile< boost::spirit::qi::domain >(subject),
        i
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_IMAG_PARSER_HPP
