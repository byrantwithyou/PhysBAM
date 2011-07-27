//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// RATIONAL_REAL<T> is similar to boost::spirit::qi::real_parser<T>, except it
// can additionally parse quotients, such as "1/3".
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_RATIONAL_REAL_PARSER_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_RATIONAL_REAL_PARSER_HPP

#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/spirit/home/phoenix/core/reference.hpp>
#include <boost/spirit/home/phoenix/operator/arithmetic.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/qi/action/action.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/detail/assign_to.hpp>
#include <boost/spirit/home/qi/domain.hpp>
#include <boost/spirit/home/qi/meta_compiler.hpp>
#include <boost/spirit/home/qi/numeric/real.hpp>
#include <boost/spirit/home/qi/numeric/real_policies.hpp>
#include <boost/spirit/home/qi/operator/optional.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/parse.hpp>
#include <boost/spirit/home/qi/parser.hpp>
#include <boost/spirit/home/qi/skip_over.hpp>
#include <boost/spirit/home/support/argument.hpp>
#include <boost/spirit/home/support/info.hpp>
#include <boost/spirit/home/support/meta_compiler.hpp>
#include <boost/spirit/home/support/terminal.hpp>
#include <boost/type_traits/integral_constant.hpp>

namespace PhysBAM
{

// Define a tag identifying RATIONAL_REAL<T> terminal objects.
namespace Tag { template< class T > struct RATIONAL_REAL { }; }

// Define the RATIONAL_REAL<T> terminal.
template< class T, class T_POLICIES = boost::spirit::qi::real_policies<T> >
struct RATIONAL_REAL
    : boost::spirit::terminal<
          boost::spirit::tag::stateful_tag< T_POLICIES, Tag::RATIONAL_REAL<T> >
      >
{
    typedef boost::spirit::terminal<
        boost::spirit::tag::stateful_tag< T_POLICIES, Tag::RATIONAL_REAL<T> >
    > terminal_;
    RATIONAL_REAL()
    { }
    RATIONAL_REAL(const T_POLICIES& policies)
        : terminal_(policies)
    { }
};

// Define RATIONAL_REAL_PARSER<T> and LITERAL_RATIONAL_REAL_PARSER<T>, the
// parsers behind RATIONAL_REAL<T> terminal objects.
template< class T, class T_POLICIES = boost::spirit::qi::real_policies<T> >
struct RATIONAL_REAL_PARSER
    : boost::spirit::qi::primitive_parser< RATIONAL_REAL_PARSER<T> >
{
    template< class Context, class Iterator >
    struct attribute
    { typedef T type; };

    template< class Iterator, class Context, class Skipper >
    bool parse(
        Iterator& first, const Iterator& last,
        Context& /*context*/, const Skipper& skipper,
        T& attribute) const
    {
        namespace phx = boost::phoenix;
        namespace qi  = boost::spirit::qi;
#ifdef __GNUC__
        static const qi::real_parser< T, T_POLICIES > real = qi::real_parser< T, T_POLICIES >();
#else // #ifdef __GNUC__
        static const qi::real_parser< T, T_POLICIES > real;
#endif // #ifdef __GNUC__
        qi::skip_over(first, last, skipper);
        return qi::parse(
            first, last,
            real [ phx::ref(attribute) = qi::_1 ] >>
            -('/' >> real) [ phx::ref(attribute) /= qi::_1 ]
        );
    }

    template< class Iterator, class Context, class Skipper, class Attribute >
    bool parse(
        Iterator& first, const Iterator& last,
        Context& context, const Skipper& skipper,
        Attribute& attribute) const
    {
        T temp;
        if(!parse(first, last, context, skipper, temp))
            return false;
        boost::spirit::traits::assign_to(temp, attribute);
        return true;
    }

    template< class Context >
    boost::spirit::info what(const Context& /*context*/) const
    { return boost::spirit::info("RATIONAL_REAL"); }
};

template< class T, class T_POLICIES = boost::spirit::qi::real_policies<T> >
struct LITERAL_RATIONAL_REAL_PARSER
    : boost::spirit::qi::primitive_parser< LITERAL_RATIONAL_REAL_PARSER< T, T_POLICIES > >
{
    template< class U >
    explicit LITERAL_RATIONAL_REAL_PARSER(const U& x)
        : m_x(x)
    { }

    template< class Context, class Iterator >
    struct attribute
    { typedef T type; };

    template< class Iterator, class Context, class Skipper, class Attribute >
    bool parse(
        Iterator& first, const Iterator& last,
        Context& context, const Skipper& skipper,
        Attribute& attribute) const
    {
        T temp;
        const bool parse_result = RATIONAL_REAL_PARSER< T, T_POLICIES >().parse(
            first, last,
            context, skipper,
            temp
        );
        if(!parse_result || !(temp == m_x))
            return false;
        boost::spirit::traits::assign_to(temp, attribute);
        return true;
    }

    template< class Context >
    boost::spirit::info what(const Context& /*context*/) const
    { return boost::spirit::info("LITERAL_RATIONAL_REAL"); }

private:
    const T m_x;
};

} // namespace PhysBAM

namespace boost { namespace spirit {

// Make Boost.Spirit aware of PhysBAM::RATIONAL_REAL<T> terminal objects; these
// are identifiable by their PhysBAM::Tag::RATIONAL_REAL<T> tags.
template< class Policies, class T >
struct use_terminal<
    qi::domain,
    tag::stateful_tag< Policies, PhysBAM::Tag::RATIONAL_REAL<T> >
>
    : boost::true_type
{ };

template< class Policies, class T, class A0 >
struct use_terminal<
    qi::domain,
    terminal_ex<
        tag::stateful_tag< Policies, PhysBAM::Tag::RATIONAL_REAL<T> >,
        boost::fusion::vector1< A0 >
    >
>
    : boost::true_type
{ };

template< class Policies, class T >
struct use_lazy_terminal<
    qi::domain,
    tag::stateful_tag< Policies, PhysBAM::Tag::RATIONAL_REAL<T> >,
    1
>
    : boost::true_type
{ };

namespace qi
{

// Define the factory functions used by the Boost.Spirit[.Qi] meta-compiler to
// transform PhysBAM::RATIONAL_REAL<T> terminal objects into
// PhysBAM::RATIONAL_REAL_PARSER<T> parser objects.
template< class Policies, class T, class Modifiers >
struct make_primitive<
    tag::stateful_tag< Policies, PhysBAM::Tag::RATIONAL_REAL<T> >,
    Modifiers
>
{
    typedef PhysBAM::RATIONAL_REAL_PARSER< T, Policies > result_type;
    result_type operator()(unused_type, unused_type) const
    { return result_type(); }
};

template< class Policies, class T, class A0, class Modifiers >
struct make_primitive<
    terminal_ex<
        tag::stateful_tag< Policies, PhysBAM::Tag::RATIONAL_REAL<T> >,
        fusion::vector1<A0>
    >,
    Modifiers
>
{
    typedef PhysBAM::LITERAL_RATIONAL_REAL_PARSER< T, Policies > result_type;
    template< class Terminal >
    result_type operator()(const Terminal& terminal, unused_type) const
    { return result_type(static_cast<T>(fusion::at_c<0>(terminal.args))); }
};

} // namespace qi

} } // namespace spirit / namespace boost

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARSING_RATIONAL_REAL_PARSER_HPP
