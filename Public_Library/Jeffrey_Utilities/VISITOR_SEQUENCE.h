//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VISITOR_SEQUENCE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VISITOR_SEQUENCE_HPP

#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/container/vector/vector10.hpp>
#include <boost/fusion/mpl.hpp>
#include <boost/fusion/support/is_sequence.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <boost/utility/enable_if.hpp>

namespace PhysBAM
{

template< class T_SEQUENCE >
struct VISITOR_SEQUENCE
{
    BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_SEQUENCE >));
    typedef T_SEQUENCE SEQUENCE_TYPE;

    template< class T_SEQUENCE2 >
    explicit VISITOR_SEQUENCE(const T_SEQUENCE2& visitors);

    typedef void result_type;

    void operator()() const;
    template< class T1 >
    void operator()(T1& x1) const;
    template< class T1 >
    void operator()(const T1& x1) const;
    template< class T1, class T2 >
    void operator()(T1& x1, T2& x2) const;
    template< class T1, class T2 >
    void operator()(T1& x1, const T2& x2) const;
    template< class T1, class T2 >
    void operator()(const T1& x1, T2& x2) const;
    template< class T1, class T2 >
    void operator()(const T1& x1, const T2& x2) const;
    template< class T1, class T2, class T3 >
    void operator()(const T1& x1, const T2& x2, const T3& x3) const;

private:
    typename boost::fusion::result_of::as_vector< T_SEQUENCE >::type m_visitors;
};

template< class T_SEQUENCE >
struct MAKE_VISITOR_SEQUENCE_RESULT
{
    BOOST_MPL_ASSERT((boost::fusion::traits::is_sequence< T_SEQUENCE >));
    BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_SEQUENCE >));
    typedef VISITOR_SEQUENCE<
        typename boost::mpl::copy<
            T_SEQUENCE,
            boost::mpl::back_inserter< boost::mpl::vector0<> >
        >::type
    > type;
};

template< class T_SEQUENCE >
struct MAKE_VISITOR_SEQUENCE_RESULT< T_SEQUENCE& >
    : MAKE_VISITOR_SEQUENCE_RESULT< T_SEQUENCE >
{ };

template< class T_SEQUENCE >
struct MAKE_VISITOR_SEQUENCE_RESULT< const T_SEQUENCE >
    : MAKE_VISITOR_SEQUENCE_RESULT< T_SEQUENCE >
{ };

template< class T_SEQUENCE >
inline typename boost::lazy_enable_if<
    boost::fusion::traits::is_sequence< T_SEQUENCE >,
    MAKE_VISITOR_SEQUENCE_RESULT< T_SEQUENCE >
>::type
Make_Visitor_Sequence(const T_SEQUENCE& visitors)
{ return typename MAKE_VISITOR_SEQUENCE_RESULT< T_SEQUENCE >::type(visitors); }

template< class T_VISITOR1, class T_VISITOR2 >
inline VISITOR_SEQUENCE< boost::mpl::vector2< T_VISITOR1, T_VISITOR2 > >
Make_Visitor_Sequence(const T_VISITOR1& visitor1, const T_VISITOR2& visitor2)
{
    return VISITOR_SEQUENCE< boost::mpl::vector2< T_VISITOR1, T_VISITOR2 > >(
        boost::fusion::make_vector(visitor1, visitor2)
    );
}

template< class T_VISITOR1, class T_VISITOR2, class T_VISITOR3 >
inline VISITOR_SEQUENCE< boost::mpl::vector3< T_VISITOR1, T_VISITOR2, T_VISITOR3 > >
Make_Visitor_Sequence(const T_VISITOR1& visitor1, const T_VISITOR2& visitor2, const T_VISITOR3& visitor3)
{
    return VISITOR_SEQUENCE< boost::mpl::vector3< T_VISITOR1, T_VISITOR2, T_VISITOR3 > >(
        boost::fusion::make_vector(visitor1, visitor2, visitor3)
    );
}

//#####################################################################
//#####################################################################

template< class T_SEQUENCE >
template< class T_SEQUENCE2 >
inline
VISITOR_SEQUENCE< T_SEQUENCE >::
VISITOR_SEQUENCE(const T_SEQUENCE2& visitors)
    : m_visitors(visitors)
{ }

namespace Detail_VISITOR_SEQUENCE
{

template< class T1 = void, class T2 = void, class T3 = void >
struct VISIT_FUNCTION
{
    VISIT_FUNCTION(T1 x1, T2 x2, T3 x3) : m_x1(x1), m_x2(x2), m_x3(x3) { }
    typedef void result_type;
    template< class T_VISITOR >
    void operator()(const T_VISITOR& visitor) const
    { visitor(m_x1, m_x2, m_x3); }
private:
    T1 m_x1; T2 m_x2; T3 m_x3;
};

template<>
struct VISIT_FUNCTION< void, void, void >
{
    typedef void result_type;
    template< class T_VISITOR >
    void operator()(const T_VISITOR& visitor) const
    { visitor(); }
};

template< class T1 >
struct VISIT_FUNCTION< T1, void, void >
{
    explicit VISIT_FUNCTION(T1 x1) : m_x1(x1) { }
    typedef void result_type;
    template< class T_VISITOR >
    void operator()(const T_VISITOR& visitor) const
    { visitor(m_x1); }
private:
    T1 m_x1;
};

template< class T1, class T2 >
struct VISIT_FUNCTION< T1, T2, void >
{
    VISIT_FUNCTION(T1 x1, T2 x2) : m_x1(x1), m_x2(x2) { }
    typedef void result_type;
    template< class T_VISITOR >
    void operator()(const T_VISITOR& visitor) const
    { visitor(m_x1, m_x2); }
private:
    T1 m_x1; T2 m_x2;
};

} // namespace Detail_VISITOR_SEQUENCE

template< class T_SEQUENCE >
inline void
VISITOR_SEQUENCE< T_SEQUENCE >::
operator()() const
{
    typedef Detail_VISITOR_SEQUENCE::VISIT_FUNCTION<> VISIT_FUNCTION_;
    boost::fusion::for_each(m_visitors, VISIT_FUNCTION_());
}

template< class T_SEQUENCE >
template< class T1 >
inline void
VISITOR_SEQUENCE< T_SEQUENCE >::
operator()(T1& x1) const
{
    typedef Detail_VISITOR_SEQUENCE::VISIT_FUNCTION< T1& > VISIT_FUNCTION_;
    boost::fusion::for_each(m_visitors, VISIT_FUNCTION_(x1));
}

template< class T_SEQUENCE >
template< class T1 >
inline void
VISITOR_SEQUENCE< T_SEQUENCE >::
operator()(const T1& x1) const
{
    typedef Detail_VISITOR_SEQUENCE::VISIT_FUNCTION< const T1& > VISIT_FUNCTION_;
    boost::fusion::for_each(m_visitors, VISIT_FUNCTION_(x1));
}

template< class T_SEQUENCE >
template< class T1, class T2 >
inline void
VISITOR_SEQUENCE< T_SEQUENCE >::
operator()(T1& x1, T2& x2) const
{
    typedef Detail_VISITOR_SEQUENCE::VISIT_FUNCTION< T1&, T2& > VISIT_FUNCTION_;
    boost::fusion::for_each(m_visitors, VISIT_FUNCTION_(x1, x2));
}

template< class T_SEQUENCE >
template< class T1, class T2 >
inline void
VISITOR_SEQUENCE< T_SEQUENCE >::
operator()(T1& x1, const T2& x2) const
{
    typedef Detail_VISITOR_SEQUENCE::VISIT_FUNCTION< T1&, const T2& > VISIT_FUNCTION_;
    boost::fusion::for_each(m_visitors, VISIT_FUNCTION_(x1, x2));
}

template< class T_SEQUENCE >
template< class T1, class T2 >
inline void
VISITOR_SEQUENCE< T_SEQUENCE >::
operator()(const T1& x1, T2& x2) const
{
    typedef Detail_VISITOR_SEQUENCE::VISIT_FUNCTION< const T1&, T2& > VISIT_FUNCTION_;
    boost::fusion::for_each(m_visitors, VISIT_FUNCTION_(x1, x2));
}

template< class T_SEQUENCE >
template< class T1, class T2 >
inline void
VISITOR_SEQUENCE< T_SEQUENCE >::
operator()(const T1& x1, const T2& x2) const
{
    typedef Detail_VISITOR_SEQUENCE::VISIT_FUNCTION< const T1&, const T2& > VISIT_FUNCTION_;
    boost::fusion::for_each(m_visitors, VISIT_FUNCTION_(x1, x2));
}

template< class T_SEQUENCE >
template< class T1, class T2, class T3 >
inline void
VISITOR_SEQUENCE< T_SEQUENCE >::
operator()(const T1& x1, const T2& x2, const T3& x3) const
{
    typedef Detail_VISITOR_SEQUENCE::VISIT_FUNCTION< const T1&, const T2&, const T3& > VISIT_FUNCTION_;
    boost::fusion::for_each(m_visitors, VISIT_FUNCTION_(x1, x2, x3));
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VISITOR_SEQUENCE_HPP
