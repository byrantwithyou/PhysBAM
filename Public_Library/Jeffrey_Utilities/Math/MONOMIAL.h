//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_HPP

#include <boost/math/special_functions/pow.hpp>
#include <boost/mpl/advance.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/iterator_range.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/single_view.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T_MULTI_POWER, int INV_SCALING = 1 >
struct MONOMIAL
{
    BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_MULTI_POWER >));
    static const int DIMENSION = boost::mpl::size< T_MULTI_POWER >::value;

    typedef T_MULTI_POWER MULTI_POWER_TYPE;

    template< int INTEGRATION_AXIS > class INTEGRATE_C_RESULT;
    template< int INTEGRATION_AXIS >
    static typename INTEGRATE_C_RESULT< INTEGRATION_AXIS >::type
    Integrate_C();

    template<class> struct result;
    template< class T_THIS, class T_VECTOR >
    struct result< T_THIS ( T_VECTOR ) >
    { typedef typename boost::remove_reference< T_VECTOR >::type::value_type type; };

    template< class T >
    T operator()(const VECTOR< T, DIMENSION >& x) const;
};

//#####################################################################
//#####################################################################

template< class T_MULTI_POWER, int INV_SCALING >
template< int INTEGRATION_AXIS >
class MONOMIAL< T_MULTI_POWER, INV_SCALING >::
INTEGRATE_C_RESULT
{
    static const int INV_INTEGRATION_SCALING =
        1 + boost::mpl::at_c< T_MULTI_POWER, INTEGRATION_AXIS - 1 >::type::value;
    typedef typename boost::mpl::copy<
        boost::mpl::joint_view<
            boost::mpl::joint_view<
                boost::mpl::iterator_range<
                    typename boost::mpl::begin< T_MULTI_POWER >::type,
                    typename boost::mpl::advance<
                        typename boost::mpl::begin< T_MULTI_POWER >::type,
                        boost::mpl::int_< INTEGRATION_AXIS - 1 >
                    >::type
                >,
                boost::mpl::single_view< boost::mpl::int_< INV_INTEGRATION_SCALING > >
            >,
            boost::mpl::iterator_range<
                typename boost::mpl::advance<
                    typename boost::mpl::begin< T_MULTI_POWER >::type,
                    boost::mpl::int_< INTEGRATION_AXIS >
                >::type,
                typename boost::mpl::end< T_MULTI_POWER >::type
            >
        >,
        boost::mpl::back_inserter< boost::mpl::vector0<> >
    >::type INTEGRATED_MULTI_POWER_TYPE;
public:
    typedef MONOMIAL<
        INTEGRATED_MULTI_POWER_TYPE,
        INV_SCALING * INV_INTEGRATION_SCALING
    > type;
};

template< class T_MULTI_POWER, int INV_SCALING >
template< int INTEGRATION_AXIS >
inline typename MONOMIAL< T_MULTI_POWER, INV_SCALING >::template
    INTEGRATE_C_RESULT< INTEGRATION_AXIS >::type
MONOMIAL< T_MULTI_POWER, INV_SCALING >::
Integrate_C()
{ return typename INTEGRATE_C_RESULT< INTEGRATION_AXIS >::type(); }

namespace Detail_MONOMIAL
{

template< class T_MULTI_POWER, int N = boost::mpl::size< T_MULTI_POWER >::value >
struct MONOMIAL_APPLY_DISPATCH
{
    template< class T, int D >
    static T Apply(const VECTOR<T,D>& x)
    {
        static const int POWER = boost::mpl::at_c< T_MULTI_POWER, N-1 >::type::value;
        return MONOMIAL_APPLY_DISPATCH< T_MULTI_POWER, N-1 >::Apply(x) * boost::math::pow< POWER >(x[N]);
    }
};

template< class T_MULTI_POWER >
struct MONOMIAL_APPLY_DISPATCH< T_MULTI_POWER, 0 >
{
    template< class T, int D >
    static T Apply(const VECTOR<T,D>&)
    { return static_cast<T>(1); }
};

template< class T_MULTI_POWER >
struct MONOMIAL_APPLY_DISPATCH< T_MULTI_POWER, 1 >
{
    template< class T, int D >
    static T Apply(const VECTOR<T,D>& x)
    {
        static const int POWER = boost::mpl::at_c< T_MULTI_POWER, 0 >::type::value;
        return boost::math::pow< POWER >(x[1]);
    }
};

} // namespace Detail_MONOMIAL

template< class T_MULTI_POWER, int INV_SCALING >
template< class T >
inline T
MONOMIAL< T_MULTI_POWER, INV_SCALING >::
operator()(const VECTOR< T, DIMENSION >& x) const
{ return Detail_MONOMIAL::MONOMIAL_APPLY_DISPATCH< T_MULTI_POWER >::Apply(x) / INV_SCALING; }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_HPP
