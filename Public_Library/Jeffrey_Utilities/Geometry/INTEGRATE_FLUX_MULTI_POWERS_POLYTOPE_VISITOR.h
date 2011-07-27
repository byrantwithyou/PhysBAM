//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR_HPP

#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/empty_sequence.hpp>
#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/vector/vector10.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

template<
    int MAX_POWER,
    class T_MULTI_POWER_FILTER,
    class T_INTEGRATE_FLUX_MULTI_POWER_AXIS,
    class T_RESULT_OF_MULTI_POWER_AXIS
>
struct INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR,
        (( /******/ int const, sign ))
        (( typename T_INTEGRATE_FLUX_MULTI_POWER_AXIS const, integrate_flux_multi_power_axis ))
        (( typename T_RESULT_OF_MULTI_POWER_AXIS const, result_of_multi_power_axis ))
    )
public:
    typedef void result_type;
    template< class T_VERTICES_X, class T, int D >
    void operator()(
        const int sign_,
        const T_VERTICES_X& vertices_x,
        const SUB_CUBE_POLYTOPE<T,D>& polytope) const;
};

template<
    int MAX_POWER,
    class T_MULTI_POWER_FILTER,
    class T_INTEGRATE_FLUX_MULTI_POWER_AXIS,
    class T_RESULT_OF_MULTI_POWER_AXIS
>
inline INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR<
    MAX_POWER, T_MULTI_POWER_FILTER, T_INTEGRATE_FLUX_MULTI_POWER_AXIS, T_RESULT_OF_MULTI_POWER_AXIS
>
Make_Integrate_Flux_Multi_Powers_Polytope_Visitor(
    const int sign,
    const T_INTEGRATE_FLUX_MULTI_POWER_AXIS& integrate_flux_multi_power_axis,
    const T_RESULT_OF_MULTI_POWER_AXIS& result_of_multi_power_axis)
{
    return INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR<
        MAX_POWER, T_MULTI_POWER_FILTER, T_INTEGRATE_FLUX_MULTI_POWER_AXIS, T_RESULT_OF_MULTI_POWER_AXIS
    >(sign, integrate_flux_multi_power_axis, result_of_multi_power_axis);
}

//#####################################################################
//#####################################################################

namespace Detail_INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR
{

template< int D >
struct PAIR_WITH_ALL_AXIS
{
    template< class T >
    struct apply
        : boost::mpl::transform<
              boost::mpl::range_c< int, 1, D >,
              boost::mpl::pair< T, boost::mpl::_1 >,
              boost::mpl::back_inserter< boost::mpl::vector0<> >
          >
    { };
};

template<
    class T_INTEGRATE_FLUX_MULTI_POWER_AXIS, class T_RESULT_OF_MULTI_POWER_AXIS,
    class T_VERTICES_X, class T_POLYTOPE
>
struct INTEGRATE_MULTI_POWER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INTEGRATE_MULTI_POWER,
        (( typename T_INTEGRATE_FLUX_MULTI_POWER_AXIS const, integrate_flux_multi_power_axis ))
        (( typename T_RESULT_OF_MULTI_POWER_AXIS const, result_of_multi_power_axis ))
        (( typename T_VERTICES_X const &, vertices_x ))
        (( typename T_POLYTOPE const &, polytope ))
    )
public:
    typedef void result_type;
    template< class T_MULTI_POWER, class T_AXIS >
    void operator()(boost::mpl::pair< T_MULTI_POWER, T_AXIS >) const
    {
        result_of_multi_power_axis(T_MULTI_POWER(), T_AXIS()) +=
            integrate_flux_multi_power_axis(polytope, vertices_x, T_MULTI_POWER(), T_AXIS());
    }
};

} // namespace Detail_INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR

template<
    int MAX_POWER,
    class T_MULTI_POWER_FILTER,
    class T_INTEGRATE_FLUX_MULTI_POWER_AXIS,
    class T_RESULT_OF_MULTI_POWER_AXIS
>
template< class T_VERTICES_X, class T, int D >
inline void
INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR<
    MAX_POWER, T_MULTI_POWER_FILTER, T_INTEGRATE_FLUX_MULTI_POWER_AXIS, T_RESULT_OF_MULTI_POWER_AXIS
>::
operator()(
    const int sign_,
    const T_VERTICES_X& vertices_x,
    const SUB_CUBE_POLYTOPE<T,D>& polytope) const
{
    if(sign_ != sign)
        return;
    static const int MAX_LINEAR_INDEX = STATIC_POW_C< 1 + MAX_POWER, D >::value;
    typedef STATIC_MULTI_INDEX_CUBE< D, 0, MAX_POWER > STATIC_MULTI_INDEX_CUBE_TYPE;
    typedef boost::mpl::filter_view<
        boost::mpl::transform_view<
            boost::mpl::range_c< int, 1, 1 + MAX_LINEAR_INDEX >,
            boost::mpl::quote1< STATIC_MULTI_INDEX_CUBE_TYPE::template MULTI_INDEX >
        >,
        T_MULTI_POWER_FILTER
    > MULTI_POWER_SEQUENCE_TYPE;
    typedef Detail_INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR::PAIR_WITH_ALL_AXIS<D> PAIR_WITH_ALL_AXIS_;
    typedef typename boost::mpl::fold<
        MULTI_POWER_SEQUENCE_TYPE,
        boost::mpl::empty_sequence,
        boost::mpl::joint_view<
            boost::mpl::_1,
            typename PAIR_WITH_ALL_AXIS_::template apply< boost::mpl::_2 >
        >
    >::type MULTI_POWER_AXIS_SEQUENCE_TYPE;
    typedef Detail_INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR::INTEGRATE_MULTI_POWER<
        T_INTEGRATE_FLUX_MULTI_POWER_AXIS, T_RESULT_OF_MULTI_POWER_AXIS,
        T_VERTICES_X, SUB_CUBE_POLYTOPE<T,D>
    > INTEGRATE_MULTI_POWER_AXIS_;
    boost::mpl::for_each< MULTI_POWER_AXIS_SEQUENCE_TYPE >(INTEGRATE_MULTI_POWER_AXIS_(
        integrate_flux_multi_power_axis, result_of_multi_power_axis,
        vertices_x, polytope
    ));
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_FLUX_MULTI_POWERS_POLYTOPE_VISITOR_HPP
