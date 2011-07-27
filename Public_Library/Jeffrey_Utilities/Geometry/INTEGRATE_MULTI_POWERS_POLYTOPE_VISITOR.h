//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR_HPP

#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/transform_view.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

template<
    int MAX_POWER,
    class T_MULTI_POWER_FILTER,
    class T_INTEGRATE_MULTI_POWER,
    class T_RESULT_OF_MULTI_POWER
>
struct INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR,
        (( /******/ int const, sign ))
        (( typename T_INTEGRATE_MULTI_POWER const, integrate_multi_power ))
        (( typename T_RESULT_OF_MULTI_POWER const, result_of_multi_power ))
    )
public:
    typedef void result_type;
    template< class T_VERTICES_X, class T, int D >
    void operator()(
        const int sign_,
        const T_VERTICES_X& vertices_x,
        const SUB_CUBE_POLYTOPE<T,D>& polytope) const;
};

template< int MAX_POWER, class T_MULTI_POWER_FILTER, class T_INTEGRATE_MULTI_POWER, class T_RESULT_OF_MULTI_POWER >
inline INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR<
    MAX_POWER, T_MULTI_POWER_FILTER, T_INTEGRATE_MULTI_POWER, T_RESULT_OF_MULTI_POWER
>
Make_Integrate_Monomials_Polytope_Visitor(
    const int sign,
    const T_INTEGRATE_MULTI_POWER& integrate_multi_power,
    const T_RESULT_OF_MULTI_POWER& result_of_multi_power)
{
    return INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR<
        MAX_POWER, T_MULTI_POWER_FILTER, T_INTEGRATE_MULTI_POWER, T_RESULT_OF_MULTI_POWER
    >(sign, integrate_multi_power, result_of_multi_power);
}

//#####################################################################
//#####################################################################

namespace Detail_INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR
{

template<
    class T_INTEGRATE_MULTI_POWER, class T_RESULT_OF_MULTI_POWER,
    class T_VERTICES_X, class T_POLYTOPE
>
struct INTEGRATE_MULTI_POWER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INTEGRATE_MULTI_POWER,
        (( typename T_INTEGRATE_MULTI_POWER const, integrate_multi_power ))
        (( typename T_RESULT_OF_MULTI_POWER const, result_of_multi_power ))
        (( typename T_VERTICES_X const &, vertices_x ))
        (( typename T_POLYTOPE const &, polytope ))
    )
public:
    typedef void result_type;
    template< class T_MULTI_POWER >
    void operator()(T_MULTI_POWER) const
    {
        result_of_multi_power(T_MULTI_POWER()) +=
            integrate_multi_power(polytope, vertices_x, T_MULTI_POWER());
    }
};

} // namespace Detail_INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR

template< int MAX_POWER, class T_MULTI_POWER_FILTER, class T_INTEGRATE_MULTI_POWER, class T_RESULT_OF_MULTI_POWER >
template< class T_VERTICES_X, class T, int D >
inline void
INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR<
    MAX_POWER, T_MULTI_POWER_FILTER, T_INTEGRATE_MULTI_POWER, T_RESULT_OF_MULTI_POWER
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
    typedef Detail_INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR::INTEGRATE_MULTI_POWER<
        T_INTEGRATE_MULTI_POWER, T_RESULT_OF_MULTI_POWER,
        T_VERTICES_X, SUB_CUBE_POLYTOPE<T,D>
    > INTEGRATE_MULTI_POWER_;
    boost::mpl::for_each< MULTI_POWER_SEQUENCE_TYPE >(INTEGRATE_MULTI_POWER_(
        integrate_multi_power, result_of_multi_power,
        vertices_x, polytope
    ));
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR_HPP
