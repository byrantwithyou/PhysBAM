//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR_HPP

#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/size.hpp>
#include <boost/type_traits/is_same.hpp>

#include <Jeffrey_Utilities/Geometry/INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Geometry/POLYTOPE_LOCAL_INTEGRATE_MONOMIAL_OVER_VOLUME.h>

namespace PhysBAM
{

template< class T_RESULT_OF_MULTI_POWER >
struct INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR_TRAITS
{
    template< class T_MULTI_POWER >
    struct MULTI_POWER_FILTER
    {
        BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_MULTI_POWER >));
        static const int DIMENSION = boost::mpl::size< T_MULTI_POWER >::value;
        static const int ORDER = boost::mpl::accumulate<
            T_MULTI_POWER,
            boost::mpl::int_<0>,
            boost::mpl::plus< boost::mpl::_1, boost::mpl::_2 >
        >::type::value;
        static const bool has_zero = !boost::is_same<
            typename boost::mpl::find_if<
                T_MULTI_POWER,
                boost::mpl::equal_to< boost::mpl::_1, boost::mpl::int_<0> >
            >::type,
            typename boost::mpl::end< T_MULTI_POWER >::type
        >::value;
        static const bool value = (ORDER <= DIMENSION) || has_zero;
        typedef MULTI_POWER_FILTER type;
    };
    typedef INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR<
        2,
        boost::mpl::quote1< MULTI_POWER_FILTER >,
        POLYTOPE_LOCAL_INTEGRATE_MONOMIAL_OVER_VOLUME,
        T_RESULT_OF_MULTI_POWER
    > INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR_;
};

template< class T_RESULT_OF_MULTI_POWER >
class INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR
    : public INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR_TRAITS< T_RESULT_OF_MULTI_POWER >::
          INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR_
{
    typedef INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR_TRAITS< T_RESULT_OF_MULTI_POWER >
        INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR_TRAITS_;
    typedef typename INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR_TRAITS_::
        INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR_ INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR_;
public:
    INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR(
        const int sign,
        const T_RESULT_OF_MULTI_POWER& result_of_multi_power)
        : INTEGRATE_MULTI_POWERS_POLYTOPE_VISITOR_(
              sign,
              POLYTOPE_LOCAL_INTEGRATE_MONOMIAL_OVER_VOLUME(),
              result_of_multi_power
          )
    { }
};

template< class T_RESULT_OF_MULTI_POWER >
inline INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR< T_RESULT_OF_MULTI_POWER >
Make_Integrate_Multi_Powers_Over_Volume_For_Poisson_Polytope_Visitor(
    const int sign,
    const T_RESULT_OF_MULTI_POWER& result_of_multi_power)
{
    return INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR< T_RESULT_OF_MULTI_POWER >(
        sign, result_of_multi_power
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR_HPP
