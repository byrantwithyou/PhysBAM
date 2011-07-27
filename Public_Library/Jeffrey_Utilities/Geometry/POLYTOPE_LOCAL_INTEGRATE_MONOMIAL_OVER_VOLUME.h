//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_MONOMIAL_OVER_VOLUME_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_MONOMIAL_OVER_VOLUME_HPP

#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/min_element.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>
#include <Jeffrey_Utilities/Math/GAUSSIAN_QUADRATURE_SIMPLEX_RULE.h>
#include <Jeffrey_Utilities/Math/MONOMIAL.h>

namespace PhysBAM
{

struct POLYTOPE_LOCAL_INTEGRATE_MONOMIAL_OVER_VOLUME
{
    template<class> struct result;
    template< class T_THIS, class T_POLYTOPE, class T_VERTICES_X, class T_MULTI_POWER >
    struct result< T_THIS ( T_POLYTOPE, T_VERTICES_X, T_MULTI_POWER ) >
    { typedef typename boost::remove_reference< T_POLYTOPE >::type::SCALAR_TYPE type; };

    template< class T, int D, class T_VERTICES_X, class T_MULTI_POWER >
    T operator()(
        const SUB_CUBE_POLYTOPE<T,D>& polytope,
        const T_VERTICES_X& vertices_x,
        T_MULTI_POWER) const
    {
        BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_MULTI_POWER >));
        static const int INTEGRATION_AXIS = 1 + boost::mpl::distance<
            typename boost::mpl::begin< T_MULTI_POWER >::type,
            typename boost::mpl::min_element< T_MULTI_POWER >::type
        >::value;
        typedef typename MONOMIAL< T_MULTI_POWER >::template
            INTEGRATE_C_RESULT< INTEGRATION_AXIS >::type INTEGRATED_MONOMIAL_TYPE;
        static const int QUADRATURE_RULE_ORDER = boost::mpl::accumulate<
            T_MULTI_POWER,
            boost::mpl::int_<1>,
            boost::mpl::plus< boost::mpl::_1, boost::mpl::_2 >
        >::type::value;
        typedef GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, D, QUADRATURE_RULE_ORDER > QUADRATURE_RULE_TYPE;
        return polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_BOUNDARY_ALIGNED ||
               polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_CUBE_ALIGNED ||
               polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED ?
               polytope.normal[INTEGRATION_AXIS] *
               polytope.template Integrate< QUADRATURE_RULE_TYPE >(vertices_x, INTEGRATED_MONOMIAL_TYPE()) :
               static_cast<T>(0);
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_MONOMIAL_OVER_VOLUME_HPP
