//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_FLUX_MONOMIAL_VECTOR_OVER_BOUNDARY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_FLUX_MONOMIAL_VECTOR_OVER_BOUNDARY_HPP

#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/plus.hpp>

#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>
#include <Jeffrey_Utilities/Math/GAUSSIAN_QUADRATURE_SIMPLEX_RULE.h>
#include <Jeffrey_Utilities/Math/MONOMIAL_VECTOR.h>

namespace PhysBAM
{

struct POLYTOPE_LOCAL_INTEGRATE_FLUX_MONOMIAL_VECTOR_OVER_BOUNDARY
{
    template<class> struct result;
    template< class T_THIS, class T_POLYTOPE, class T_VERTICES_X, class T_MULTI_POWER >
    struct result< T_THIS ( T_POLYTOPE, T_VERTICES_X, T_MULTI_POWER ) >
    { typedef typename boost::remove_reference< T_POLYTOPE >::type::SCALAR_TYPE type; };

    template< class T, int D, class T_VERTICES_X, class T_MULTI_POWER, class T_INT >
    T operator()(
        const SUB_CUBE_POLYTOPE<T,D>& polytope,
        const T_VERTICES_X& vertices_x,
        T_MULTI_POWER, T_INT) const
    {
        BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_MULTI_POWER >));
        static const int QUADRATURE_RULE_ORDER = boost::mpl::accumulate<
            T_MULTI_POWER,
            boost::mpl::int_<0>,
            boost::mpl::plus< boost::mpl::_1, boost::mpl::_2 >
        >::type::value;
        typedef GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, D, QUADRATURE_RULE_ORDER > QUADRATURE_RULE_TYPE;
        return polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_BOUNDARY_ALIGNED ||
               polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED ?
               polytope.template Integrate_Flux< QUADRATURE_RULE_TYPE >(
                   vertices_x,
                   MONOMIAL_VECTOR< T_MULTI_POWER, T_INT::value >()
               ) :
               static_cast<T>(0);
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_FLUX_MONOMIAL_VECTOR_OVER_BOUNDARY_HPP
