//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// TODO: Replace BASE_QUADRATURE_RULE_ORDER with a Boost.MPL
// metafunction that gives the quadrature rule.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_FLUX_FUNCTION_TIMES_MONOMIAL_OVER_BOUNDARY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_FLUX_FUNCTION_TIMES_MONOMIAL_OVER_BOUNDARY_HPP

#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/plus.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>
#include <Jeffrey_Utilities/Math/GAUSSIAN_QUADRATURE_SIMPLEX_RULE.h>
#include <Jeffrey_Utilities/Math/MONOMIAL.h>
#include <Jeffrey_Utilities/Math/MONOMIAL_TIMES_FLUX_FUNCTION.h>

namespace PhysBAM
{

template< int BASE_QUADRATURE_RULE_ORDER, class T_F >
struct POLYTOPE_LOCAL_INTEGRATE_FLUX_FUNCTION_TIMES_MONOMIAL_OVER_BOUNDARY
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POLYTOPE_LOCAL_INTEGRATE_FLUX_FUNCTION_TIMES_MONOMIAL_OVER_BOUNDARY,
        (( typename T_F const, f ))
    )
public:

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
        static const int QUADRATURE_RULE_ORDER =
            BASE_QUADRATURE_RULE_ORDER
          + boost::mpl::accumulate<
                T_MULTI_POWER,
                boost::mpl::int_<0>,
                boost::mpl::plus< boost::mpl::_1, boost::mpl::_2 >
            >::type::value;
        typedef GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, D, QUADRATURE_RULE_ORDER > QUADRATURE_RULE_TYPE;
        return polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_BOUNDARY_ALIGNED ||
               polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED ?
               polytope.template Integrate_Flux< QUADRATURE_RULE_TYPE >(
                   vertices_x,
                   Make_Monomial_Times_Flux_Function< MONOMIAL< T_MULTI_POWER > >(f)
               ) :
               static_cast<T>(0);
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_POLYTOPE_LOCAL_INTEGRATE_FLUX_FUNCTION_TIMES_MONOMIAL_OVER_BOUNDARY_HPP
