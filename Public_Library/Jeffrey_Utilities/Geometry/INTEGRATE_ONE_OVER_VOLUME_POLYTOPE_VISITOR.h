//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_ONE_OVER_VOLUME_POLYTOPE_VISITOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_ONE_OVER_VOLUME_POLYTOPE_VISITOR_HPP

#include <boost/mpl/int.hpp>
#include <boost/mpl/push_front.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>
#include <Jeffrey_Utilities/Math/GAUSSIAN_QUADRATURE_SIMPLEX_RULE.h>
#include <Jeffrey_Utilities/Math/MONOMIAL.h>
#include <Jeffrey_Utilities/Math/MONOMIAL_VECTOR.h>
#include <Jeffrey_Utilities/ZERO_MPL_VECTOR.h>

namespace PhysBAM
{

template< class T >
struct INTEGRATE_ONE_OVER_VOLUME_POLYTOPE_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INTEGRATE_ONE_OVER_VOLUME_POLYTOPE_VISITOR,
        (( /******/ int const, sign ))
        (( typename T&, result ))
    )
public:
    typedef void result_type;

    template< class T_VERTICES_X, int D >
    void operator()(const int sign_, const T_VERTICES_X& vertices_x, const SUB_CUBE_POLYTOPE<T,D>& polytope) const
    {
        typedef GAUSSIAN_QUADRATURE_SIMPLEX_RULE<T,D,0> QUADRATURE_RULE_TYPE;
        typedef typename boost::mpl::push_front<
            typename ZERO_MPL_VECTOR_C< D-1 >::type,
            boost::mpl::int_<1>
        >::type MULTI_POWER_TYPE;
        typedef MONOMIAL_VECTOR< MONOMIAL< MULTI_POWER_TYPE >, 1 > MONOMIAL_VECTOR_TYPE; 
        if(
            sign_ == sign &&
            (polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_BOUNDARY_ALIGNED ||
             polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_CUBE_ALIGNED ||
             polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED)
        )
            result += polytope.template Integrate_Flux< QUADRATURE_RULE_TYPE >(vertices_x, MONOMIAL_VECTOR_TYPE());
    }
};

template< class T >
inline INTEGRATE_ONE_OVER_VOLUME_POLYTOPE_VISITOR<T>
Make_Integrate_One_Over_Volume_Polytope_Visitor(const int sign, T& result)
{ return INTEGRATE_ONE_OVER_VOLUME_POLYTOPE_VISITOR<T>(sign, result); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_ONE_OVER_VOLUME_POLYTOPE_VISITOR_HPP
