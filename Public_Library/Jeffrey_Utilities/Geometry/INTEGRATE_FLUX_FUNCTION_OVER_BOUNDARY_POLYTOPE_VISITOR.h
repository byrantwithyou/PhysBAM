//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>

namespace PhysBAM
{

template< class T, class T_QUADRATURE_RULE, class T_F >
struct INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR,
        (( /******/ int const, sign ))
        (( typename T_F const, f ))
        (( typename T&, result ))
    )
public:
    typedef void result_type;
    template< class T_VERTICES_X, int D >
    void operator()(
        const int sign_,
        const T_VERTICES_X& vertices_x,
        const SUB_CUBE_POLYTOPE<T,D>& polytope) const
    {
        if(
            sign_ == sign &&
            (polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_BOUNDARY_ALIGNED ||
             polytope.alignment == SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED)
        )
            result += polytope.template Integrate_Flux< T_QUADRATURE_RULE >(vertices_x, f);
    }
};

template< class T_QUADRATURE_RULE, class T_F, class T >
inline INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR< T, T_QUADRATURE_RULE, T_F >
Make_Integrate_Flux_Function_Over_Boundary_Polytope_Visitor(const int sign, const T_F& f, T& result)
{ return INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR< T, T_QUADRATURE_RULE, T_F >(sign, f, result); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR_HPP
