//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Geometry/Basic_Geometry/TORUS.h>
namespace PhysBAM{
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,3> TORUS<T>::
Hessian(const TV& X) const
{
    typedef DIFF_LAYOUT<T,TV::m> LAYOUT;
    auto x=Hess_From_Var<LAYOUT,0>(X)-center;
    auto axial=x.Dot(axis);
    auto ret=hypot((x-axial*axis).Magnitude()-outer_radius,axial)-inner_radius;
    SYMMETRIC_MATRIX<T,3> sm;
    Get<0,0>(sm,ret.ddx);
    return sm;
}
template class TORUS<float>;
template class TORUS<double>;
}
