//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/AUTO_HESS.h>
#include <Geometry/Basic_Geometry/TORUS.h>
namespace PhysBAM{
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,3> TORUS<T>::
Hessian(const TV& X) const
{
    AUTO_HESS<TV,TV> x=AUTO_HESS<TV,TV>::From_Var(X)-center;
    AUTO_HESS<T,TV> axial=x.Dot(axis),radial=(x-axial*axis).Magnitude();
    AUTO_HESS<T,TV> ret=hypot(radial-outer_radius,axial)-inner_radius;
    return ret.ddx;
}
template class TORUS<float>;
template class TORUS<double>;
}
