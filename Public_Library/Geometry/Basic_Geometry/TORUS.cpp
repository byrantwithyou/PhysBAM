//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Basic_Geometry/TORUS.h>
namespace PhysBAM{
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T TORUS<T>::
Signed_Distance(const TV& X) const
{
    TV x=X-center;
    T axial=Dot_Product(x,axis);
    T radial=(x-axial*axis).Magnitude();
    return sqrt(sqr(radial-outer_radius)+sqr(axial))-inner_radius;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> TORUS<T>::
Normal(const TV& X) const
{
    TV x=X-center;
    T axial=Dot_Product(x,axis);
    TV inplane=x-axial*axis;
    VECTOR<T,2> plane_normal(axial,inplane.Normalize()-outer_radius);
    plane_normal.Normalize();
    return plane_normal.x*axis+plane_normal.y*inplane;
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,3> TORUS<T>::
Hessian(const TV& X) const
{
    TV Y=X-center;
    T a=Y.Dot(axis);
    TV dw=Y-a*axis;
    T w=dw.Normalize();
    T z=w-outer_radius,p=hypot(z,a),e=z/p;
    return e/w*((T)1-Outer_Product(axis)-Outer_Product(dw))+Outer_Product(a/p*dw-e*axis)/p;
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> TORUS<T>::
Principal_Curvatures(const TV& X) const
{
    T inplane=(X-TV::Dot_Product(X,axis)*axis).Magnitude();
    return VECTOR<T,2>(-1/inner_radius,-(inplane-outer_radius)/(inplane*inner_radius));
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > TORUS<T>::
Bounding_Box() const
{
    T x=sqr(axis.x),y=sqr(axis.y),z=sqr(axis.z);
    TV corner_offset=sqrt(TV(y+z,x+z,x+y))*outer_radius+inner_radius;
    return RANGE<TV>(center-corner_offset,center+corner_offset);
}
template class TORUS<float>;
template class TORUS<double>;
}
