//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CYLINDER
//##################################################################### 
#include <Core/Math_Tools/cube.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
using namespace PhysBAM;
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> CYLINDER<T>::
Normal(const TV& location,const int aggregate) const 
{
    assert(aggregate >= 1 && aggregate <= 3);
    if(aggregate == 1) return (location-plane1.x0).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Normalized(); // CYLINDER
    else if(aggregate == 2) return plane1.Normal();
    else return plane2.Normal();
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool CYLINDER<T>::
Inside(const TV& location,const T thickness_over_two) const 
{
    if(!plane1.Inside(location,thickness_over_two) || !plane2.Inside(location,thickness_over_two)) return false;
    return (location-plane1.x0).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() <= sqr(radius-thickness_over_two);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool CYLINDER<T>::
Lazy_Inside(const TV& location) const 
{
    if(!plane1.Lazy_Inside(location) || !plane2.Lazy_Inside(location)) return false;
    return (location-plane1.x0).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() <= sqr(radius);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool CYLINDER<T>::
Outside(const TV& location,const T thickness_over_two) const  
{
    if(plane1.Outside(location,thickness_over_two) || plane2.Outside(location,thickness_over_two)) return true;
    return (location-plane1.x0).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() >= sqr(radius+thickness_over_two);
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool CYLINDER<T>::
Lazy_Outside(const TV& location) const  
{
    if(plane1.Lazy_Outside(location) || plane2.Lazy_Outside(location)) return true;
    return (location-plane1.x0).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() >= sqr(radius);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool CYLINDER<T>::
Boundary(const TV& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> CYLINDER<T>::
Surface(const TV& location) const
{
    TV v=location-plane1.x0;T axial_distance=-TV::Dot_Product(v,plane1.normal);
    TV radial_direction=v+axial_distance*plane1.normal;T radial_distance=radial_direction.Normalize();
    T radial_distance_minus_radius=radial_distance-radius;
    if(radial_distance_minus_radius>0 || (radial_distance_minus_radius>-axial_distance && radial_distance_minus_radius>axial_distance-height)) // closest point is on infinite CYLINDER
        return plane1.x0-clamp(axial_distance,(T)0,height)*plane1.normal+radius*radial_direction;
    if(axial_distance < height-axial_distance) return location+axial_distance*plane1.normal; // closest point is on plane1
    else return location-(height-axial_distance)*plane1.normal; // closest point is on plane2
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T CYLINDER<T>::
Signed_Distance(const TV& location) const
{  
    TV v=location-plane1.x0;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane_distance=max(plane1_distance,-height-plane1_distance);
    T cylinder_distance=(v-plane1_distance*plane1.normal).Magnitude()-radius;
    return cylinder_distance>0 && plane_distance>0?sqrt(sqr(cylinder_distance)+sqr(plane_distance)):max(cylinder_distance,plane_distance);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> CYLINDER<T>::
Normal(const TV& location) const
{
    TV v=location-plane1.x0;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane2_distance=-height-plane1_distance;
    TV infinite_cylinder_normal=v-plane1_distance*plane1.normal;
    T cylinder_distance=infinite_cylinder_normal.Normalize()-radius;
    if(plane1_distance>=plane2_distance){
        if(cylinder_distance>0 && plane1_distance>0){
            T magnitude=sqrt(sqr(cylinder_distance)+sqr(plane1_distance));
            return cylinder_distance/magnitude*infinite_cylinder_normal+plane1_distance/magnitude*plane1.normal;}
        else if(cylinder_distance>plane1_distance) return infinite_cylinder_normal;
        return plane1.normal;}
    if(cylinder_distance>0 && plane2_distance>0){
        T magnitude=sqrt(sqr(cylinder_distance)+sqr(plane2_distance));
        return cylinder_distance/magnitude*infinite_cylinder_normal+plane2_distance/magnitude*plane2.normal;}
    else if(cylinder_distance>plane2_distance) return infinite_cylinder_normal;
    return plane2.normal;
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,3> CYLINDER<T>::
Hessian(const TV& X) const
{
    TV v=X-plane1.x0;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane_distance=max(plane1_distance,-height-plane1_distance);
    TV infinite_cylinder_normal=v-plane1_distance*plane1.normal;
    T norm_infinite_cylinder_normal=infinite_cylinder_normal.Normalize();
    T cylinder_distance=norm_infinite_cylinder_normal-radius;

    bool edge=(cylinder_distance>0 && plane_distance>0),face=plane_distance>cylinder_distance;
    if(!edge && face) return SYMMETRIC_MATRIX<T,3>();

    TV diff_cylinder_distance=infinite_cylinder_normal-plane1.normal.Dot(infinite_cylinder_normal)*plane1.normal;
    T norm2_v=v.Magnitude_Squared();
    T Hd=cube(norm_infinite_cylinder_normal);
    SYMMETRIC_MATRIX<T,3> H0=norm2_v-SYMMETRIC_MATRIX<T,3>::Outer_Product(v);
    SYMMETRIC_MATRIX<T,3> H1=norm2_v*SYMMETRIC_MATRIX<T,3>::Outer_Product(plane1.normal);
    SYMMETRIC_MATRIX<T,3> H2=plane1_distance*SYMMETRIC_MATRIX<T,3>::Symmetric_Outer_Product(plane1.normal,v);
    SYMMETRIC_MATRIX<T,3> hess_cylinder_distance=(H0-H1+H2-sqr(plane1_distance))/Hd;
    if(!edge) return hess_cylinder_distance;

    TV diff_plane_distance=(plane1_distance>-height-plane1_distance)?plane1.normal:-plane1.normal;
    T a=sqrt(sqr(cylinder_distance)+sqr(plane_distance));
    TV da=(cylinder_distance*diff_cylinder_distance+plane_distance*diff_plane_distance)/a;
    SYMMETRIC_MATRIX<T,3> A0=SYMMETRIC_MATRIX<T,3>::Outer_Product(diff_cylinder_distance);
    SYMMETRIC_MATRIX<T,3> A1=cylinder_distance*hess_cylinder_distance;
    SYMMETRIC_MATRIX<T,3> A2=SYMMETRIC_MATRIX<T,3>::Outer_Product(plane1.normal);
    SYMMETRIC_MATRIX<T,3> A3=SYMMETRIC_MATRIX<T,3>::Outer_Product(da);
    return (A0+A1+A2-A3)/a;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > CYLINDER<T>::
Bounding_Box() const
{
    TV C=(plane1.x0+plane2.x0)/2,u=plane1.normal,w=height/2*abs(u);
    return RANGE<TV>(C).Thickened(w+(Outer_Product(u,u)-1).Column_Magnitudes()*radius);
}
//#####################################################################
namespace PhysBAM{
template class CYLINDER<float>;
template class CYLINDER<double>;
}
