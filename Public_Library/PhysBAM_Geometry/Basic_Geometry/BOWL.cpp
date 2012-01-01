//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOWL
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOWL.h>
namespace PhysBAM{
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > BOWL<T>::
Bounding_Box() const
{
    return CYLINDER<T>(TV(0,0,0),TV(0,depth+thickness,0),outer_radius).Bounding_Box();
}
//#####################################################################
// Function Compute_Helper
//#####################################################################
template<class T> void BOWL<T>::
Compute_Helper(const TV& X,HELPER& h) const
{
    h.radial=X-TV(0,X.y,0);
    h.radius=h.radial.Normalize();
    h.dX=VECTOR<T,2>(h.radius-hole_radius,X.y);
    if( (h.dX.x>=0 && h.dX.y>=0) || (h.dX.x<=0 && h.dX.y<=0) )
    {
        h.dr=h.dX.Magnitude()*((h.dX.x>=0 && h.dX.y>=0)?(T)1:(T)(-1))-(depth+height)*0.5;
        h.signed_distance=abs(h.dr)-thickness*0.5;
        if (h.dr>0)
        {
            h.c1=-1/height;
            h.c2=-h.dX.x/(height*h.radius);
        }
        else
        {
            h.c1=abs(h.dX.x)/(depth*h.radius);
            h.c2=1/depth;
        }
    }
    else
    {
        T max_value = max(h.dX.x,h.dX.y);
        T min_value = min(h.dX.x,h.dX.y);
        if (max_value<depth) h.signed_distance=VECTOR<T,2>(min_value,max_value-depth).Magnitude();
        else if (max_value>height) h.signed_distance=VECTOR<T,2>(min_value,max_value-height).Magnitude();
        else h.signed_distance=-min_value;

        if (h.dX.x>h.dX.y)
        {
            h.c1=0;
            h.c2=0;
        }
        else
        {
            h.c1=-1/hole_radius;
            h.c2=0;
        }
    }
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T BOWL<T>::
Signed_Distance(const TV& X) const
{
    HELPER h;
    Compute_Helper(X,h);
    return h.signed_distance;
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> BOWL<T>::
Surface(const TV& X,const HELPER& h) const
{
    VECTOR<T,2> slice_surface;
    if( (h.dX.x>=0 && h.dX.y>=0) || (h.dX.x<=0 && h.dX.y<=0) )
    {
        VECTOR<T,2> dX_Normalized=h.dX;
        dX_Normalized.Normalize();
        slice_surface=dX_Normalized*((h.dr>=0)?height:depth);
    }
    else
    {
        if (h.dX.x>h.dX.y)
        {
            if (h.dX.x<depth) slice_surface=VECTOR<T,2>(depth,0);
            else if (h.dX.x>height) slice_surface=VECTOR<T,2>(height,0);
            else slice_surface=VECTOR<T,2>(h.dX.x,0);
        }
        else
        {
            if (h.dX.y<depth) slice_surface=VECTOR<T,2>(0,depth);
            else if (h.dX.y>height) slice_surface=VECTOR<T,2>(0,height);
            else slice_surface=VECTOR<T,2>(0,h.dX.y);
        }
    }
    return h.radial*slice_surface.x+TV(0,slice_surface.y,0);
}
//#####################################################################
// Function Suface
//#####################################################################
template<class T> VECTOR<T,3> BOWL<T>::
Surface(const TV& X) const 
{
    HELPER h;
    Compute_Helper(X,h);
    return Surface(X,h);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> BOWL<T>::
Normal(const TV& X,const HELPER& h) const
{
    VECTOR<T,2> slice_normal;
    if( (h.dX.x>=0 && h.dX.y>=0) || (h.dX.x<=0 && h.dX.y<=0) )
    {
        slice_normal=(h.dX*((h.dr>=0||(h.dX.x<=0 && h.dX.y<=0))?(T)1:(T)(-1)));
        slice_normal.Normalize();
    }
    else
    {
        slice_normal=((h.dX.x>h.dX.y)?VECTOR<T,2>(0,-1):VECTOR<T,2>(-1,0));
    }
    return h.radial*slice_normal.x+TV(0,slice_normal.y,0);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> BOWL<T>::
Normal(const TV& X) const 
{
    HELPER h;
    Compute_Helper(X,h);
    return Normal(X,h);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> BOWL<T>::
Normal(const TV& X,const int aggregate) const 
{
    return Normal(X);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> BOWL<T>::
Principal_Curvatures(const TV& X) const
{
    HELPER h;
    Compute_Helper(X,h);
    return VECTOR<T,2>(h.c1,h.c2);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool BOWL<T>::
Lazy_Inside(const TV& X) const 
{
    return Signed_Distance(X)<0;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool BOWL<T>::
Lazy_Outside(const TV& X) const  
{
    return !Lazy_Inside(X);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool BOWL<T>::
Inside(const TV& X,const T thickness_over_two) const 
{
    return Signed_Distance(X)<=-thickness_over_two;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool BOWL<T>::
Outside(const TV& X,const T thickness_over_two) const  
{
    return Signed_Distance(X)>=thickness_over_two;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool BOWL<T>::
Boundary(const TV& X,const T thickness_over_two) const
{
    T sd=Signed_Distance(X);
    return (sd<thickness_over_two && sd>-thickness_over_two);
}
//#####################################################################
// Function Name
//#####################################################################
template<class T> std::string BOWL<T>::
Name()
{
    return "BOWL<T>";
}
//#####################################################################
template class BOWL<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOWL<double>;
#endif
}
