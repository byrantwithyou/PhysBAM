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
    h.radius=h.radial.Magnitude();
    h.dX=VECTOR<T,2>(h.radius-hole_radius,X.y);
    if( (h.dX.x>=0 && h.dX.y>=0) || (h.dX.x<=0 && h.dX.y<=0) )
    {
        h.dr=h.dX.Magnitude()*((h.dX.x>=0 && h.dX.y>=0)?(T)1:(T)(-1))-(depth+height)*0.5;
        h.signed_distance=abs(h.dr)-thickness*0.5;
    }
    else
    {
        T max_value = max(h.dX.x,h.dX.y);
        T min_value = min(h.dX.x,h.dX.y);
        if (max_value<depth) h.signed_distance=VECTOR<T,2>(min_value,max_value-depth).Magnitude();
        else if (max_value>height) h.signed_distance=VECTOR<T,2>(min_value,max_value-height).Magnitude();
        else h.signed_distance=-min_value;
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
// TODO
    return TV();
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
//TODO
    return TV();
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
    return VECTOR<T,2>();//VECTOR<T,1>(h.ui?1/s:-1/s);
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
