//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Geoffrey Irving, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONE
//##################################################################### 
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/CONE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

namespace{
template<class T,class TV>
struct HELPER
{
    T p,cr,ch,x,y,yy;
    TV dy,d;
    int cs;

    HELPER(const CONE<T>& cone,const TV& X)
    {
        d=cone.dir;
        TV u=X-cone.base;
        x=u.Dot(d);
        TV w=u-d*x;
        dy=w;
        yy=y=dy.Normalize();
        if(x<=0 && y<=cone.radius){p=-x;cs=0;return;}
        T yc=y-cone.radius;
        if(yc*cone.radius>cone.height*x){y=yc;p=hypot(x,y);cs=1;return;}
        T xc=x-cone.height;
        if(xc*cone.height>cone.radius*y){x=xc;p=hypot(x,y);cs=1;return;}
        T h=hypot(cone.radius,cone.height);
        cr=cone.radius/h;
        ch=cone.height/h;
        T z=xc*cr+y*ch;
        if(-x>=z){p=-x;cs=0;return;}
        p=z;
        cs=2;
    }

    TV N() const
    {
        if(cs==0) return -d;
        if(cs==1) return x/p*d+y/p*dy;
        return d*cr+dy*ch;
    }

    SYMMETRIC_MATRIX<T,3> H() const
    {
        if(cs==0) return {};
        SYMMETRIC_MATRIX<T,3> ddy=((T)1-Outer_Product(d)-Outer_Product(dy))/yy;
        if(cs==1) return y/p*ddy+Outer_Product(x/p*dy-y/p*d)/p;
        return ddy*ch;
    }
};
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > CONE<T>::
Bounding_Box() const
{
    TV w=sqrt((T)1-sqr(dir))*radius;
    RANGE<TV> range(-w,w);
    range.Enlarge_To_Include_Point(dir*height);
    return range+base;
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T CONE<T>::
Signed_Distance(const TV& X) const
{
    HELPER<T,TV> h(*this,X);
    return h.p;
}
//#####################################################################
// Function Suface
//#####################################################################
template<class T> VECTOR<T,3> CONE<T>::
Surface(const TV& X) const 
{
    HELPER<T,TV> h(*this,X);
    return X-h.p*h.N();
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> CONE<T>::
Normal(const TV& X) const 
{
    HELPER<T,TV> h(*this,X);
    return h.N();
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> CONE<T>::
Normal(const TV& X,const int aggregate) const 
{
    return Normal(X);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,3> CONE<T>::
Hessian(const TV& X) const
{
    HELPER<T,TV> h(*this,X);
    return h.H();
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> CONE<T>::
Principal_Curvatures(const TV& X) const
{
    HELPER<T,TV> h(*this,X);
    return Compute_Principal_Curvatures(h.N(),h.H());
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool CONE<T>::
Lazy_Inside(const TV& X) const 
{
    TV u=X-base;
    T x=u.Dot(dir);
    if(x<0) return false;
    T y=(u-dir*x).Magnitude();
    return y*height<=radius*(height-x);
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool CONE<T>::
Lazy_Outside(const TV& X) const  
{
    return !Lazy_Inside(X);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool CONE<T>::
Inside(const TV& X,const T thickness_over_two) const 
{
    return Signed_Distance(X)<-thickness_over_two;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool CONE<T>::
Outside(const TV& X,const T thickness_over_two) const  
{
    return Signed_Distance(X)>thickness_over_two;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool CONE<T>::
Boundary(const TV& X,const T thickness_over_two) const
{
    T d=Signed_Distance(X);
    return d<=thickness_over_two && d>=-thickness_over_two;
}
//#####################################################################
// Function Name
//#####################################################################
template<class T> std::string CONE<T>::
Name()
{
    return "CONE<T>";
}
//#####################################################################
template class CONE<float>;
template class CONE<double>;
}
