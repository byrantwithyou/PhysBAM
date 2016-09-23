//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Geoffrey Irving, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONE
//##################################################################### 
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/AUTO_DIFF.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Auto_Diff/AUTO_NO_DIFF.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/CONE.h>
namespace PhysBAM{
template<class DT,class DTV,class T,class TV>
inline DT Signed_Distance_Helper(const CONE<T>& cone,const TV& X)
{
    DTV u=DTV::From_Var(X)-cone.base;
    DT x=u.Dot(cone.dir);
    DT y=(u-cone.dir*x).Magnitude();

    if(x.x<=0 && y.x<=cone.radius) return -x;
    if((y.x-cone.radius)*cone.radius>cone.height*x.x) return sqrt(sqr(x)+sqr(y-cone.radius));
    if((x.x-cone.height)*cone.height>cone.radius*y.x) return sqrt(sqr(x-cone.height)+sqr(y));
    return max(-x,((x-cone.height)*cone.radius+y*cone.height)/sqrt(sqr(cone.radius)+sqr(cone.height)));
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
    return Signed_Distance_Helper<AUTO_NO_DIFF<T,TV>,AUTO_NO_DIFF<TV,TV> >(*this,X).x;
}
//#####################################################################
// Function Suface
//#####################################################################
template<class T> VECTOR<T,3> CONE<T>::
Surface(const TV& X) const 
{
    AUTO_DIFF<T,TV> phi=Signed_Distance_Helper<AUTO_DIFF<T,TV>,AUTO_DIFF<TV,TV> >(*this,X);
    return X-phi.x*phi.dx;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> CONE<T>::
Normal(const TV& X) const 
{
    return Signed_Distance_Helper<AUTO_DIFF<T,TV>,AUTO_DIFF<TV,TV> >(*this,X).dx;
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
    return Signed_Distance_Helper<AUTO_HESS<T,TV>,AUTO_HESS<TV,TV> >(*this,X).ddx;
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> CONE<T>::
Principal_Curvatures(const TV& X) const
{
    return ::PhysBAM::Principal_Curvatures(Signed_Distance_Helper<AUTO_HESS<T,TV>,AUTO_HESS<TV,TV> >(*this,X));
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
