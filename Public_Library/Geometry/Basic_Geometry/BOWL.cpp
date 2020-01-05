//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOWL
//##################################################################### 
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/AUTO_DIFF.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/BOWL.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
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
// Function Signed_Distance
//#####################################################################
template<class T> T BOWL<T>::
Signed_Distance(const TV& Z) const
{
    TV X=frame*Z;
    T radius=::std::hypot(X.x,X.z);
    T dX_x=radius-hole_radius,dX_y=X.y,dX_mag=::std::hypot(dX_x,dX_y);

    if((dX_x>=0 && dX_y>=0) || (dX_x<=0 && dX_y<=0)){
        T dr=dX_mag*((dX_x>=0 && dX_y>=0)?(T)1:(T)(-1))-(depth+height)*0.5;
        T dshell=abs(dr)-thickness*(T)0.5;
        if(dshell<=0 && abs(dX_y)<abs(dshell)) return -abs(dX_y);
        else return dshell;
    }
    T max_value=max(dX_x,dX_y);
    T min_value=min(dX_x,dX_y);
    if(max_value<depth) return ::std::hypot(min_value,max_value-depth);
    if(max_value>height) return ::std::hypot(min_value,max_value-height);
    return -min_value;
}
//#####################################################################
// Function Compute_Helper
//#####################################################################
template<class T,class TV> static AUTO_DIFF<T,TV>
Compute_Diff_Helper(const BOWL<T>& b,const TV& X) PHYSBAM_FLATTEN;
template<class T,class TV> static AUTO_DIFF<T,TV>
Compute_Diff_Helper(const BOWL<T>& b,const TV& X)
{
    AUTO_DIFF<T,TV> radius=hypot(AUTO_DIFF<T,TV>::From_Var(X,0),AUTO_DIFF<T,TV>::From_Var(X,2));
    AUTO_DIFF<T,TV> dX_x=radius-b.hole_radius,dX_y=AUTO_DIFF<T,TV>::From_Var(X,1),dX_mag=hypot(dX_x,dX_y);

    if((dX_x.x>=0 && dX_y.x>=0) || (dX_x.x<=0 && dX_y.x<=0)){
        AUTO_DIFF<T,TV> dr=dX_mag*((dX_x.x>=0 && dX_y.x>=0)?(T)1:(T)(-1))-(b.depth+b.height)*0.5;
        AUTO_DIFF<T,TV> dshell=abs(dr)-b.thickness*(T)0.5;
        if(dshell.x<=0 && abs(dX_y.x)<abs(dshell.x)) return -abs(dX_y);
        else return dshell;}
    AUTO_DIFF<T,TV> max_value=max(dX_x,dX_y);
    AUTO_DIFF<T,TV> min_value=min(dX_x,dX_y);
    if(max_value.x<b.depth) return hypot(min_value,max_value-b.depth);
    if(max_value.x>b.height) return hypot(min_value,max_value-b.height);
    return -min_value;
}
//#####################################################################
// Function Compute_Helper
//#####################################################################
template<class T,class TV> static AUTO_HESS<T,TV>
Compute_Hess_Helper(const BOWL<T>& b,const TV& X) PHYSBAM_FLATTEN;
template<class T,class TV> static AUTO_HESS<T,TV>
Compute_Hess_Helper(const BOWL<T>& b,const TV& X)
{
    AUTO_HESS<T,TV> radius=hypot(AUTO_HESS<T,TV>::From_Var(X,0),AUTO_HESS<T,TV>::From_Var(X,2));
    AUTO_HESS<T,TV> dX_x=radius-b.hole_radius,dX_y=AUTO_HESS<T,TV>::From_Var(X,1),dX_mag=hypot(dX_x,dX_y);

    if((dX_x.x>=0 && dX_y.x>=0) || (dX_x.x<=0 && dX_y.x<=0)){
        AUTO_HESS<T,TV> dr=dX_mag*((dX_x.x>=0 && dX_y.x>=0)?(T)1:(T)(-1))-(b.depth+b.height)*0.5;
        AUTO_HESS<T,TV> dshell=abs(dr)-b.thickness*(T)0.5;
        if(dshell.x<=0 && abs(dX_y.x)<abs(dshell.x)) return -abs(dX_y);
        else return dshell;}
    AUTO_HESS<T,TV> max_value=max(dX_x,dX_y);
    AUTO_HESS<T,TV> min_value=min(dX_x,dX_y);
    if(max_value.x<b.depth) return hypot(min_value,max_value-b.depth);
    if(max_value.x>b.height) return hypot(min_value,max_value-b.height);
    return -min_value;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> BOWL<T>::
Normal(const TV& Z) const
{
    TV X=frame*Z;
    return frame.r.Rotate(Compute_Diff_Helper(*this,X).dx);
}
//#####################################################################
// Function Suface
//#####################################################################
template<class T> VECTOR<T,3> BOWL<T>::
Surface(const TV& Z) const 
{
    TV X=frame*Z;
    AUTO_DIFF<T,TV> diff=Compute_Diff_Helper(*this,X);
    return Z-frame.r.Rotate(diff.dx)*diff.x;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> BOWL<T>::
Normal(const TV& Z,const int aggregate) const 
{
    return Normal(Z);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,3> BOWL<T>::
Hessian(const TV& Z) const
{
    TV X=frame*Z;
    auto H=Compute_Hess_Helper(*this,X).ddx;
    auto R=frame.r.Rotation_Matrix();
    return SYMMETRIC_MATRIX<T,TV::m>::Conjugate(R,H);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> BOWL<T>::
Principal_Curvatures(const TV& Z) const
{
    TV X=frame*Z;
    return ::PhysBAM::Principal_Curvatures(Compute_Hess_Helper(*this,X));
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool BOWL<T>::
Lazy_Inside(const TV& Z) const 
{
    return Signed_Distance(Z)<0;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool BOWL<T>::
Lazy_Outside(const TV& Z) const  
{
    return !Lazy_Inside(Z);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool BOWL<T>::
Inside(const TV& Z,const T thickness_over_two) const 
{
    return Signed_Distance(Z)<=-thickness_over_two;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool BOWL<T>::
Outside(const TV& Z,const T thickness_over_two) const  
{
    return Signed_Distance(Z)>=thickness_over_two;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool BOWL<T>::
Boundary(const TV& Z,const T thickness_over_two) const
{
    T sd=Signed_Distance(Z);
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
template class BOWL<double>;
}
