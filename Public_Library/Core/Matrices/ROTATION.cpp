//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROTATION
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/Robust_Functions.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/TWIST.h>
namespace PhysBAM{
using ::std::sin;
using ::std::cos;
//#####################################################################
// Function From_Rotated_Vector
//#####################################################################
template<class T>  ROTATION<VECTOR<T,2> > ROTATION<VECTOR<T,2> >::
From_Rotated_Vector(const TV& v1,const TV& v2)
{
    return ROTATION<TV>(std::complex<T>(v1.x,-v1.y)*std::complex<T>(v2.x,v2.y)).Normalized();
}
//#####################################################################
// Function Scale_Angle
//#####################################################################
template<class T>  ROTATION<VECTOR<T,2> > ROTATION<VECTOR<T,2> >::
Scale_Angle(const T a) const
{
    return ROTATION<TV>::From_Rotation_Vector(a*Rotation_Vector());
}
//#####################################################################
// Function Spherical_Linear_Interpolation
//#####################################################################
template<class T>  ROTATION<VECTOR<T,2> > ROTATION<VECTOR<T,2> >::
Spherical_Linear_Interpolation(const ROTATION<TV>& r1,const ROTATION<TV>& r2,const T t)
{
    return r1*(r1.Inverse()*r2).Scale_Angle(t);
}
//#####################################################################
// Function Average_Rotation
//#####################################################################
template<class T>  ROTATION<VECTOR<T,2> > ROTATION<VECTOR<T,2> >::
Average_Rotation(const ARRAY<ROTATION<TV> >& rotations)
{
    std::complex<T> sum;
    for(int i=0;i<rotations.m;i++) sum+=rotations(i).c;
    return ROTATION<TV>(sum/abs(sum));
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> ROTATION<VECTOR<T,3> >::
ROTATION(const T angle,const TV& direction)
    :q(cos((T).5*angle),direction)
{
    q.v.Normalize();
    q.v*=sin((T).5*angle);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> ROTATION<VECTOR<T,3> >::
ROTATION(const MATRIX<T,3>& A) // matches A with a quaternion
{
    T trace=1+A(0,0)+A(1,1)+A(2,2);// trace=4*cos^2(theta/1)
    if(trace>1){q.s=(T).5*sqrt(trace);q.v.x=A(2,1)-A(1,2);q.v.y=A(0,2)-A(2,0);q.v.z=A(1,0)-A(0,1);q.v*=(T).25/q.s;}
    else{int i=(A(0,0)>A(1,1)) ? 0:1;i=(A(i,i)>A(2,2)) ? i:2; // set i to be the index of the dominating diagonal term
        switch(i){
            case 0:q.v.x=T(.5)*sqrt(1+A(0,0)-A(1,1)-A(2,2));q.v.y=T(.25)*(A(1,0)+A(0,1))/q.v.x;q.v.z=T(.25)*(A(0,2)+A(2,0))/q.v.x;q.s=T(.25)*(A(2,1)-A(1,2))/q.v.x;break;
            case 1:q.v.y=T(.5)*sqrt(1-A(0,0)+A(1,1)-A(2,2));q.v.x=T(.25)*(A(1,0)+A(0,1))/q.v.y;q.v.z=T(.25)*(A(2,1)+A(1,2))/q.v.y;q.s=T(.25)*(A(0,2)-A(2,0))/q.v.y;break;
            case 2:default:q.v.z=T(.5)*sqrt(1-A(0,0)-A(1,1)+A(2,2));q.v.x=T(.25)*(A(0,2)+A(2,0))/q.v.z;q.v.y=T(.25)*(A(2,1)+A(1,2))/q.v.z;q.s=T(.25)*(A(1,0)-A(0,1))/q.v.z;break;}}
    Normalize();
}
//#####################################################################
// Function Rotation_Vector
//#####################################################################
template<class T> VECTOR<T,3> ROTATION<VECTOR<T,3> >::
Rotation_Vector() const
{
    return 2*atan2_y_x_over_y(q.v.Magnitude(),abs(q.s))*(q.s<0?-q.v:q.v);
}
//#####################################################################
// Function From_Rotation_Vector
//#####################################################################
template<class T>  ROTATION<VECTOR<T,3> > ROTATION<VECTOR<T,3> >::
From_Rotation_Vector(const TV& v)
{
    T magnitude=v.Magnitude();
    ROTATION<TV> r;
    r.q.s=cos((T).5*magnitude);r.q.v=(T).5*sinc((T).5*magnitude)*v;
    return r;
}
//#####################################################################
// Function From_Euler_Angles
//#####################################################################
template<class T>  ROTATION<VECTOR<T,3> > ROTATION<VECTOR<T,3> >::
From_Euler_Angles(const T euler_angle_x,const T euler_angle_y,const T euler_angle_z) // rotation about fixed axes x, then y, then z
{
    T half_x_angle=(T).5*euler_angle_x,half_y_angle=(T).5*euler_angle_y,half_z_angle=(T).5*euler_angle_z;
    T cx=cos(half_x_angle),sx=sin(half_x_angle),cy=cos(half_y_angle),sy=sin(half_y_angle),cz=cos(half_z_angle),sz=sin(half_z_angle);
    return ROTATION<TV>(cx*cy*cz+sx*sy*sz,sx*cy*cz-cx*sy*sz,cx*sy*cz+sx*cy*sz,cx*cy*sz-sx*sy*cz);
}
//#####################################################################
// Function Euler_Angles
//#####################################################################
template<class T>  void ROTATION<VECTOR<T,3> >::
Euler_Angles(T& euler_angle_x,T& euler_angle_y,T& euler_angle_z) const
{
    MATRIX<T,3> R(Rotation_Matrix());
    T cos_beta_squared=sqr(R.x[0])+sqr(R.x[1]);
    if(cos_beta_squared<1e-30){
        euler_angle_x=0;
        euler_angle_z=-atan2(R.x[3],R.x[4]);
        euler_angle_y=R.x[2]<0?(T).5*(T)pi:-(T).5*(T)pi;}
    else{
        euler_angle_x=atan2(R.x[5],R.x[8]); // between -pi and pi
        euler_angle_y=atan2(-R.x[2],sqrt(cos_beta_squared)); // between -pi/2 and pi/2
        euler_angle_z=atan2(R.x[1],R.x[0]);} // between -pi and pi
}
//#####################################################################
// Function Rotated_X_Axis
//#####################################################################
template<class T>  VECTOR<T,3> ROTATION<VECTOR<T,3> >::
Rotated_X_Axis() const // Q*(1,0,0)
{
    T vy2=sqr(q.v.y),vz2=sqr(q.v.z),vxvy=q.v.x*q.v.y,vxvz=q.v.x*q.v.z,svy=q.s*q.v.y,svz=q.s*q.v.z;
    return TV(1-2*(vy2+vz2),2*(vxvy+svz),2*(vxvz-svy));
}
//#####################################################################
// Function Rotated_Y_Axis
//#####################################################################
template<class T>  VECTOR<T,3> ROTATION<VECTOR<T,3> >::
Rotated_Y_Axis() const // Q*(0,1,0)
{
    T vx2=sqr(q.v.x),vz2=sqr(q.v.z),vxvy=q.v.x*q.v.y,vyvz=q.v.y*q.v.z,svx=q.s*q.v.x,svz=q.s*q.v.z;
    return TV(2*(vxvy-svz),1-2*(vx2+vz2),2*(vyvz+svx));
}
//#####################################################################
// Function Rotated_Z_Axis
//#####################################################################
template<class T>  VECTOR<T,3> ROTATION<VECTOR<T,3> >::
Rotated_Z_Axis() const // Q*(0,0,1)
{
    T vx2=sqr(q.v.x),vy2=sqr(q.v.y),vxvz=q.v.x*q.v.z,vyvz=q.v.y*q.v.z,svx=q.s*q.v.x,svy=q.s*q.v.y;
    return TV(2*(vxvz+svy),2*(vyvz-svx),1-2*(vx2+vy2));
}
//#####################################################################
// Function Get_Rotated_Frame
//#####################################################################
template<class T> void ROTATION<VECTOR<T,3> >::
Get_Rotated_Frame(TV& x_axis,TV& y_axis,TV& z_axis) const
{
    MATRIX<T,3> M=Rotation_Matrix();
    x_axis=M.Column(0);
    y_axis=M.Column(1);
    z_axis=M.Column(2);
}
//#####################################################################
// Function Get_Angle_Axis
//#####################################################################
template<class T> void ROTATION<VECTOR<T,3> >::
Get_Angle_Axis(T& angle,TV& axis) const
{
    axis=q.s<0?-q.v:q.v;
    angle=2*atan2(axis.Normalize(),abs(q.s));
}
//#####################################################################
// Function Angle
//#####################################################################
template<class T> T ROTATION<VECTOR<T,3> >::
Angle() const
{
    return 2*atan2(q.v.Magnitude(),abs(q.s));
}
//#####################################################################
// Function Get_Axis
//#####################################################################
template<class T>  VECTOR<T,3> ROTATION<VECTOR<T,3> >::
Get_Axis() const
{
    return (q.s<0?-q.v:q.v).Normalized();
}
//#####################################################################
// Function Rotation_Matrix
//#####################################################################
template<class T>  MATRIX<T,3> ROTATION<VECTOR<T,3> >::
Rotation_Matrix() const // 18 mult and 12 add/sub
{
    T vx2=sqr(q.v.x),vy2=sqr(q.v.y),vz2=sqr(q.v.z),vxvy=q.v.x*q.v.y,vxvz=q.v.x*q.v.z,vyvz=q.v.y*q.v.z,svx=q.s*q.v.x,svy=q.s*q.v.y,svz=q.s*q.v.z;
    return MATRIX<T,3>(1-2*(vy2+vz2),2*(vxvy+svz),2*(vxvz-svy),2*(vxvy-svz),1-2*(vx2+vz2),2*(vyvz+svx),2*(vxvz+svy),2*(vyvz-svx),1-2*(vx2+vy2));
}
//#####################################################################
// Function Scale_Angle
//#####################################################################
template<class T>  ROTATION<VECTOR<T,3> > ROTATION<VECTOR<T,3> >::
Scale_Angle(const T a) const
{
    T angle;TV axis;
    Get_Angle_Axis(angle,axis);
    return ROTATION<TV>(a*angle,axis);
}
//#####################################################################
// Function Spherical_Linear_Interpolation
//#####################################################################
template<class T>  ROTATION<VECTOR<T,3> > ROTATION<VECTOR<T,3> >::
Spherical_Linear_Interpolation(const ROTATION<TV>& r1,const ROTATION<TV>& r2,const T t)
{
    return r1*(r1.Inverse()*r2).Scale_Angle(t);
}
//#####################################################################
// Function Average_Rotation
//#####################################################################
template<class T>  ROTATION<VECTOR<T,3> > ROTATION<VECTOR<T,3> >::
Average_Rotation(const ARRAY<ROTATION<TV> >& rotations)
{
    if(rotations.m==0) return ROTATION<TV>();
    ARRAY<ROTATION<TV> > r(rotations);
    for(int i=1;i<r.m;i+=2) r.Append(Spherical_Linear_Interpolation(r(i),r(i+1),(T).5));
    return r.Last();
}
//#####################################################################
// Function From_Rotated_Vector
//#####################################################################
template<class T>  ROTATION<VECTOR<T,3> > ROTATION<VECTOR<T,3> >::
From_Rotated_Vector(const TV& initial_vector,const TV& final_vector)
{
    TV u=initial_vector.Normalized(),v=final_vector.Normalized(),w=u.Cross(v);
    T d=u.Dot(v);
    if(d>=-(T).5) return QUATERNION<T>(w.Prepend(1+d).Normalized());
    T m=w.Normalize();
    w.Project_Orthogonal_To_Unit_Direction(u);
    T epsilon=std::numeric_limits<T>::epsilon();
    if(w.Magnitude_Squared()<epsilon)
        w=u.Orthogonal_Vector().Projected_Orthogonal_To_Unit_Direction(v).Normalized();
    return QUATERNION<T>(w.Prepend(m/(1-d)).Normalized());
}
//#####################################################################
// Function Random_Fill_Uniform_
//#####################################################################
template<class T> void
Random_Fill_Uniform_(RANDOM_NUMBERS<T>& rand,ROTATION<VECTOR<T,1> >& r)
{
}
//#####################################################################
// Function Random_Fill_Uniform_
//#####################################################################
template<class T> void
Random_Fill_Uniform_(RANDOM_NUMBERS<T>& rand,ROTATION<VECTOR<T,2> >& r)
{
    VECTOR<T,2> v=rand.template Get_Direction<VECTOR<T,2> >();
    r.c=std::complex<T>(v.x,v.y);
}
//#####################################################################
// Function Random_Fill_Uniform_
//#####################################################################
template<class T> void
Random_Fill_Uniform_(RANDOM_NUMBERS<T>& rand,ROTATION<VECTOR<T,3> >& r)
{
    r.q=QUATERNION<T>(rand.template Get_Direction<VECTOR<T,4> >());
}
//#####################################################################
// Function Random_Fill_Uniform
//#####################################################################
template<class T,int d> void
Random_Fill_Uniform(RANDOM_NUMBERS<T>& rand,ROTATION<VECTOR<T,d> >& r)
{
    Random_Fill_Uniform_(rand,r);
}
//#####################################################################
template class ROTATION<VECTOR<float,2> >;
template class ROTATION<VECTOR<float,3> >;
template class ROTATION<VECTOR<double,2> >;
template class ROTATION<VECTOR<double,3> >;
template void Random_Fill_Uniform<double,1>(RANDOM_NUMBERS<double>&,
    ROTATION<VECTOR<double,1> >&);
template void Random_Fill_Uniform<double,2>(RANDOM_NUMBERS<double>&,
    ROTATION<VECTOR<double,2> >&);
template void Random_Fill_Uniform<double,3>(RANDOM_NUMBERS<double>&,
    ROTATION<VECTOR<double,3> >&);
template void Random_Fill_Uniform<float,1>(RANDOM_NUMBERS<float>&,
    ROTATION<VECTOR<float,1> >&);
template void Random_Fill_Uniform<float,2>(RANDOM_NUMBERS<float>&,
    ROTATION<VECTOR<float,2> >&);
template void Random_Fill_Uniform<float,3>(RANDOM_NUMBERS<float>&,
    ROTATION<VECTOR<float,3> >&);
}
