//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Math_Tools/Robust_Functions.h>
#include <Core/Vectors/TWIST.h>
#include <Hybrid_Methods/Forces/MOVE_RIGID_BODY_DIFF.h>
namespace PhysBAM{
namespace{
template<class T>
void Move_Rigid_Body_Diff_Helper(TENSOR<T,3>& t,VECTOR<T,3> u)
{
    T th=u.Normalize();
    MATRIX<T,3> K=MATRIX<T,3>::Cross_Product_Matrix(u);
    SYMMETRIC_MATRIX<T,3> A=(T)1-Outer_Product(u);
    t=Tensor_Product<2>(cos(th)*K+sin(th)*K*K,u);
    t+=Contract<2,1>(PERMUTATION_TENSOR<T>(-sinc(th)),A);
    t+=Twice_Symmetric_Part<0,1>(Tensor_Product<1>(A,one_minus_cos_x_over_x(th)*u));
}
template<class T>
void Move_Rigid_Body_Diff_Helper(TENSOR<T,2,2,1>& t,VECTOR<T,1> u)
{
    T c=cos(u.x),s=sin(u.x);
    t(0,0,0)=-s;
    t(0,1,0)=-c;
    t(1,0,0)=c;
    t(1,1,0)=-s;
}
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MOVE_RIGID_BODY_DIFF<TV>::
Compute(const FRAME<TV>& frame0,const TWIST<TV>& dt_twist)
{
    frame.t=frame0.t+dt_twist.linear;
    frame.r=ROTATION<TV>::From_Rotation_Vector(dt_twist.angular)*frame0.r;
    Move_Rigid_Body_Diff_Helper(tensor,dt_twist.angular);
    tensor=Contract<1,0>(tensor,frame0.r.Rotation_Matrix());
}
//#####################################################################
// Function Frame_Times
//#####################################################################
template<class TV> TV MOVE_RIGID_BODY_DIFF<TV>::
Frame_Times(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m>& dZdL,
    MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const
{
    dZdv=frame.r.Rotation_Matrix();
    dZdL=MATRIX<T,TV::m>()+1;
    dZdA=Contract<1>(tensor,v);
    return frame*v;
}
//#####################################################################
// Function Frame_Inverse_Times
//#####################################################################
template<class TV> TV MOVE_RIGID_BODY_DIFF<TV>::
Frame_Inverse_Times(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m>& dZdL,
    MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const
{
    dZdv=frame.r.Rotation_Matrix().Transposed();
    dZdL=-frame.r.Rotation_Matrix().Transposed();
    dZdA=Contract<0>(tensor,v-frame.t);
    return frame.Inverse_Times(v);
}
template class MOVE_RIGID_BODY_DIFF<VECTOR<float,2> >;
template class MOVE_RIGID_BODY_DIFF<VECTOR<float,3> >;
template class MOVE_RIGID_BODY_DIFF<VECTOR<double,2> >;
template class MOVE_RIGID_BODY_DIFF<VECTOR<double,3> >;
}
