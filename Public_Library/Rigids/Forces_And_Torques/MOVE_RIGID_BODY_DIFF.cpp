//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Math_Tools/Robust_Functions.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/TWIST.h>
#include <Rigids/Forces_And_Torques/MOVE_RIGID_BODY_DIFF.h>
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
    twist=dt_twist;
    frame.t=frame0.t+dt_twist.linear;
    frame.r=ROTATION<TV>::From_Rotation_Vector(dt_twist.angular)*frame0.r;
    Move_Rigid_Body_Diff_Helper(tensor,dt_twist.angular);
    tensor=Contract<1,0>(tensor,frame0.r.Rotation_Matrix());
    R=frame.r.Rotation_Matrix();
}
//#####################################################################
// Function Frame_Times
//#####################################################################
template<class TV> TV MOVE_RIGID_BODY_DIFF<TV>::
Frame_Times(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m>& dZdL,
    MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const
{
    dZdv=R;
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
    dZdv=R.Transposed();
    dZdL=-R.Transposed();
    dZdA=Contract<0>(tensor,v-frame.t);
    return frame.Inverse_Times(v);
}
//#####################################################################
// Function Rotate
//#####################################################################
template<class TV> TV MOVE_RIGID_BODY_DIFF<TV>::
Rotate(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const
{
    dZdv=R;
    dZdA=Contract<1>(tensor,v);
    return R*v;
}
//#####################################################################
// Function Inverse_Rotate
//#####################################################################
template<class TV> TV MOVE_RIGID_BODY_DIFF<TV>::
Inverse_Rotate(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const
{
    dZdv=R.Transposed();
    dZdA=Contract<0>(tensor,v);
    return R.Transpose_Times(v);
}
//#####################################################################
// Function Test_Move_Rigid_Body_Diff
//#####################################################################
template<class TV> void
Test_Move_Rigid_Body_Diff()
{
    typedef typename TV::SCALAR T;
    T eps=1e-6;
    
    RANDOM_NUMBERS<T> rand;

    FRAME<TV> frame0;
    rand.Fill_Uniform(frame0,-1,1);
    MATRIX<T,TV::m> inertia_tmp;
    rand.Fill_Uniform(inertia_tmp,-1,1);
    SYMMETRIC_MATRIX<T,TV::m> inertia0=inertia_tmp.Outer_Product_Matrix();

    TWIST<TV> w[2],dw;
    rand.Fill_Uniform(w[0],-1,1);
    rand.Fill_Uniform(dw,-eps,eps);
    w[1]=w[0]+dw;

    TV vec[2],dvec;
    rand.Fill_Uniform(vec[0],-1,1);
    rand.Fill_Uniform(dvec,-eps,eps);
    vec[1]=vec[0]+dvec;
    T dt;
    rand.Fill_Uniform(dt,.1,1);

    TV Y[2];
    MATRIX<T,TV::m> dYdV[2];
    MATRIX<T,TV::m> dYdL[2];
    MATRIX<T,TV::m,TV::SPIN::m> dYdA[2];
    TV Z[2];
    MATRIX<T,TV::m> dZdV[2];
    MATRIX<T,TV::m> dZdL[2];
    MATRIX<T,TV::m,TV::SPIN::m> dZdA[2];
    TV A[2];
    MATRIX<T,TV::m> dAdV[2];
    MATRIX<T,TV::m,TV::SPIN::m> dAdA[2];
    TV B[2];
    MATRIX<T,TV::m> dBdV[2];
    MATRIX<T,TV::m,TV::SPIN::m> dBdA[2];
    MOVE_RIGID_BODY_DIFF<TV> mr[2];
    for(int s=0;s<2;s++)
    {
        mr[s].Compute(frame0,dt*w[s]);
        Y[s]=mr[s].Frame_Times(vec[s],dYdV[s],dYdL[s],dYdA[s]);
        Z[s]=mr[s].Frame_Inverse_Times(vec[s],dZdV[s],dZdL[s],dZdA[s]);
        A[s]=mr[s].Rotate(vec[s],dAdV[s],dAdA[s]);
        B[s]=mr[s].Inverse_Rotate(vec[s],dBdV[s],dBdA[s]);
    }
    TV dY0=Y[1]-Y[0];
    dY0/=eps;
    TV dY1;
    dY1+=dYdV[0]*dvec+dYdL[0]*dw.linear*dt+dYdA[0]*dw.angular*dt;
    dY1+=dYdV[1]*dvec+dYdL[1]*dw.linear*dt+dYdA[1]*dw.angular*dt;
    dY1/=2*eps;
    LOG::printf("%P %P %P\n",dY0.Magnitude(),dY1.Magnitude(),(dY0-dY1).Magnitude());

    TV dZ0=Z[1]-Z[0];
    dZ0/=eps;
    TV dZ1;
    dZ1+=dZdV[0]*dvec+dZdL[0]*dw.linear*dt+dZdA[0]*dw.angular*dt;
    dZ1+=dZdV[1]*dvec+dZdL[1]*dw.linear*dt+dZdA[1]*dw.angular*dt;
    dZ1/=2*eps;
    LOG::printf("%P %P %P\n",dZ0.Magnitude(),dZ1.Magnitude(),(dZ0-dZ1).Magnitude());

    TV dA0=A[1]-A[0];
    dA0/=eps;
    TV dA1;
    dA1+=dAdV[0]*dvec+dAdA[0]*dw.angular*dt;
    dA1+=dAdV[1]*dvec+dAdA[1]*dw.angular*dt;
    dA1/=2*eps;
    LOG::printf("%P %P %P\n",dA0.Magnitude(),dA1.Magnitude(),(dA0-dA1).Magnitude());

    TV dB0=B[1]-B[0];
    dB0/=eps;
    TV dB1;
    dB1+=dBdV[0]*dvec+dBdA[0]*dw.angular*dt;
    dB1+=dBdV[1]*dvec+dBdA[1]*dw.angular*dt;
    dB1/=2*eps;
    LOG::printf("%P %P %P\n",dB0.Magnitude(),dB1.Magnitude(),(dB0-dB1).Magnitude());
}
template class MOVE_RIGID_BODY_DIFF<VECTOR<float,2> >;
template class MOVE_RIGID_BODY_DIFF<VECTOR<float,3> >;
template class MOVE_RIGID_BODY_DIFF<VECTOR<double,2> >;
template class MOVE_RIGID_BODY_DIFF<VECTOR<double,3> >;
template void Test_Move_Rigid_Body_Diff<VECTOR<float,2> >();
template void Test_Move_Rigid_Body_Diff<VECTOR<float,3> >();
template void Test_Move_Rigid_Body_Diff<VECTOR<double,2> >();
template void Test_Move_Rigid_Body_Diff<VECTOR<double,3> >();
}
