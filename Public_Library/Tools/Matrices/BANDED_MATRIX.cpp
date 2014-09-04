//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
#include <cassert>
namespace PhysBAM{
namespace{
template<class T,class T2> void
Givens_Transpose_Times(const VECTOR<T,2>& g,T2& x,T2& y)
{
    T2 xx=g.x*x+g.y*y;
    y=-g.y*x+g.x*y;
    x=xx;
}
}
//#####################################################################
// Function Givens_Shift_Diagonal
//#####################################################################
template<class T,int w> VECTOR<T,2> BANDED_MATRIX<T,w>::
Givens(ROW& x,ROW& y)
{
    TV gt(x(0),y(0));
    T z=gt.Normalize();
    x(0)=z;
    for(int i=1;i<w;i++){
        Givens_Transpose_Times(gt,x(i),y(i));
        y(i-1)=y(i);}
    y(w-1)=0;
    return gt;
}
//#####################################################################
// Function Givens_Shift_Diagonal
//#####################################################################
template<class T,int w> template<class T2> void BANDED_MATRIX<T,w>::
Givens_Shift_Diagonal(ARRAY<T2>& u)
{
    assert(diagonal_column>0);
    A(0)=A(0).Remove_Index(0).Append(0);
    for(int i=0;i<u.m-1;i++)
    {
        LOG::printf("%P %P  ->  ",A(i),A(i+1));
        Givens_Transpose_Times(Givens(A(i),A(i+1)),u(i),u(i+1));
        LOG::printf("%P %P\n",A(i),A(i+1));
    }
    diagonal_column--;
}
//#####################################################################
// Function Givens_Shift_Diagonal
//#####################################################################
template<class T,int w> template<class T2> void BANDED_MATRIX<T,w>::
Backsolve(ARRAY<T2>& u) const
{
    assert(diagonal_column==0);
    VECTOR<T2,w> r;
    for(int i=u.m-1;i>=0;i--){
        LOG::printf("%P %P\n",A(i),r);
        VECTOR<T2,w-1> q=r.Remove_Index(w-1);
        T2 rhs=u(i);
        const ROW& ar=A(i);
        for(int j=0;j<w-1;j++)
            rhs-=ar(j+1)*q(j);
        rhs/=A(i)(0);
        u(i)=rhs;
        r=q.Prepend(rhs);}
}
//#####################################################################
// Function QR_Solve
//#####################################################################
template<class T,int w> template<class T2> void BANDED_MATRIX<T,w>::
QR_Solve(ARRAY<T2>& u)
{
    while(diagonal_column>0)
        Givens_Shift_Diagonal(u);
    Backsolve(u);
}
template class BANDED_MATRIX<float,4>;
template class BANDED_MATRIX<double,4>;
template void BANDED_MATRIX<float,4>::QR_Solve<VECTOR<float,2> >(ARRAY<VECTOR<float,2>,int>&);
template void BANDED_MATRIX<double,4>::QR_Solve<VECTOR<double,2> >(ARRAY<VECTOR<double,2>,int>&);
}
