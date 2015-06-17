//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/Is_NaN.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <cassert>
namespace PhysBAM{
namespace{
//#####################################################################
// Function Givens_Transpose_Times
//#####################################################################
template<class T,class T2> void
Givens_Transpose_Times(const VECTOR<T,2>& g,T2& x,T2& y)
{
    T2 xx=g.x*x+g.y*y;
    y=-g.y*x+g.x*y;
    x=xx;
}
template<class T,class T2> void
Givens_Transpose_Times(const VECTOR<T,2>& g,ARRAY<T2>& x,ARRAY<T2>& y)
{
    ARRAY<T2> xx(x.m);
    xx.Copy(g.x,x,g.y,y);
    y.Copy(g.x,y,-g.y,x);
    x=xx;
}
}
//#####################################################################
// Function Givens
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
    for(int i=0;i<u.m-1;i++){
        VECTOR<T,2> g=Givens(A(i),A(i+1));
        Givens_Transpose_Times(g,u(i),u(i+1));
        if(wrap) Givens_Transpose_Times(g,extra(i),extra(i+1));}
    diagonal_column--;
}
template<class T2,int d> void 
Init_Helper(VECTOR<T2,d>& a,const ARRAY<T2>& u)
{}
template<class T2,int d> void 
Init_Helper(VECTOR<ARRAY<T2>,d>& a,const ARRAY<ARRAY<T2>>& u)
{
    for(int i=0;i<d;i++)
        a(i).Resize(u(0).m);
}
//#####################################################################
// Function Backsolve
//#####################################################################
template<class T,int w> template<class T2> void BANDED_MATRIX<T,w>::
Backsolve(ARRAY<T2>& u) const
{
    assert(diagonal_column==0);
    VECTOR<T2,w> r;
    Init_Helper(r,u);
    for(int i=u.m-1;i>=0;i--){
        VECTOR<T2,w-1> q=r.Remove_Index(w-1);
        const ROW& ar=A(i);
        for(int j=0;j<w-1;j++) u(i)-=ar(j+1)*q(j);
        u(i)/=A(i)(0);
        r=q.Prepend(u(i));}
}
//#####################################################################
// Function QR_Solve
//#####################################################################
template<class T,int w> template<class T2> void BANDED_MATRIX<T,w>::
QR_Solve(ARRAY<T2>& u)
{
    int odc=diagonal_column;
    if(wrap){
        // Move the top-right corner to the rhs
        extra.Resize(u.m);
        for(int i=0;i<extra.m;i++)
            extra(i).Resize(w-1);
        for(int i=0;i<odc;i++){
            for(int j=0;j<odc-i;j++){
                extra(i)(j+i+w-1-odc)=A(i)(j);
                A(i)(j)=0;}}
        // Move the bottom-left corner to the rhs
        for(int i=A.m-w+1+odc;i<A.m;i++){
            for(int j=odc+A.m-i;j<w;j++){
                extra(i)(j+i-A.m-odc)=A(i)(j);
                A(i)(j)=0;}}}

    while(diagonal_column>0)
        Givens_Shift_Diagonal(u);
    Backsolve(u);

    if(wrap){
        Backsolve(extra);
        MATRIX<T,w-1> ICE;
        for(int i=0;i<w-1;i++){
            for(int j=0;j<odc;j++)
                ICE(j,i)=extra(j)(i)+((i==j)?1:0);
            for(int j=odc;j<w-1;j++)
                ICE(j,i)=extra(j+u.m-w+1)(i)+((i==j)?1:0);}
        ICE.Invert();
        VECTOR<T2,w-1> a;
        Init_Helper(a,u);
        for(int i=0;i<w-1;i++){
            for(int j=0;j<odc;j++)
                a(i)+=ICE(i,j)*u(j);
            for(int j=odc;j<w-1;j++)
                a(i)+=ICE(i,j)*u(j+u.m-w+1);}
        for(int i=0;i<u.m;i++)
            for(int j=0;j<w-1;j++)
                u(i)-=extra(i)(j)*a(j);}
}

template class BANDED_MATRIX<float,4>;
template class BANDED_MATRIX<double,4>;
template class BANDED_MATRIX<float,3>;
template class BANDED_MATRIX<double,3>;
template class BANDED_MATRIX<float,2>;
template class BANDED_MATRIX<double,2>;
template void BANDED_MATRIX<float,4>::QR_Solve<VECTOR<float,2> >(ARRAY<VECTOR<float,2>,int>&);
template void BANDED_MATRIX<double,4>::QR_Solve<VECTOR<double,2> >(ARRAY<VECTOR<double,2>,int>&);
template void BANDED_MATRIX<float,3>::QR_Solve<VECTOR<float,2> >(ARRAY<VECTOR<float,2>,int>&);
template void BANDED_MATRIX<double,3>::QR_Solve<VECTOR<double,2> >(ARRAY<VECTOR<double,2>,int>&);
template void BANDED_MATRIX<double,2>::QR_Solve<double>(ARRAY<double>&);
template void BANDED_MATRIX<double,3>::QR_Solve<double>(ARRAY<double>&);
template void BANDED_MATRIX<double,4>::QR_Solve<double>(ARRAY<double>&);
template void BANDED_MATRIX<float,3>::QR_Solve<ARRAY<VECTOR<float,2>> >(ARRAY<ARRAY<VECTOR<float,2>,int>,int>&);
template void BANDED_MATRIX<double,3>::QR_Solve<ARRAY<VECTOR<double,2>> >(ARRAY<ARRAY<VECTOR<double,2>,int>,int>&);
template void BANDED_MATRIX<float,3>::QR_Solve<ARRAY<VECTOR<float,3>> >(ARRAY<ARRAY<VECTOR<float,3>,int>,int>&);
template void BANDED_MATRIX<double,3>::QR_Solve<ARRAY<VECTOR<double,3>> >(ARRAY<ARRAY<VECTOR<double,3>,int>,int>&);
}
