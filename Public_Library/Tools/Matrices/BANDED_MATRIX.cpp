//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
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
    for(int i=0;i<u.m-1;i++)
        Givens_Transpose_Times(Givens(A(i),A(i+1)),u(i),u(i+1));
    diagonal_column--;
}
//#####################################################################
// Function Backsolve
//#####################################################################
template<class T,int w> template<class T2> void BANDED_MATRIX<T,w>::
Backsolve(ARRAY<T2>& u) const
{
    assert(diagonal_column==0);
    VECTOR<T2,w> r;
    for(int i=u.m-1;i>=0;i--){
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
    // Check the bottom-left corner.
    for(int r=A.m+diagonal_column-w+1;r<A.m;r++){
        for(int i=A.m+diagonal_column-r;i<w;i++){
            if(abs(A(r)(i))>1e-14){
                ROW temp;
                for(int j=i;j<w;j++){temp(j-i)=A(r)(j);A(r)(j)=0;}
                for(int c=-diagonal_column+r+i-A.m,R=c+diagonal_column;R<r;c++,R++){
                    // leftmost entry of temp corresponds to slot A[r][c].
                    T multiplier=-temp(0)/A(R)(0); // THIS CAN BE DIVISION BY 0.
                    temp+=multiplier*A(R);
                    for(int k=0;k<w-1;k++) temp(k)=temp(k+1);
                    temp(w-1)=0;
                    u(r)+=multiplier*u(R);}
                for(int k=0;k<w;k++) A(r)(k)+=temp(k);}}}

    // Now check the top-right corner.
    int DC=w-diagonal_column-1;
    for(int r=A.m+DC-w+1;r<A.m;r++){
        for(int i=A.m+DC-r;i<w;i++){
            if(abs(A(A.m-1-r)(w-1-i))>1e-14){
                ROW temp;
                for(int j=i;j<w;j++){temp(j-i)=A(A.m-1-r)(w-1-j);A(A.m-1-r)(w-1-j)=0;}
                for(int c=-DC+r+i-A.m,R=c+DC;R<r;c++,R++){
                    T multiplier=-temp(0)/A(A.m-1-R)(w-1);
                    for(int k=0;k<w-1;k++) temp(k)=temp(k+1)+multiplier*A(A.m-1-R)(w-2-k);
                    u(A.m-1-r)+=multiplier*u(A.m-1-R);}
                for(int k=0;k<w;k++) A(A.m-1-r)(w-1-k)+=temp(k);}}}

    while(diagonal_column>0)
        Givens_Shift_Diagonal(u);
    Backsolve(u);
}
template<class T,int w> template<class T2> void BANDED_MATRIX<T,w>::
QR_Solve(ARRAY<ARRAY<T2>>& u)
{
    // Check the bottom-left corner.
    for(int r=A.m+diagonal_column-w+1;r<A.m;r++){
        for(int i=A.m+diagonal_column-r;i<w;i++){
            if(abs(A(r)(i))>1e-14){
                ROW temp;
                for(int j=i;j<w;j++){temp(j-i)=A(r)(j);A(r)(j)=0;}
                for(int c=-diagonal_column+r+i-A.m,R=c+diagonal_column;R<r;c++,R++){
                    // leftmost entry of temp corresponds to slot A[r][c].
                    T multiplier=-temp(0)/A(R)(0); // THIS CAN BE DIVISION BY 0.
                    temp+=multiplier*A(R);
                    for(int k=0;k<w-1;k++) temp(k)=temp(k+1);
                    temp(w-1)=0;
                    u(r).Copy(multiplier,u(R),u(r));}
                for(int k=0;k<w;k++) A(r)(k)+=temp(k);}}}

    // Now check the top-right corner.
    int DC=w-diagonal_column-1;
    for(int r=A.m+DC-w+1;r<A.m;r++){
        for(int i=A.m+DC-r;i<w;i++){
            if(abs(A(A.m-1-r)(w-1-i))>1e-14){
                ROW temp;
                for(int j=i;j<w;j++){temp(j-i)=A(A.m-1-r)(w-1-j);A(A.m-1-r)(w-1-j)=0;}
                for(int c=-DC+r+i-A.m,R=c+DC;R<r;c++,R++){
                    T multiplier=-temp(0)/A(A.m-1-R)(w-1);
                    for(int k=0;k<w-1;k++) temp(k)=temp(k+1)+multiplier*A(A.m-1-R)(w-2-k); 
                    u(A.m-1-r).Copy(multiplier,u(A.m-1-R),u(A.m-1-r));}
                for(int k=0;k<w;k++) A(A.m-1-r)(w-1-k)+=temp(k);}}}

    while(diagonal_column>0)
        Givens_Shift_Diagonal(u);
    Backsolve(u);
}
template class BANDED_MATRIX<float,4>;
template class BANDED_MATRIX<double,4>;
template class BANDED_MATRIX<float,3>;
template class BANDED_MATRIX<double,3>;
template void BANDED_MATRIX<float,4>::QR_Solve<VECTOR<float,2> >(ARRAY<VECTOR<float,2>,int>&);
template void BANDED_MATRIX<double,4>::QR_Solve<VECTOR<double,2> >(ARRAY<VECTOR<double,2>,int>&);
template void BANDED_MATRIX<float,3>::QR_Solve<VECTOR<float,2> >(ARRAY<VECTOR<float,2>,int>&);
template void BANDED_MATRIX<double,3>::QR_Solve<VECTOR<double,2> >(ARRAY<VECTOR<double,2>,int>&);
template void BANDED_MATRIX<double,3>::QR_Solve<double>(ARRAY<double>&);
template void BANDED_MATRIX<double,4>::QR_Solve<double>(ARRAY<double>&);
template void BANDED_MATRIX<float,3>::QR_Solve<VECTOR<float,2> >(ARRAY<ARRAY<VECTOR<float,2>,int>,int>&);
template void BANDED_MATRIX<double,3>::QR_Solve<VECTOR<double,2> >(ARRAY<ARRAY<VECTOR<double,2>,int>,int>&);
template void BANDED_MATRIX<float,3>::QR_Solve<VECTOR<float,3> >(ARRAY<ARRAY<VECTOR<float,3>,int>,int>&);
template void BANDED_MATRIX<double,3>::QR_Solve<VECTOR<double,3> >(ARRAY<ARRAY<VECTOR<double,3>,int>,int>&);
}
