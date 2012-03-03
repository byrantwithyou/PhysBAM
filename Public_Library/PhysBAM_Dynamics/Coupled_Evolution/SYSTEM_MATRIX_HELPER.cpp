//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
#include <map>
using namespace PhysBAM;
//#####################################################################
// Function Add_Matrix
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Add_Matrix(const SYSTEM_MATRIX_BASE<T>& base,bool trans,int dr,int dc)
{
    New_Block();
    base.Add_Raw_Matrix(data);
    if(trans) Transpose();
    if(dr || dc) Shift(dr,dc);
}
//#####################################################################
// Function Transpose_Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Transpose()
{
    for(int i=start;i<data.m;i++) exchange(data(i).x,data(i).y);
}
//#####################################################################
// Function Transpose_Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Scale(T s)
{
    for(int i=start;i<data.m;i++) data(i).z*=s;
}
//#####################################################################
// Function Shift_Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Shift(int dr,int dc)
{
    for(int i=start;i<data.m;i++){data(i).x+=dr;data(i).y+=dc;}
}
//#####################################################################
// Function Compact
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Compact(int rows, T tol)
{
    ARRAY<std::map<int,T> > maps(rows);
    for(int i=0;i<data.m;i++)
        maps(data(i).x)[data(i).y]+=data(i).z;

    int k=0;
    for(int i=0;i<maps.m;i++)
        for(typename std::map<int,T>::iterator it=maps(i).begin();it!=maps(i).end();it++)
            if(abs(it->second)>=tol)
                data(k++)=TRIPLE<int,int,T>(i,it->first,it->second);
    data.Resize(k);
    compacted=true;
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Set_Matrix(int m,int n,SPARSE_MATRIX_FLAT_MXN<T>& M, T tol)
{
    if(!compacted) Compact(m, tol);
    M.Reset(n);
    M.m=m;
    M.A.Remove_All();
    M.offsets.Resize(m+1);
    for(int i=0;i<data.m;i++) M.offsets(data(i).x+1)++;
    for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);
    for(int i=0;i<data.m;i++) M.A.Append(SPARSE_MATRIX_ENTRY<T>(data(i).y,data(i).z));
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Set_Matrix(int n,SPARSE_MATRIX_FLAT_NXN<T>& M, T tol)
{
    if(!compacted) Compact(n, tol);
    M.Reset();
    M.n=n;
    M.A.Remove_All();
    M.offsets.Resize(n+1);
    for(int i=0;i<data.m;i++) M.offsets(data(i).x+1)++;
    for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);
    for(int i=0;i<data.m;i++) M.A.Append(SPARSE_MATRIX_ENTRY<T>(data(i).y,data(i).z));
}
//#####################################################################
// Function Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Add_Matrix(const SPARSE_MATRIX_FLAT_MXN<T>& M,bool trans,int dr,int dc)
{
    New_Block();
    for(int i=0;i<M.m;i++){
        int s=M.offsets(i),e=M.offsets(i+1);
        for(int j=s;j<e;j++)
            data.Append(TRIPLE<int,int,T>(i,M.A(j).j,M.A(j).a));}
    if(trans) Transpose();
    if(dr || dc) Shift(dr,dc);
}
//#####################################################################
// Function Base_To_Matrix
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Base_To_Matrix(int m,int n,const SYSTEM_MATRIX_BASE<T>& base,SPARSE_MATRIX_FLAT_MXN<T>& M,bool tranpose)
{
    SYSTEM_MATRIX_HELPER<T> helper;
    helper.Add_Matrix(base);
    if(tranpose) helper.Transpose();
    helper.Set_Matrix(m,n,M);
}
template struct SYSTEM_MATRIX_HELPER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template struct SYSTEM_MATRIX_HELPER<double>;
#endif
