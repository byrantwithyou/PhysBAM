//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
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
// Function Add_Helper
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Add_Helper(const SYSTEM_MATRIX_HELPER<T>& helper)
{
    New_Block();
    data.Append_Elements(helper.data);
}
//#####################################################################
// Function Transpose_Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Add_Transpose(INTERVAL<int> range)
{
    for(int i=range.min_corner;i<range.max_corner;i++) data.Append(TRIPLE<int,int,T>(data(i).y,data(i).x,data(i).z));
}
//#####################################################################
// Function Transpose_Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Transpose(INTERVAL<int> range)
{
    for(int i=range.min_corner;i<range.max_corner;i++) exchange(data(i).x,data(i).y);
}
//#####################################################################
// Function Transpose_Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Scale(T s,INTERVAL<int> range)
{
    for(int i=range.min_corner;i<range.max_corner;i++) data(i).z*=s;
}
//#####################################################################
// Function Shift_Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Shift(int dr,int dc,INTERVAL<int> range)
{
    for(int i=range.min_corner;i<range.max_corner;i++){data(i).x+=dr;data(i).y+=dc;}
}
//#####################################################################
// Function Compact
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Compact(int rows)
{
    ARRAY<std::map<int,T> > maps(rows);
    for(int i=0;i<data.m;i++)
        maps(data(i).x)[data(i).y]+=data(i).z;

    int k=0;
    for(int i=0;i<maps.m;i++)
        for(auto it:maps(i))
            data(k++)=TRIPLE<int,int,T>(i,it.first,it.second);
    data.Resize(k);
    compacted=true;
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Set_Matrix(int m,int n,SPARSE_MATRIX_FLAT_MXN<T>& M,bool sorted) const
{
    M.Reset(n);
    M.m=m;
    M.A.Resize(data.m);
    M.offsets.Resize(m+1,init_all,0);
    for(auto i:data) M.offsets(i.x+1)++;
    int start=0;
    for(int i=0;i<M.offsets.m;i++){
        int k=M.offsets(i);
        M.offsets(i)=start;
        start+=k;}
    for(auto i:data)
        M.A(M.offsets(i.x+1)++)={i.y,i.z};
    if(!sorted) M.Sort_Entries();
}
//#####################################################################
// Function Add
//#####################################################################
template<class T> void SYSTEM_MATRIX_HELPER<T>::
Add_Matrix(const SPARSE_MATRIX_FLAT_MXN<T>& M,bool trans,int dr,int dc)
{
    New_Block();
    M.For_Each([&](int i,int j,T a){data.Append({i,j,a});});
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
namespace PhysBAM{
template struct SYSTEM_MATRIX_HELPER<float>;
template struct SYSTEM_MATRIX_HELPER<double>;
}
