//#####################################################################
// Copyright 2007-2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//###################################################################
#ifndef __TENSOR_TESTS_DEFINITIONS__
#define __TENSOR_TESTS_DEFINITIONS__
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Tensors/TENSOR.h>
#include "TENSOR_TESTS.h"
using namespace PhysBAM;
//#####################################################################
// Function Assert_Zero
//#####################################################################
template<class T> template<class T_VECTOR> bool TENSOR_TESTS<T>::
Assert_Zero_Vector(const T_VECTOR& a,const char* str,T tolerance) const
{
    bool b=a.Max_Abs()<=tolerance;
    if(!b){
        status=false;
        LOG::cout<<"FAIL: "<<str<<"  "<<a.Max_Abs()<<"  "<<tolerance<<std::endl;}
    return b;
}
//#####################################################################
// Function Assert_Equal
//#####################################################################
template<class T> template<class T_VECTOR0,class T_VECTOR1> bool TENSOR_TESTS<T>::
Assert_Equal_Vector(const T_VECTOR0& a,const T_VECTOR1& b,const char* str,T tolerance) const
{
    bool k=Assert_Zero_Vector(a-b,str,tolerance);
    return k;
}
//#####################################################################
// Function Assert_Zero
//#####################################################################
template<class T> template<class T_MATRIX> bool TENSOR_TESTS<T>::
Assert_Zero_Matrix(const T_MATRIX& a,const char* str,T tolerance) const
{
    bool b=a.Max_Abs()<=tolerance;
    if(!b){
        status=false;
        LOG::cout<<"FAIL: "<<str<<"  "<<a.Max_Abs()<<"  "<<tolerance<<std::endl;}
    return b;
}
//#####################################################################
// Function Assert_Equal
//#####################################################################
template<class T> template<class T_MATRIX0,class T_MATRIX1> bool TENSOR_TESTS<T>::
Assert_Equal_Matrix(const T_MATRIX0& a,const T_MATRIX1& b,const char* str,T tolerance) const
{
    bool k=Assert_Zero_Matrix(a-b,str,tolerance);
    return k;
}
//#####################################################################
// Function Assert_Zero
//#####################################################################
template<class T> template<class T_TENSOR> bool TENSOR_TESTS<T>::
Assert_Zero_Tensor(const T_TENSOR& a,const char* str,T tolerance) const
{
    bool b=true;
    T mx=0;
    for(int i=0;i<a.x.m;i++) if(a.x(i).Max_Abs()>tolerance){b=false;mx=std::max(mx,a.x(i).Max_Abs());}
    if(!b){
        status=false;
        LOG::cout<<"FAIL: "<<str<<"  "<<mx<<"  "<<tolerance<<std::endl;}
    return b;
}
//#####################################################################
// Function Assert_Equal
//#####################################################################
template<class T> template<class T_TENSOR0,class T_TENSOR1> bool TENSOR_TESTS<T>::
Assert_Equal_Tensor(const T_TENSOR0& a,const T_TENSOR1& b,const char* str,T tolerance) const
{
    bool k=Assert_Zero_Tensor(a-b,str,tolerance);
    return k;
}
template<class T> void TENSOR_TESTS<T>::
Test_0()
{
    Test_Contract(Fill_pt());
    Test_Contract_T<0,1>(Fill_pt());
    Test_Contract_T<1,0>(Fill_pt());
    Test_Contract_T<0,2>(Fill_pt());
    Test_Contract_T<2,0>(Fill_pt());
    Test_Contract_T<1,2>(Fill_pt());
    Test_Contract_T<2,1>(Fill_pt());

    Test_Transposed<0,1>(Fill_pt());
    Test_Transposed<1,0>(Fill_pt());
    Test_Transposed<0,2>(Fill_pt());
    Test_Transposed<2,0>(Fill_pt());
    Test_Transposed<1,2>(Fill_pt());
    Test_Transposed<2,1>(Fill_pt());

    Test_Twice_Symmetric_Part<0,1>(Fill_pt());
    Test_Twice_Symmetric_Part<1,0>(Fill_pt());
    Test_Twice_Symmetric_Part<0,2>(Fill_pt());
    Test_Twice_Symmetric_Part<2,0>(Fill_pt());
    Test_Twice_Symmetric_Part<1,2>(Fill_pt());
    Test_Twice_Symmetric_Part<2,1>(Fill_pt());

    Test_Contract_TM<0,0>(Fill_pt(),Fill_sm<3>());
    Test_Contract_TM<1,0>(Fill_pt(),Fill_sm<3>());
    Test_Contract_TM<2,0>(Fill_pt(),Fill_sm<3>());
    Test_Contract_TM<0,1>(Fill_pt(),Fill_sm<3>());
    Test_Contract_TM<1,1>(Fill_pt(),Fill_sm<3>());
    Test_Contract_TM<2,1>(Fill_pt(),Fill_sm<3>());

    Test_Contract_TM<0,0>(Fill_pt(),Fill_cm<3>());
    Test_Contract_TM<1,0>(Fill_pt(),Fill_cm<3>());
    Test_Contract_TM<2,0>(Fill_pt(),Fill_cm<3>());
    Test_Contract_TM<0,1>(Fill_pt(),Fill_cm<3>());
    Test_Contract_TM<1,1>(Fill_pt(),Fill_cm<3>());
    Test_Contract_TM<2,1>(Fill_pt(),Fill_cm<3>());

    Test_Contract_TM<0,0>(Fill_pt(),Fill_im<3>());
    Test_Contract_TM<1,0>(Fill_pt(),Fill_im<3>());
    Test_Contract_TM<2,0>(Fill_pt(),Fill_im<3>());
    Test_Contract_TM<0,1>(Fill_pt(),Fill_im<3>());
    Test_Contract_TM<1,1>(Fill_pt(),Fill_im<3>());
    Test_Contract_TM<2,1>(Fill_pt(),Fill_im<3>());

    Test_Add_Sub(Fill_pt(),Fill_gt<3,3,3>());
    Test_Add_Sub(Fill_pt(),Fill_zt<3,3,3>());
    Test_Add_Sub(Fill_pt(),Fill_pt());
}
template<class T> template<int m> void TENSOR_TESTS<T>::
Test_1()
{
    Test_Contract_TM<0,0>(Fill_pt(),Fill_gm<3,m>());
    Test_Contract_TM<1,0>(Fill_pt(),Fill_gm<3,m>());
    Test_Contract_TM<2,0>(Fill_pt(),Fill_gm<3,m>());
    Test_Contract_TM<0,1>(Fill_pt(),Fill_gm<m,3>());
    Test_Contract_TM<1,1>(Fill_pt(),Fill_gm<m,3>());
    Test_Contract_TM<2,1>(Fill_pt(),Fill_gm<m,3>());

    Test_Contract_TM<0,0>(Fill_pt(),Fill_zm<3,m>());
    Test_Contract_TM<1,0>(Fill_pt(),Fill_zm<3,m>());
    Test_Contract_TM<2,0>(Fill_pt(),Fill_zm<3,m>());
    Test_Contract_TM<0,1>(Fill_pt(),Fill_zm<m,3>());
    Test_Contract_TM<1,1>(Fill_pt(),Fill_zm<m,3>());
    Test_Contract_TM<2,1>(Fill_pt(),Fill_zm<m,3>());

    Test_Add_Sub(Fill_pt(),Fill_st<m-1,3,3>());
    Test_Add_Sub(Fill_pt(),Fill_vt<m-1,3,3>());
    Test_Add_Sub(Fill_pt(),Fill_it<m-1,3>());
}
template<class T> template<int m,int n> void TENSOR_TESTS<T>::
Test_2()
{
    Test_Contract(Fill_it<m-1,n>());
    Test_Contract_T<0,1>(Fill_gt<m,m,n>());
    Test_Contract_T<1,0>(Fill_gt<m,m,n>());
    Test_Contract_T<0,2>(Fill_gt<m,n,m>());
    Test_Contract_T<2,0>(Fill_gt<m,n,m>());
    Test_Contract_T<1,2>(Fill_gt<m,n,n>());
    Test_Contract_T<2,1>(Fill_gt<m,n,n>());

    Test_Contract_T<0,1>(Fill_zt<m,m,n>());
    Test_Contract_T<1,0>(Fill_zt<m,m,n>());
    Test_Contract_T<0,2>(Fill_zt<m,n,m>());
    Test_Contract_T<2,0>(Fill_zt<m,n,m>());
    Test_Contract_T<1,2>(Fill_zt<m,n,n>());
    Test_Contract_T<2,1>(Fill_zt<m,n,n>());

    Test_Contract_T<0,1>(Fill_st<n-1,m,m>());
    Test_Contract_T<1,0>(Fill_st<n-1,m,m>());
    Test_Contract_T<0,2>(Fill_st<n-1,m,m>());
    Test_Contract_T<2,0>(Fill_st<n-1,m,m>());
    Test_Contract_T<1,2>(Fill_st<n-1,m,m>());
    Test_Contract_T<2,1>(Fill_st<n-1,m,m>());

    Test_Contract_T<0,1>(Fill_vt<n-1,m,m>());
    Test_Contract_T<1,0>(Fill_vt<n-1,m,m>());
    Test_Contract_T<0,2>(Fill_vt<n-1,m,m>());
    Test_Contract_T<2,0>(Fill_vt<n-1,m,m>());
    Test_Contract_T<1,2>(Fill_vt<n-1,m,m>());
    Test_Contract_T<2,1>(Fill_vt<n-1,m,m>());

    Test_Contract_T<0,1>(Fill_it<n-1,m>());
    Test_Contract_T<1,0>(Fill_it<n-1,m>());
    Test_Contract_T<0,2>(Fill_it<n-1,m>());
    Test_Contract_T<2,0>(Fill_it<n-1,m>());
    Test_Contract_T<1,2>(Fill_it<n-1,m>());
    Test_Contract_T<2,1>(Fill_it<n-1,m>());

    Test_Transposed<0,1>(Fill_it<m-1,n>());
    Test_Transposed<1,0>(Fill_it<m-1,n>());
    Test_Transposed<0,2>(Fill_it<m-1,n>());
    Test_Transposed<2,0>(Fill_it<m-1,n>());
    Test_Transposed<1,2>(Fill_it<m-1,n>());
    Test_Transposed<2,1>(Fill_it<m-1,n>());

    Test_Twice_Symmetric_Part<m-1,m%3>(Fill_st<m-1,m,m>());
    Test_Twice_Symmetric_Part<m-1,m%3>(Fill_vt<m-1,m,m>());

    Test_Twice_Symmetric_Part<m-1,(m+1)%3>(Fill_st<m-1,m,m>());
    Test_Twice_Symmetric_Part<m-1,(m+1)%3>(Fill_vt<m-1,m,m>());

    Test_Twice_Symmetric_Part<0,1>(Fill_gt<m,m,n>());
    Test_Twice_Symmetric_Part<1,0>(Fill_gt<m,m,n>());
    Test_Twice_Symmetric_Part<0,2>(Fill_gt<m,n,m>());
    Test_Twice_Symmetric_Part<2,0>(Fill_gt<m,n,m>());
    Test_Twice_Symmetric_Part<1,2>(Fill_gt<m,n,n>());
    Test_Twice_Symmetric_Part<2,1>(Fill_gt<m,n,n>());

    Test_Twice_Symmetric_Part<0,1>(Fill_zt<m,m,n>());
    Test_Twice_Symmetric_Part<1,0>(Fill_zt<m,m,n>());
    Test_Twice_Symmetric_Part<0,2>(Fill_zt<m,n,m>());
    Test_Twice_Symmetric_Part<2,0>(Fill_zt<m,n,m>());
    Test_Twice_Symmetric_Part<1,2>(Fill_zt<m,n,n>());
    Test_Twice_Symmetric_Part<2,1>(Fill_zt<m,n,n>());

    Test_Twice_Symmetric_Part<0,1>(Fill_st<2,m,n>());
    Test_Twice_Symmetric_Part<1,0>(Fill_st<2,m,n>());
    Test_Twice_Symmetric_Part<0,2>(Fill_st<1,m,n>());
    Test_Twice_Symmetric_Part<2,0>(Fill_st<1,m,n>());
    Test_Twice_Symmetric_Part<1,2>(Fill_st<0,m,n>());
    Test_Twice_Symmetric_Part<2,1>(Fill_st<0,m,n>());

    Test_Twice_Symmetric_Part<0,1>(Fill_vt<2,m,n>());
    Test_Twice_Symmetric_Part<1,0>(Fill_vt<2,m,n>());
    Test_Twice_Symmetric_Part<0,2>(Fill_vt<1,m,n>());
    Test_Twice_Symmetric_Part<2,0>(Fill_vt<1,m,n>());
    Test_Twice_Symmetric_Part<1,2>(Fill_vt<0,m,n>());
    Test_Twice_Symmetric_Part<2,1>(Fill_vt<0,m,n>());

    Test_Twice_Symmetric_Part<0,1>(Fill_it<m-1,n>());
    Test_Twice_Symmetric_Part<1,0>(Fill_it<m-1,n>());
    Test_Twice_Symmetric_Part<0,2>(Fill_it<m-1,n>());
    Test_Twice_Symmetric_Part<2,0>(Fill_it<m-1,n>());
    Test_Twice_Symmetric_Part<1,2>(Fill_it<m-1,n>());
    Test_Twice_Symmetric_Part<2,1>(Fill_it<m-1,n>());

    Test_Add_Sub(Fill_st<0,m,n>(),Fill_gt<m,n,n>());
    Test_Add_Sub(Fill_st<1,m,n>(),Fill_gt<n,m,n>());
    Test_Add_Sub(Fill_st<2,m,n>(),Fill_gt<n,n,m>());
    Test_Add_Sub(Fill_st<0,m,n>(),Fill_zt<m,n,n>());
    Test_Add_Sub(Fill_st<1,m,n>(),Fill_zt<n,m,n>());
    Test_Add_Sub(Fill_st<2,m,n>(),Fill_zt<n,n,m>());

    Test_Add_Sub(Fill_vt<0,m,n>(),Fill_gt<m,n,n>());
    Test_Add_Sub(Fill_vt<1,m,n>(),Fill_gt<n,m,n>());
    Test_Add_Sub(Fill_vt<2,m,n>(),Fill_gt<n,n,m>());
    Test_Add_Sub(Fill_vt<0,m,n>(),Fill_zt<m,n,n>());
    Test_Add_Sub(Fill_vt<1,m,n>(),Fill_zt<n,m,n>());
    Test_Add_Sub(Fill_vt<2,m,n>(),Fill_zt<n,n,m>());

    Test_Add_Sub(Fill_it<m-1,n>(),Fill_gt<n,n,n>());
    Test_Add_Sub(Fill_it<m-1,n>(),Fill_zt<n,n,n>());
}
template<class T> template<int m,int n,int p> void TENSOR_TESTS<T>::
Test_3()
{
    Test_Contract(Fill_gt<m,n,p>());
    Test_Contract(Fill_st<m-1,n,p>());
    Test_Contract(Fill_zt<m,n,p>());
    Test_Contract(Fill_vt<m-1,n,p>());
    Test_Contract_T<(m-1)==0?1:0,(m-1)==2?1:2>(Fill_st<m-1,n,p>());
    Test_Contract_T<(m-1)==0?1:0,(m-1)==2?1:2>(Fill_vt<m-1,n,p>());
    Test_Tensor_Product<m-1>(Fill_sm<n>(),Fill_gv<p>());
    Test_Tensor_Product<m-1>(Fill_sm<n>(),Fill_zv<p>());
    Test_Tensor_Product<m-1>(Fill_cm<n>(),Fill_gv<p>());
    Test_Tensor_Product<m-1>(Fill_cm<n>(),Fill_zv<p>());
    Test_Tensor_Product<m-1>(Fill_im<n>(),Fill_gv<p>());
    Test_Tensor_Product<m-1>(Fill_im<n>(),Fill_zv<p>());
    Test_Transposed<0,1>(Fill_gt<m,n,p>());
    Test_Transposed<1,0>(Fill_gt<m,n,p>());
    Test_Transposed<0,2>(Fill_gt<m,n,p>());
    Test_Transposed<2,0>(Fill_gt<m,n,p>());
    Test_Transposed<1,2>(Fill_gt<m,n,p>());
    Test_Transposed<2,1>(Fill_gt<m,n,p>());
    Test_Transposed<0,1>(Fill_st<m-1,n,p>());
    Test_Transposed<1,0>(Fill_st<m-1,n,p>());
    Test_Transposed<0,2>(Fill_st<m-1,n,p>());
    Test_Transposed<2,0>(Fill_st<m-1,n,p>());
    Test_Transposed<1,2>(Fill_st<m-1,n,p>());
    Test_Transposed<2,1>(Fill_st<m-1,n,p>());
    Test_Transposed<0,1>(Fill_zt<m,n,p>());
    Test_Transposed<1,0>(Fill_zt<m,n,p>());
    Test_Transposed<0,2>(Fill_zt<m,n,p>());
    Test_Transposed<2,0>(Fill_zt<m,n,p>());
    Test_Transposed<1,2>(Fill_zt<m,n,p>());
    Test_Transposed<2,1>(Fill_zt<m,n,p>());
    Test_Transposed<0,1>(Fill_vt<m-1,n,p>());
    Test_Transposed<1,0>(Fill_vt<m-1,n,p>());
    Test_Transposed<0,2>(Fill_vt<m-1,n,p>());
    Test_Transposed<2,0>(Fill_vt<m-1,n,p>());
    Test_Transposed<1,2>(Fill_vt<m-1,n,p>());
    Test_Transposed<2,1>(Fill_vt<m-1,n,p>());

    Test_Contract_TM<0,0>(Fill_gt<m,n,p>(),Fill_cm<m>());
    Test_Contract_TM<1,0>(Fill_gt<m,n,p>(),Fill_cm<n>());
    Test_Contract_TM<2,0>(Fill_gt<m,n,p>(),Fill_cm<p>());
    Test_Contract_TM<0,1>(Fill_gt<m,n,p>(),Fill_cm<m>());
    Test_Contract_TM<1,1>(Fill_gt<m,n,p>(),Fill_cm<n>());
    Test_Contract_TM<2,1>(Fill_gt<m,n,p>(),Fill_cm<p>());

    Test_Contract_TM<0,0>(Fill_zt<m,n,p>(),Fill_cm<m>());
    Test_Contract_TM<1,0>(Fill_zt<m,n,p>(),Fill_cm<n>());
    Test_Contract_TM<2,0>(Fill_zt<m,n,p>(),Fill_cm<p>());
    Test_Contract_TM<0,1>(Fill_zt<m,n,p>(),Fill_cm<m>());
    Test_Contract_TM<1,1>(Fill_zt<m,n,p>(),Fill_cm<n>());
    Test_Contract_TM<2,1>(Fill_zt<m,n,p>(),Fill_cm<p>());

    Test_Contract_TM<m-1,0>(Fill_st<m-1,n,p>(),Fill_cm<n>());
    Test_Contract_TM<m==1?1:0,0>(Fill_st<m-1,n,p>(),Fill_cm<p>());
    Test_Contract_TM<m==3?1:2,0>(Fill_st<m-1,n,p>(),Fill_cm<p>());
    Test_Contract_TM<m-1,1>(Fill_st<m-1,n,p>(),Fill_cm<n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_st<m-1,n,p>(),Fill_cm<p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_st<m-1,n,p>(),Fill_cm<p>());

    Test_Contract_TM<m-1,0>(Fill_vt<m-1,n,p>(),Fill_cm<n>());
    Test_Contract_TM<m==1?1:0,0>(Fill_vt<m-1,n,p>(),Fill_cm<p>());
    Test_Contract_TM<m==3?1:2,0>(Fill_vt<m-1,n,p>(),Fill_cm<p>());
    Test_Contract_TM<m-1,1>(Fill_vt<m-1,n,p>(),Fill_cm<n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_vt<m-1,n,p>(),Fill_cm<p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_vt<m-1,n,p>(),Fill_cm<p>());

    Test_Contract_TM<p-1,0>(Fill_it<m-1,n>(),Fill_cm<n>());
    Test_Contract_TM<p-1,1>(Fill_it<m-1,n>(),Fill_cm<n>());

    Test_Contract_TM<0,0>(Fill_gt<m,n,p>(),Fill_im<m>());
    Test_Contract_TM<1,0>(Fill_gt<m,n,p>(),Fill_im<n>());
    Test_Contract_TM<2,0>(Fill_gt<m,n,p>(),Fill_im<p>());
    Test_Contract_TM<0,1>(Fill_gt<m,n,p>(),Fill_im<m>());
    Test_Contract_TM<1,1>(Fill_gt<m,n,p>(),Fill_im<n>());
    Test_Contract_TM<2,1>(Fill_gt<m,n,p>(),Fill_im<p>());

    Test_Contract_TM<0,0>(Fill_zt<m,n,p>(),Fill_im<m>());
    Test_Contract_TM<1,0>(Fill_zt<m,n,p>(),Fill_im<n>());
    Test_Contract_TM<2,0>(Fill_zt<m,n,p>(),Fill_im<p>());
    Test_Contract_TM<0,1>(Fill_zt<m,n,p>(),Fill_im<m>());
    Test_Contract_TM<1,1>(Fill_zt<m,n,p>(),Fill_im<n>());
    Test_Contract_TM<2,1>(Fill_zt<m,n,p>(),Fill_im<p>());

    Test_Contract_TM<m-1,0>(Fill_st<m-1,n,p>(),Fill_im<n>());
    Test_Contract_TM<m==1?1:0,0>(Fill_st<m-1,n,p>(),Fill_im<p>());
    Test_Contract_TM<m==3?1:2,0>(Fill_st<m-1,n,p>(),Fill_im<p>());
    Test_Contract_TM<m-1,1>(Fill_st<m-1,n,p>(),Fill_im<n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_st<m-1,n,p>(),Fill_im<p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_st<m-1,n,p>(),Fill_im<p>());

    Test_Contract_TM<m-1,0>(Fill_vt<m-1,n,p>(),Fill_im<n>());
    Test_Contract_TM<m==1?1:0,0>(Fill_vt<m-1,n,p>(),Fill_im<p>());
    Test_Contract_TM<m==3?1:2,0>(Fill_vt<m-1,n,p>(),Fill_im<p>());
    Test_Contract_TM<m-1,1>(Fill_vt<m-1,n,p>(),Fill_im<n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_vt<m-1,n,p>(),Fill_im<p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_vt<m-1,n,p>(),Fill_im<p>());

    Test_Contract_TM<p-1,0>(Fill_it<m-1,n>(),Fill_im<n>());
    Test_Contract_TM<p-1,1>(Fill_it<m-1,n>(),Fill_im<n>());

    Test_Contract_TM<0,0>(Fill_gt<m,n,p>(),Fill_sm<m>());
    Test_Contract_TM<1,0>(Fill_gt<m,n,p>(),Fill_sm<n>());
    Test_Contract_TM<2,0>(Fill_gt<m,n,p>(),Fill_sm<p>());
    Test_Contract_TM<0,1>(Fill_gt<m,n,p>(),Fill_sm<m>());
    Test_Contract_TM<1,1>(Fill_gt<m,n,p>(),Fill_sm<n>());
    Test_Contract_TM<2,1>(Fill_gt<m,n,p>(),Fill_sm<p>());

    Test_Contract_TM<0,0>(Fill_zt<m,n,p>(),Fill_sm<m>());
    Test_Contract_TM<1,0>(Fill_zt<m,n,p>(),Fill_sm<n>());
    Test_Contract_TM<2,0>(Fill_zt<m,n,p>(),Fill_sm<p>());
    Test_Contract_TM<0,1>(Fill_zt<m,n,p>(),Fill_sm<m>());
    Test_Contract_TM<1,1>(Fill_zt<m,n,p>(),Fill_sm<n>());
    Test_Contract_TM<2,1>(Fill_zt<m,n,p>(),Fill_sm<p>());

    Test_Contract_TM<m-1,0>(Fill_st<m-1,n,p>(),Fill_sm<n>());
    Test_Contract_TM<m==1?1:0,0>(Fill_st<m-1,n,p>(),Fill_sm<p>());
    Test_Contract_TM<m==3?1:2,0>(Fill_st<m-1,n,p>(),Fill_sm<p>());
    Test_Contract_TM<m-1,1>(Fill_st<m-1,n,p>(),Fill_sm<n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_st<m-1,n,p>(),Fill_sm<p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_st<m-1,n,p>(),Fill_sm<p>());

    Test_Contract_TM<m-1,0>(Fill_vt<m-1,n,p>(),Fill_sm<n>());
    Test_Contract_TM<m==1?1:0,0>(Fill_vt<m-1,n,p>(),Fill_sm<p>());
    Test_Contract_TM<m==3?1:2,0>(Fill_vt<m-1,n,p>(),Fill_sm<p>());
    Test_Contract_TM<m-1,1>(Fill_vt<m-1,n,p>(),Fill_sm<n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_vt<m-1,n,p>(),Fill_sm<p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_vt<m-1,n,p>(),Fill_sm<p>());

    Test_Contract_TM<p-1,0>(Fill_it<m-1,n>(),Fill_sm<n>());
    Test_Contract_TM<p-1,1>(Fill_it<m-1,n>(),Fill_sm<n>());

    Test_Add_Sub(Fill_gt<m,n,p>(),Fill_gt<m,n,p>());
    Test_Add_Sub(Fill_gt<m,n,p>(),Fill_zt<m,n,p>());
    Test_Add_Sub(Fill_zt<m,n,p>(),Fill_zt<m,n,p>());
    
    Test_Add_Sub(Fill_st<m-1,n,p>(),Fill_st<m-1,n,p>());
    Test_Add_Sub(Fill_st<m-1,n,p>(),Fill_vt<m-1,n,p>());
    Test_Add_Sub(Fill_vt<m-1,n,p>(),Fill_vt<m-1,n,p>());
                            
    Test_Add_Sub(Fill_it<m-1,n>(),Fill_it<p-1,n>());
    Test_Add_Sub(Fill_it<m-1,n>(),Fill_st<p-1,n,n>());
    Test_Add_Sub(Fill_it<m-1,n>(),Fill_vt<p-1,n,n>());
                                                    
    Test_Add_Sub(Fill_st<m-1,n,n>(),Fill_st<p-1,n,n>());
    Test_Add_Sub(Fill_st<m-1,n,n>(),Fill_vt<p-1,n,n>());
    Test_Add_Sub(Fill_vt<m-1,n,n>(),Fill_vt<p-1,n,n>());
}
template<class T> template<int m,int n,int p,int q> void TENSOR_TESTS<T>::
Test_4()
{
    Test_Tensor_Product<m-1>(Fill_gm<n,p>(),Fill_gv<q>());
    Test_Tensor_Product<m-1>(Fill_gm<n,p>(),Fill_zv<q>());
    Test_Tensor_Product<m-1>(Fill_zm<n,p>(),Fill_gv<q>());
    Test_Tensor_Product<m-1>(Fill_zm<n,p>(),Fill_zv<q>());

    Test_Contract_TM<0,0>(Fill_gt<m,n,p>(),Fill_gm<m,q>());
    Test_Contract_TM<1,0>(Fill_gt<m,n,p>(),Fill_gm<n,q>());
    Test_Contract_TM<2,0>(Fill_gt<m,n,p>(),Fill_gm<p,q>());
    Test_Contract_TM<0,1>(Fill_gt<m,n,p>(),Fill_gm<q,m>());
    Test_Contract_TM<1,1>(Fill_gt<m,n,p>(),Fill_gm<q,n>());
    Test_Contract_TM<2,1>(Fill_gt<m,n,p>(),Fill_gm<q,p>());

    Test_Contract_TM<0,0>(Fill_zt<m,n,p>(),Fill_gm<m,q>());
    Test_Contract_TM<1,0>(Fill_zt<m,n,p>(),Fill_gm<n,q>());
    Test_Contract_TM<2,0>(Fill_zt<m,n,p>(),Fill_gm<p,q>());
    Test_Contract_TM<0,1>(Fill_zt<m,n,p>(),Fill_gm<q,m>());
    Test_Contract_TM<1,1>(Fill_zt<m,n,p>(),Fill_gm<q,n>());
    Test_Contract_TM<2,1>(Fill_zt<m,n,p>(),Fill_gm<q,p>());

    Test_Contract_TM<m-1,0>(Fill_st<m-1,n,p>(),Fill_gm<n,q>());
    Test_Contract_TM<m==1?1:0,0>(Fill_st<m-1,n,p>(),Fill_gm<p,q>());
    Test_Contract_TM<m==3?1:2,0>(Fill_st<m-1,n,p>(),Fill_gm<p,q>());
    Test_Contract_TM<m-1,1>(Fill_st<m-1,n,p>(),Fill_gm<q,n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_st<m-1,n,p>(),Fill_gm<q,p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_st<m-1,n,p>(),Fill_gm<q,p>());

    Test_Contract_TM<m-1,0>(Fill_vt<m-1,n,p>(),Fill_gm<n,q>());
    Test_Contract_TM<m==1?1:0,0>(Fill_vt<m-1,n,p>(),Fill_gm<p,q>());
    Test_Contract_TM<m==3?1:2,0>(Fill_vt<m-1,n,p>(),Fill_gm<p,q>());
    Test_Contract_TM<m-1,1>(Fill_vt<m-1,n,p>(),Fill_gm<q,n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_vt<m-1,n,p>(),Fill_gm<q,p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_vt<m-1,n,p>(),Fill_gm<q,p>());

    //////Test_Contract_TM<p-1,0>(Fill_it<m-1,n>(),Fill_gm<n,q>());
    //////Test_Contract_TM<p-1,1>(Fill_it<m-1,n>(),Fill_gm<q,n>());

    Test_Contract_TM<0,0>(Fill_gt<m,n,p>(),Fill_zm<m,q>());
    Test_Contract_TM<1,0>(Fill_gt<m,n,p>(),Fill_zm<n,q>());
    Test_Contract_TM<2,0>(Fill_gt<m,n,p>(),Fill_zm<p,q>());
    Test_Contract_TM<0,1>(Fill_gt<m,n,p>(),Fill_zm<q,m>());
    Test_Contract_TM<1,1>(Fill_gt<m,n,p>(),Fill_zm<q,n>());
    Test_Contract_TM<2,1>(Fill_gt<m,n,p>(),Fill_zm<q,p>());
    
    Test_Contract_TM<0,0>(Fill_zt<m,n,p>(),Fill_zm<m,q>());
    Test_Contract_TM<1,0>(Fill_zt<m,n,p>(),Fill_zm<n,q>());
    Test_Contract_TM<2,0>(Fill_zt<m,n,p>(),Fill_zm<p,q>());
    Test_Contract_TM<0,1>(Fill_zt<m,n,p>(),Fill_zm<q,m>());
    Test_Contract_TM<1,1>(Fill_zt<m,n,p>(),Fill_zm<q,n>());
    Test_Contract_TM<2,1>(Fill_zt<m,n,p>(),Fill_zm<q,p>());

    Test_Contract_TM<m-1,0>(Fill_st<m-1,n,p>(),Fill_zm<n,q>());
    Test_Contract_TM<m==1?1:0,0>(Fill_st<m-1,n,p>(),Fill_zm<p,q>());
    Test_Contract_TM<m==3?1:2,0>(Fill_st<m-1,n,p>(),Fill_zm<p,q>());
    Test_Contract_TM<m-1,1>(Fill_st<m-1,n,p>(),Fill_zm<q,n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_st<m-1,n,p>(),Fill_zm<q,p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_st<m-1,n,p>(),Fill_zm<q,p>());

    Test_Contract_TM<m-1,0>(Fill_vt<m-1,n,p>(),Fill_zm<n,q>());
    Test_Contract_TM<m==1?1:0,0>(Fill_vt<m-1,n,p>(),Fill_zm<p,q>());
    Test_Contract_TM<m==3?1:2,0>(Fill_vt<m-1,n,p>(),Fill_zm<p,q>());
    Test_Contract_TM<m-1,1>(Fill_vt<m-1,n,p>(),Fill_zm<q,n>());
    Test_Contract_TM<m==1?1:0,1>(Fill_vt<m-1,n,p>(),Fill_zm<q,p>());
    Test_Contract_TM<m==3?1:2,1>(Fill_vt<m-1,n,p>(),Fill_zm<q,p>());

    Test_Contract_TM<p-1,0>(Fill_it<m-1,n>(),Fill_zm<n,q>());
    Test_Contract_TM<p-1,1>(Fill_it<m-1,n>(),Fill_zm<q,n>());
}
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> template<class T_TENSOR> void TENSOR_TESTS<T>::
Test_Contract(const T_TENSOR& t)
{
    Test_Contract_TV<0>(t,Fill_gv<T_TENSOR::m>());
    Test_Contract_TV<1>(t,Fill_gv<T_TENSOR::n>());
    Test_Contract_TV<2>(t,Fill_gv<T_TENSOR::p>());
    Test_Contract_TV<0>(t,Fill_zv<T_TENSOR::m>());
    Test_Contract_TV<1>(t,Fill_zv<T_TENSOR::n>());
    Test_Contract_TV<2>(t,Fill_zv<T_TENSOR::p>());
}
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> template<int s,class T_TENSOR,class T_VECTOR> void TENSOR_TESTS<T>::
Test_Contract_TV(const T_TENSOR& t,const T_VECTOR& v)
{
    TENSOR<T,T_TENSOR::m,T_TENSOR::n,T_TENSOR::p> gt;
    gt+=t;
    VECTOR<T,T_VECTOR::m> gv;
    gv+=v;
    auto u=Contract<s>(t,v);
    auto gu=Contract<s>(gt,gv);
    MATRIX<T,decltype(gu)::m,decltype(gu)::n> m0,m1;
    m0+=u;
    m1+=gu;
    if(!Assert_Equal_Matrix(m0,m1,"Contract TV",1e-10)) LOG::cout<<"TYPE: "<<s<<" "<<typeid(t).name()<<" "<<typeid(v).name()<<std::endl;
}
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> template<int r,int s,class T_TENSOR> void TENSOR_TESTS<T>::
Test_Contract_T(const T_TENSOR& t)
{
    TENSOR<T,T_TENSOR::m,T_TENSOR::n,T_TENSOR::p> gt;
    gt+=t;
    auto u=Contract<r,s>(t);
    auto gu=Contract<r,s>(gt);
    VECTOR<T,decltype(gu)::m> m0,m1;
    m0+=u;
    m1+=gu;
    if(!Assert_Equal_Vector(m0,m1,"Contract T",1e-10)) LOG::cout<<"TYPE: "<<r<<" "<<s<<" "<<typeid(t).name()<<std::endl;
}
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> template<int r,int s,class T_TENSOR,class T_MATRIX> void TENSOR_TESTS<T>::
Test_Contract_TM(const T_TENSOR& t,const T_MATRIX& m)
{
    TENSOR<T,T_TENSOR::m,T_TENSOR::n,T_TENSOR::p> gt;
    gt+=t;
    MATRIX<T,T_MATRIX::m,T_MATRIX::n> gm;
    gm+=m;
    auto u=Contract<r,s>(t,m);
    auto gu=Contract<r,s>(gt,m);
    TENSOR<T,decltype(gu)::m,decltype(gu)::n,decltype(gu)::p> m0,m1;
    m0+=u;
    m1+=gu;
    if(!Assert_Equal_Tensor(m0,m1,"Contract TM",1e-10)) LOG::cout<<"TYPE: "<<r<<" "<<s<<" "<<typeid(t).name()<<" "<<typeid(m).name()<<std::endl;
}
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> template<int s,class T_MATRIX,class T_VECTOR> void TENSOR_TESTS<T>::
Test_Tensor_Product(const T_MATRIX& m,const T_VECTOR& v)
{
    MATRIX<T,T_MATRIX::m,T_MATRIX::n> gm;
    gm+=m;
    VECTOR<T,T_VECTOR::m> gv;
    gv+=v;
    auto u=Tensor_Product<s>(m,v);
    auto gu=Tensor_Product<s>(gm,gv);
    TENSOR<T,decltype(gu)::m,decltype(gu)::n,decltype(gu)::p> m0,m1;
    m0+=u;
    m1+=gu;
    if(!Assert_Equal_Tensor(m0,m1,"Tensor Product",1e-10)) LOG::cout<<"TYPE: "<<s<<" "<<typeid(m).name()<<" "<<typeid(v).name()<<std::endl;
}
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> template<int r,int s,class T_TENSOR> void TENSOR_TESTS<T>::
Test_Transposed(const T_TENSOR& m)
{
    TENSOR<T,T_TENSOR::m,T_TENSOR::n,T_TENSOR::p> gm;
    gm+=m;
    auto u=Transposed<r,s>(m);
    auto gu=Transposed<r,s>(gm);
    TENSOR<T,decltype(gu)::m,decltype(gu)::n,decltype(gu)::p> m0,m1;
    m0+=u;
    m1+=gu;
    if(!Assert_Equal_Tensor(m0,m1,"Transposed",1e-10)) LOG::cout<<"TYPE: "<<s<<" "<<typeid(m).name()<<std::endl;
}
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> template<int r,int s,class T_TENSOR> void TENSOR_TESTS<T>::
Test_Twice_Symmetric_Part(const T_TENSOR& m)
{
    TENSOR<T,T_TENSOR::m,T_TENSOR::n,T_TENSOR::p> gm;
    gm+=m;
    auto u=Twice_Symmetric_Part<r,s>(m);
    auto gu=Twice_Symmetric_Part<r,s>(gm);
    SYMMETRIC_TENSOR<T,3-r-s,decltype(gu)::um,decltype(gu)::un> m0,m1;
    m0+=u;
    m1+=gu;
    if(!Assert_Equal_Tensor(m0,m1,"Twice_Symmetric_Part",1e-10)) LOG::cout<<"TYPE: "<<s<<" "<<typeid(m).name()<<std::endl;
}
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> template<class T_TENSOR0,class T_TENSOR1> void TENSOR_TESTS<T>::
Test_Add_Sub(const T_TENSOR0& a,const T_TENSOR1& b)
{
    {
        TENSOR<T,T_TENSOR0::m,T_TENSOR0::n,T_TENSOR0::p> ga,gb,m0,m1;
        ga+=a;
        gb+=b;
        auto u=a+b;
        auto gu=ga+gb;
        m0+=u;
        m1+=gu;
        if(!Assert_Equal_Tensor(m0,m1,"a+b",1e-10)) LOG::cout<<"TYPE: "<<typeid(a).name()<<" "<<typeid(b).name()<<std::endl;
    }
    {
        TENSOR<T,T_TENSOR0::m,T_TENSOR0::n,T_TENSOR0::p> ga,gb,m0,m1;
        ga+=a;
        gb+=b;
        auto u=a-b;
        auto gu=ga-gb;
        m0+=u;
        m1+=gu;
        if(!Assert_Equal_Tensor(m0,m1,"a-b",1e-10)) LOG::cout<<"TYPE: "<<typeid(a).name()<<" "<<typeid(b).name()<<std::endl;
    }
    {
        TENSOR<T,T_TENSOR0::m,T_TENSOR0::n,T_TENSOR0::p> ga,gb,m0,m1;
        ga+=a;
        gb+=b;
        auto u=b+a;
        auto gu=gb+ga;
        m0+=u;
        m1+=gu;
        if(!Assert_Equal_Tensor(m0,m1,"b+a",1e-10)) LOG::cout<<"TYPE: "<<typeid(a).name()<<" "<<typeid(b).name()<<std::endl;
    }
    {
        TENSOR<T,T_TENSOR0::m,T_TENSOR0::n,T_TENSOR0::p> ga,gb,m0,m1;
        ga+=a;
        gb+=b;
        auto u=b-a;
        auto gu=gb-ga;
        m0+=u;
        m1+=gu;
        if(!Assert_Equal_Tensor(m0,m1,"b-a",1e-10)) LOG::cout<<"TYPE: "<<typeid(a).name()<<" "<<typeid(b).name()<<std::endl;
    }
}
//#####################################################################
// Function Fill
//#####################################################################
template<class T> template<int m,int n,int p> TENSOR<T,m,n,p> TENSOR_TESTS<T>::
Fill_gt()
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<t.x.m;i++) rand.Fill_Uniform(t.x(i),-1,1);
    return t;
}
//#####################################################################
// Function Fill
//#####################################################################
template<class T> template<int m,int n,int p> SYMMETRIC_TENSOR<T,m,n,p> TENSOR_TESTS<T>::
Fill_st()
{
    SYMMETRIC_TENSOR<T,m,n,p> t;
    for(int i=0;i<t.x.m;i++) rand.Fill_Uniform(t.x(i),-1,1);
    return t;
}
//#####################################################################
// Function Fill_zt
//#####################################################################
template<class T> template<int m,int n,int p> ZERO_TENSOR<T,m,n,p> TENSOR_TESTS<T>::
Fill_zt()
{
    return ZERO_TENSOR<T,m,n,p>();
}
//#####################################################################
// Function Fill_vt
//#####################################################################
template<class T> template<int u,int m,int n> VEC_ID_TENSOR<T,u,m,n> TENSOR_TESTS<T>::
Fill_vt()
{
    VEC_ID_TENSOR<T,u,m,n> v;
    rand.Fill_Uniform(v.v,-1,1);
    return v;
}
//#####################################################################
// Function Fill_it
//#####################################################################
template<class T> template<int u,int m> VEC_ID_SYM_TENSOR<T,u,m> TENSOR_TESTS<T>::
Fill_it()
{
    VEC_ID_SYM_TENSOR<T,u,m> v;
    rand.Fill_Uniform(v.v,-1,1);
    return v;
}
//#####################################################################
// Function Fill_pt
//#####################################################################
template<class T> PERMUTATION_TENSOR<T> TENSOR_TESTS<T>::
Fill_pt()
{
    PERMUTATION_TENSOR<T> v;
    rand.Fill_Uniform(v.x,-1,1);
    return v;
}

//#####################################################################
// Function Fill_gm
//#####################################################################
template<class T> template<int m,int n> MATRIX<T,m,n> TENSOR_TESTS<T>::
Fill_gm()
{
    MATRIX<T,m,n> v;
    rand.Fill_Uniform(v,-1,1);
    return v;
}
//#####################################################################
// Function Fill_sm
//#####################################################################
template<class T> template<int m> SYMMETRIC_MATRIX<T,m> TENSOR_TESTS<T>::
Fill_sm()
{
    SYMMETRIC_MATRIX<T,m> v;
    rand.Fill_Uniform(v,-1,1);
    return v;
}
//#####################################################################
// Function Fill_zm
//#####################################################################
template<class T> template<int m,int n> ZERO_MATRIX<T,m,n> TENSOR_TESTS<T>::
Fill_zm()
{
    return ZERO_MATRIX<T,m,n>();
}
//#####################################################################
// Function Fill_cm
//#####################################################################
template<class T> template<int m> SCALE_MATRIX<T,m> TENSOR_TESTS<T>::
Fill_cm()
{
    SCALE_MATRIX<T,m> v;
    rand.Fill_Uniform(v.x,-1,1);
    return v;
}
//#####################################################################
// Function Fill_im
//#####################################################################
template<class T> template<int m> IDENTITY_MATRIX<T,m> TENSOR_TESTS<T>::
Fill_im()
{
    return IDENTITY_MATRIX<T,m>();
}

//#####################################################################
// Function Fill_gv
//#####################################################################
template<class T> template<int m> VECTOR<T,m> TENSOR_TESTS<T>::
Fill_gv()
{
    VECTOR<T,m> v;
    rand.Fill_Uniform(v,-1,1);
    return v;
}
//#####################################################################
// Function Fill_zv
//#####################################################################
template<class T> template<int m> ZERO_VECTOR<T,m> TENSOR_TESTS<T>::
Fill_zv()
{
    return ZERO_VECTOR<T,m>();
}
#endif
