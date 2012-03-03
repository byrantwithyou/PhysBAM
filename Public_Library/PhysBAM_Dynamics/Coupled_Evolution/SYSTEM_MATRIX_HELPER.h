//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_MATRIX_HELPER
//#####################################################################
#ifndef __SYSTEM_MATRIX_HELPER__
#define __SYSTEM_MATRIX_HELPER__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{

template<class T> struct SYSTEM_MATRIX_HELPER;
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class T> class SPARSE_MATRIX_FLAT_NXN;

template<class T>
struct SYSTEM_MATRIX_BASE
{
    virtual void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const=0;
};

template<class T>
struct SYSTEM_MATRIX_HELPER:public NONCOPYABLE
{
    ARRAY<TRIPLE<int,int,T> > data;
    int start;
    bool compacted;

    SYSTEM_MATRIX_HELPER()
        :start(0), compacted(false)
    {}

    void New_Block()
    {start=data.m;compacted=false;}

    void Add_Matrix(const SYSTEM_MATRIX_BASE<T>& base,bool trans=false,int dr=0,int dc=0);
    void Add_Matrix(const SPARSE_MATRIX_FLAT_MXN<T>& M,bool trans=false,int dr=0,int dc=0);
    void Transpose();
    void Scale(T s);
    void Shift(int dr,int dc);
    void Compact(int rows, T tol=0);
    void Set_Matrix(int m,int n,SPARSE_MATRIX_FLAT_MXN<T>& M, T tol=0);
    void Set_Matrix(int n,SPARSE_MATRIX_FLAT_NXN<T>& M, T tol=0);

    static void Base_To_Matrix(int m,int n,const SYSTEM_MATRIX_BASE<T>& base,SPARSE_MATRIX_FLAT_MXN<T>& M,bool tranpose=false);
};
}
#endif
