//#####################################################################
// Copyright 2007, Jon Gretarsson, Avi Robinson-Mosher, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_ROW
//#####################################################################
#ifndef __SPARSE_MATRIX_ROW__
#define __SPARSE_MATRIX_ROW__

#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
namespace PhysBAM{

template<class T>
class SPARSE_MATRIX_ROW
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;
    int n;
    ARRAY_VIEW<const SPARSE_MATRIX_ENTRY<T> > A;

    SPARSE_MATRIX_ROW(const SPARSE_MATRIX_FLAT_MXN<T>& M,int r)
        :n(M.n),A(M.offsets(r+1)-M.offsets(r),&M.A(M.offsets(r)))
    {}

    T Dot(const SPARSE_MATRIX_ROW& r) const
    {
        assert(n==r.n);
        T x=0;
        int i=0,j=0;
        while(i<A.m && j<r.A.m){
            int a=A(i).j,b=r.A(j).j;
            if(a<b) i++;
            else if(b<a) j++;
            else x+=A(i++).a*r.A(j++).a;}
        return x;
    }

    T Dot(const ARRAY<T>& r) const
    {
        assert(n==r.m);
        T x=0;
        for(int i=0;i<A.m;i++)
            x+=A(i).a*r(A(i).j);
        return x;
    }

    static void Add_Multiple(T a,const SPARSE_MATRIX_ROW& v,ARRAY<T>& w)
    {
        for(int i=0;i<v.A.m;i++)
            w(v.A(i).j)+=a*v.A(i).a;
    }
};
}
#endif
