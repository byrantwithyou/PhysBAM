//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __XFORM__
#define __XFORM__
#include <Core/Matrices/MATRIX.h>
namespace PhysBAM{

template<class TV>
struct XFORM
{
    typedef typename TV::SCALAR T;
    MATRIX<T,TV::m> M;
    TV b;

    XFORM(const TV& b=TV()) :M(MATRIX<T,TV::m>::Identity_Matrix()),b(b) {}
    XFORM(const MATRIX<T,TV::m>& M,const TV& b) :M(M),b(b) {}

    XFORM operator*(const XFORM& x) const {return {M*x.M,M*x.b+b};}
    TV operator*(const TV& x) const {return M*x+b;}
    XFORM Inverse() const {auto inv=M.Inverse();return {inv,-inv*b};}
};

template<int d,class T>
MATRIX<T,d> To_Dim(const MATRIX<T,2>& M)
{
    MATRIX<T,d> r;
    r(d-1,d-1)=1;
    r.Set_Submatrix(0,0,M);
    return r;
}

template<class T>
VECTOR<T,2> xform(const XFORM<VECTOR<T,2> >& A,const VECTOR<T,2>& v) {return A*v;}

template<class T>
VECTOR<T,3> xform(const XFORM<VECTOR<T,2> >& A,const VECTOR<T,3>& v) {return (A*v.Remove_Index(2)).Append(v(2));}

template<class T>
VECTOR<T,2> xform(const MATRIX<T,2>& A,const VECTOR<T,2>& v) {return A*v;}

template<class T>
VECTOR<T,3> xform(const MATRIX<T,2>& A,const VECTOR<T,3>& v) {return (A*v.Remove_Index(2)).Append(v(2));}

template<class T>
VECTOR<T,2> xform_trans(const MATRIX<T,2>& A,const VECTOR<T,2>& v) {return A.Transpose_Times(v);}

template<class T>
VECTOR<T,3> xform_trans(const MATRIX<T,2>& A,const VECTOR<T,3>& v) {return (A.Transpose_Times(A*v.Remove_Index(2))).Append(v(2));}

}
#endif
