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
}
#endif
