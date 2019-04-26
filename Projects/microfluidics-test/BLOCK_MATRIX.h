//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BLOCK_MATRIX__
#define __BLOCK_MATRIX__
#include <Core/Matrices/MATRIX_MXN.h>
#include "COMMON.h"

namespace PhysBAM{

template<class TV>
struct BLOCK_MATRIX
{
    typedef typename TV::SCALAR T;
    DOF_COUNTS nr,nc;
    MATRIX_MXN<T> M;

    void Resize()
    {
        M.Resize(TV::m*nr.v+TV::m*nr.e+nr.p,TV::m*nc.v+TV::m*nc.e+nc.p);
    }

    MATRIX<T,TV::m> Get_M(int r,int c) const
    {
        MATRIX<T,TV::m> m;
        M.Get_Submatrix(r,c,m);
        return m;
    }
    TV Get_C(int r,int c) const
    {
        TV u;
        for(int i=0;i<TV::m;i++) u(i)=M(r+i,c);
        return u;
    }
    TV Get_R(int r,int c) const
    {
        TV u;
        for(int i=0;i<TV::m;i++) u(i)=M(r,c+i);
        return u;
    }

    void Add_M(int r,int c,MATRIX<T,TV::m> m) {M.Add_To_Submatrix(r,c,m);}
    void Add_C(int r,int c,TV u) {for(int i=0;i<TV::m;i++) M(r+i,c)+=u(i);}
    void Add_R(int r,int c,TV u) {for(int i=0;i<TV::m;i++) M(r,c+i)+=u(i);}

    MATRIX<T,TV::m> Get_vv(int r,int c) const {return Get_M(TV::m*r,TV::m*c);}
    MATRIX<T,TV::m> Get_ve(int r,int c) const {return Get_M(TV::m*r,TV::m*(c+nc.v));}
    MATRIX<T,TV::m> Get_ev(int r,int c) const {return Get_M(TV::m*(r+nr.v),TV::m*c);}
    MATRIX<T,TV::m> Get_ee(int r,int c) const {return Get_M(TV::m*(r+nr.v),TV::m*(c+nc.v));}
    MATRIX<T,TV::m> Get_uu(int r,int er,int c,int ec) const {return Get_M(TV::m*(r+er*nr.v),TV::m*(c+ec*nc.v));}
    TV Get_vp(int r,int c) const {return Get_C(TV::m*r,TV::m*(nc.v+nc.e)+c);}
    TV Get_pv(int r,int c) const {return Get_R(TV::m*(nr.v+nr.e)+r,TV::m*c);}
    TV Get_ep(int r,int c) const {return Get_C(TV::m*(r+nr.v),TV::m*(nc.v+nc.e)+c);}
    TV Get_pe(int r,int c) const {return Get_R(TV::m*(nr.v+nr.e)+r,TV::m*(c+nc.v));}
    TV Get_up(int r,int er,int c) const {return Get_C(TV::m*(r+er*nr.v),TV::m*(nc.v+nc.e)+c);}
    TV Get_pu(int r,int c,int ec) const {return Get_R(TV::m*(nr.v+nr.e)+r,TV::m*(c+ec*nc.v));}

    void Add_vv(int r,int c,MATRIX<T,TV::m> m) {Add_M(TV::m*r,TV::m*c,m);}
    void Add_ve(int r,int c,MATRIX<T,TV::m> m) {Add_M(TV::m*r,TV::m*(c+nc.v),m);}
    void Add_ev(int r,int c,MATRIX<T,TV::m> m) {Add_M(TV::m*(r+nr.v),TV::m*c,m);}
    void Add_ee(int r,int c,MATRIX<T,TV::m> m) {Add_M(TV::m*(r+nr.v),TV::m*(c+nc.v),m);}
    void Add_uu(int r,int er,int c,int ec,MATRIX<T,TV::m> m) {Add_M(TV::m*(r+er*nr.v),TV::m*(c+ec*nc.v),m);}
    void Add_vp(int r,int c,TV u) {Add_C(TV::m*r,TV::m*(nc.v+nc.e)+c,u);}
    void Add_pv(int r,int c,TV u) {Add_R(TV::m*(nr.v+nr.e)+r,TV::m*c,u);}
    void Add_ep(int r,int c,TV u) {Add_C(TV::m*(r+nr.v),TV::m*(nc.v+nc.e)+c,u);}
    void Add_pe(int r,int c,TV u) {Add_R(TV::m*(nr.v+nr.e)+r,TV::m*(c+nc.v),u);}
    void Add_up(int r,int er,int c,TV u) {Add_C(TV::m*(r+er*nr.v),TV::m*(nc.v+nc.e)+c,u);}
    void Add_pu(int r,int c,int ec,TV u) {Add_R(TV::m*(nr.v+nr.e)+r,TV::m*(c+ec*nc.v),u);}
};

}
#endif
