//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BLOCK_MATRIX__
#define __BLOCK_MATRIX__
#include <Core/Matrices/MATRIX_MXN.h>

namespace PhysBAM{

template<class T>
struct BLOCK_MATRIX
{
    typedef VECTOR<T,2> TV;
    int nv_r,ne_r,np_r,nv_c,ne_c,np_c;
    MATRIX_MXN<T> M;

    void Resize()
    {
        M.Resize(2*nv_r+2*ne_r+np_r,2*nv_c+2*ne_c+np_c);
    }

    MATRIX<T,2> Get_M(int r,int c) const
    {
        MATRIX<T,2> m;
        M.Get_Submatrix(r,c,m);
        return m;
    }
    TV Get_C(int r,int c) const
    {
        TV u;
        for(int i=0;i<2;i++) u(i)=M(r+i,c);
        return u;
    }
    TV Get_R(int r,int c) const
    {
        TV u;
        for(int i=0;i<2;i++) u(i)=M(r,c+i);
        return u;
    }

    void Add_M(int r,int c,MATRIX<T,2> m) {M.Add_To_Submatrix(r,c,m);}
    void Add_C(int r,int c,TV u) {for(int i=0;i<2;i++) M(r+i,c)+=u(i);}
    void Add_R(int r,int c,TV u) {for(int i=0;i<2;i++) M(r,c+i)+=u(i);}

    MATRIX<T,2> Get_vv(int r,int c) const {return Get_M(2*r,2*c);}
    MATRIX<T,2> Get_ve(int r,int c) const {return Get_M(2*r,2*(c+nv_c));}
    MATRIX<T,2> Get_ev(int r,int c) const {return Get_M(2*(r+nv_r),2*c);}
    MATRIX<T,2> Get_ee(int r,int c) const {return Get_M(2*(r+nv_r),2*(c+nv_c));}
    MATRIX<T,2> Get_uu(int r,int er,int c,int ec) const {return Get_M(2*(r+er*nv_r),2*(c+ec*nv_c));}
    TV Get_vp(int r,int c) const {return Get_C(2*r,2*(nv_c+ne_c)+c);}
    TV Get_pv(int r,int c) const {return Get_R(2*(nv_r+ne_r)+r,2*c);}
    TV Get_ep(int r,int c) const {return Get_C(2*(r+nv_r),2*(nv_c+ne_c)+c);}
    TV Get_pe(int r,int c) const {return Get_R(2*(nv_r+ne_r)+r,2*(c+nv_c));}
    TV Get_up(int r,int er,int c) const {return Get_C(2*(r+er*nv_r),2*(nv_c+ne_c)+c);}
    TV Get_pu(int r,int c,int ec) const {return Get_R(2*(nv_r+ne_r)+r,2*(c+ec*nv_c));}

    void Add_vv(int r,int c,MATRIX<T,2> m) {Add_M(2*r,2*c,m);}
    void Add_ve(int r,int c,MATRIX<T,2> m) {Add_M(2*r,2*(c+nv_c),m);}
    void Add_ev(int r,int c,MATRIX<T,2> m) {Add_M(2*(r+nv_r),2*c,m);}
    void Add_ee(int r,int c,MATRIX<T,2> m) {Add_M(2*(r+nv_r),2*(c+nv_c),m);}
    void Add_uu(int r,int er,int c,int ec,MATRIX<T,2> m) {Add_M(2*(r+er*nv_r),2*(c+ec*nv_c),m);}
    void Add_vp(int r,int c,TV u) {Add_C(2*r,2*(nv_c+ne_c)+c,u);}
    void Add_pv(int r,int c,TV u) {Add_R(2*(nv_r+ne_r)+r,2*c,u);}
    void Add_ep(int r,int c,TV u) {Add_C(2*(r+nv_r),2*(nv_c+ne_c)+c,u);}
    void Add_pe(int r,int c,TV u) {Add_R(2*(nv_r+ne_r)+r,2*(c+nv_c),u);}
    void Add_up(int r,int er,int c,TV u) {Add_C(2*(r+er*nv_r),2*(nv_c+ne_c)+c,u);}
    void Add_pu(int r,int c,int ec,TV u) {Add_R(2*(nv_r+ne_r)+r,2*(c+ec*nv_c),u);}
};

}
#endif