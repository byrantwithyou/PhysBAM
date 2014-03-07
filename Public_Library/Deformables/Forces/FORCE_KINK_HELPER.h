//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FORCE_KINK_HELPER
//#####################################################################
#ifndef __FORCE_KINK_HELPER__
#define __FORCE_KINK_HELPER__

#include <Tools/Matrices/MATRIX.h>
namespace PhysBAM{

template<class TV,int d>
class FORCE_KINK_HELPER
{
private:
    typedef typename TV::SCALAR T;typedef VECTOR<TV,d> SV;typedef MATRIX<MATRIX<T,TV::m>,d> SM;
public:

    T t,a,b;
    SV x,z,dc,dc_hat,g,da,db,k,dt,s;
    SM dz,ddc,ddt,dk;

    FORCE_KINK_HELPER()
    {}

    ~FORCE_KINK_HELPER()
    {}

    static void Add_Scalar_Identity(SM& s,T sc)
    {for(int i=0;i<d;i++) s(i,i)+=sc;}

    static void Add_Scalar_Outer_Product(SM& s,T a,const SV& u,const SV& v)
    {for(int i=0;i<d;i++){TV au=a*u(i);for(int j=0;j<d;j++) s(i,j)+=MATRIX<T,TV::m>::Outer_Product(au,v(j));}}

    static T Dot(const SV& a,const SV& b)
    {T r=0;for(int i=0;i<d;i++) r+=a(i).Dot(b(i));return r;}

    static SV Transpose_Times(const SM& m,const SV& u)
    {SV r;for(int i=0;i<d;i++) for(int j=0;j<d;j++) r(j)+=m(i,j).Transpose_Times(u(i));return r;}

    static SV Times(const SM& m,const SV& u)
    {SV r;for(int i=0;i<d;i++) for(int j=0;j<d;j++) r(i)+=m(i,j)*u(j);return r;}

    static void Add_Times(SM& o,const SM& m,const SM& n)
    {for(int i=0;i<d;i++) for(int j=0;j<d;j++) for(int k=0;k<d;k++) o(i,j)+=m(i,k)*n(k,j);}

    static SV Scale(T a,const SV& u)
    {SV r;for(int i=0;i<d;i++) r(i)=a*u(i);return r;}

    static void Add_Scale(SM& o,T a,const SM& m)
    {for(int i=0;i<d;i++) for(int j=0;j<d;j++) o(i,j)+=a*m(i,j);}

    SV Contract_ddz_01(const SV& u,const SV& v)
    {
        SV u_ddc=Transpose_Times(ddc,u),v_ddt=Transpose_Times(ddt,v);
        return Scale(Dot(u_ddc,v),dt)+Scale(Dot(dt,v),u_ddc)+Scale(Dot(dc,u),v_ddt);
    }

    template<class F>
    void Compute(const SV& x_,F func,int use_hess)
    {
        x=x_;
        func(x,dc,&ddc,&t); // Computes c_{,i}, c_{,ij}, and optionally \tau.
        if(!use_hess){
            z=x+Scale(t,dc);
            T c_hat=func(z,dc_hat,0,0);
            LOG::printf("ch %g\n",c_hat);
            Add_Scale(dz,t,ddc);}
        else{
            t=0;
            z=x;
            dc_hat=dc;}
        g=Transpose_Times(ddc,dc_hat);
        k=dc_hat;
        if(!use_hess) k+=Scale(t,g);

        LOG::printf("%P %P\n",dc_hat,dc);
        a=Dot(dc_hat,dc);
        b=1/a;
        dt=Scale(-b,k);
        Add_Scalar_Identity(dz,1);
        Add_Scalar_Outer_Product(dz,1,dc,dt);

        if(!use_hess) return;

        s=Transpose_Times(ddc,dc);
        da=Transpose_Times(dz,s)+g;
        db=Scale(-b*b,da);
        Add_Times(dk,ddc,dz);
        Add_Scalar_Outer_Product(dk,1,g,dt);
        Add_Scale(ddt,-b,dk);
        Add_Scalar_Outer_Product(ddt,-1,k,db);
    }
//#####################################################################
};
}
#endif
