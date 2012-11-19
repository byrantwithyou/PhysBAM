//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>

using namespace PhysBAM;

//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MARCHING_CUBES_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& bv_input,KRYLOV_VECTOR_BASE<T>& bv_result) const PHYSBAM_OVERRIDE
{
    const ARRAY<TV>& x_input=debug_cast<const VECTOR_T&>(bv_input).x;
    ARRAY<TV>& x_result=debug_cast<VECTOR_T&>(bv_result).x;
    x_result.Fill(TV());
    for(int b=0;b<blocks.m;b++){
        const BLOCK& block=blocks(b);
        for(int j=0;j<TV::m+1;j++)
        for(int k=0;k<TV::m+1;k++){
            int a=block.index(j),b=block.index(k);
            if(a>=0 && b>=0) x_result(a)+=block.matrix(j)(k)*x_input(b);}}
}

//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MARCHING_CUBES_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) const PHYSBAM_OVERRIDE
{
    const VECTOR_T& v1=debug_cast<const MARCHING_CUBES_VECTOR<TV>&>(bv1);
    const VECTOR_T& v2=debug_cast<const MARCHING_CUBES_VECTOR<TV>&>(bv2);
    return v1.x.Dot(v2.x);
}

//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MARCHING_CUBES_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE
{
    return sqrt(Inner_Product(bv,bv));
}

//#####################################################################
// Function Bending_Energy
//#####################################################################
template<class TV,class TM> typename TV::SCALAR
Bending_Energy(const VECTOR<TV,3>& x,VECTOR<TV,3>& DE,VECTOR<VECTOR<TM,3>,3>& DDE)
{
    PHYSBAM_FATAL_ERROR();
}

//#####################################################################
// Function Bending_Energy
//#####################################################################
template<class TV,class TM> typename TV::SCALAR
Bending_Energy(const VECTOR<TV,4>& x,VECTOR<TV,4>& DE,VECTOR<VECTOR<TM,4>,4>& DDE)
{
    typedef typename TV::SCALAR T;
    typedef DIAGONAL_MATRIX<T,TV::m> TDM;
    typedef SYMMETRIC_MATRIX<T,TV::m> TSM;

    // ENERGY

    const TV a=x(0)-x(1);
    const TV b=x(3)-x(1);
    const TV d=x(2)-x(1);

    const TV mu=TV::Cross_Product(d,b);
    const TV nu=TV::Cross_Product(a,d);

    const T one_over_d_norm=(T)1/d.Magnitude();
    const T one_over_mu_norm=(T)1/mu.Magnitude();
    const T one_over_nu_norm=(T)1/nu.Magnitude();

    const TV d_hat=d*one_over_d_norm;
    const TV mu_hat=mu*one_over_mu_norm;
    const TV nu_hat=nu*one_over_nu_norm;

    const T s=TV::Triple_Product(nu_hat,mu_hat,d_hat);
    const T c=TV::Dot_Product(nu_hat,mu_hat);

    const T theta=atan2(s,c);
    const T E=(T).5*sqr(theta);

    // GRADIENT

    VECTOR<VECTOR<TV,3>,4> Da,Db,Dd;
    VECTOR<VECTOR<TV,3>,4> Dmu,Dnu;
    VECTOR<VECTOR<TV,3>,4> Dd_hat,Dmu_hat,Dnu_hat;
    VECTOR<TV,4> Ds,Dc,Dtheta;

    for(int i=0;i<3;i++){
        Da(0)(i)(i)=1;Da(1)(i)(i)=-1;
        Db(3)(i)(i)=1;Db(1)(i)(i)=-1;
        Dd(2)(i)(i)=1;Dd(1)(i)(i)=-1;}

    for(int mi=0;mi<4;mi++)
    for(int mj=0;mj<3;mj++)
    {
        Dmu(mi)(mj)=TV::Cross_Product(Dd(mi)(mj),b)+TV::Cross_Product(d,Db(mi)(mj));
        Dnu(mi)(mj)=TV::Cross_Product(Da(mi)(mj),d)+TV::Cross_Product(a,Dd(mi)(mj));
        
        Dd_hat(mi)(mj)=(TDM::Identity_Matrix()-TSM::Outer_Product(d_hat))*Dd(mi)(mj)*one_over_d_norm;
        Dmu_hat(mi)(mj)=(TDM::Identity_Matrix()-TSM::Outer_Product(mu_hat))*Dmu(mi)(mj)*one_over_mu_norm;
        Dnu_hat(mi)(mj)=(TDM::Identity_Matrix()-TSM::Outer_Product(nu_hat))*Dnu(mi)(mj)*one_over_nu_norm;

        Ds(mi)(mj)=
            TV::Triple_Product(Dnu_hat(mi)(mj),mu_hat,d_hat)+
            TV::Triple_Product(nu_hat,Dmu_hat(mi)(mj),d_hat)+
            TV::Triple_Product(nu_hat,mu_hat,Dd_hat(mi)(mj));

        Dc(mi)(mj)=
            TV::Dot_Product(Dnu_hat(mi)(mj),mu_hat)+
            TV::Dot_Product(nu_hat,Dmu_hat(mi)(mj));

        Dtheta(mi)(mj)=c*Ds(mi)(mj)-s*Dc(mi)(mj);
        DE(mi)(mj)=theta*Dtheta(mi)(mj);
    }

    // HESSIAN
    
    VECTOR<VECTOR<VECTOR<VECTOR<TV,3>,3>,4>,4> DDmu,DDnu;
    VECTOR<VECTOR<VECTOR<VECTOR<TV,3>,3>,4>,4> DDd_hat,DDmu_hat,DDnu_hat;
    VECTOR<VECTOR<TM,4>,4> DDs,DDc,DDtheta;

    for(int mi=0;mi<4;mi++)for(int ni=0;ni<4;ni++)
    for(int mj=0;mj<3;mj++)for(int nj=0;nj<3;nj++)
    {
        DDmu(mi)(ni)(mj)(nj)=TV::Cross_Product(Dd(mi)(mj),Db(ni)(nj))+TV::Cross_Product(Dd(ni)(nj),Db(mi)(mj));
        DDnu(mi)(ni)(mj)(nj)=TV::Cross_Product(Da(mi)(mj),Dd(ni)(nj))+TV::Cross_Product(Da(ni)(nj),Dd(mi)(mj));
        
        DDd_hat(mi)(ni)(mj)(nj)=(
            d_hat.Dot(Dd(ni)(nj))*Dd_hat(mi)(mj)+
            d_hat.Dot(Dd(mi)(mj))*Dd_hat(ni)(nj)+
            Dd_hat(ni)(nj).Dot(Dd(mi)(mj))*d_hat)
            *(-one_over_d_norm);
        
        DDmu_hat(mi)(ni)(mj)(nj)=(
            mu_hat.Dot(Dmu(ni)(nj))*Dmu_hat(mi)(mj)+
            mu_hat.Dot(Dmu(mi)(mj))*Dmu_hat(ni)(nj)+
            Dmu_hat(ni)(nj).Dot(Dmu(mi)(mj))*mu_hat-
            (TDM::Identity_Matrix()-TSM::Outer_Product(mu_hat))*DDmu(mi)(ni)(mj)(nj))
            *(-one_over_mu_norm);

        DDnu_hat(mi)(ni)(mj)(nj)=(
            nu_hat.Dot(Dnu(ni)(nj))*Dnu_hat(mi)(mj)+
            nu_hat.Dot(Dnu(mi)(mj))*Dnu_hat(ni)(nj)+
            Dnu_hat(ni)(nj).Dot(Dnu(mi)(mj))*nu_hat-
            (TDM::Identity_Matrix()-TSM::Outer_Product(nu_hat))*DDnu(mi)(ni)(mj)(nj))
            *(-one_over_nu_norm);
        
        DDs(mi)(ni)(mj,nj)=
            TV::Triple_Product(DDnu_hat(mi)(ni)(mj)(nj),mu_hat,d_hat)+
            TV::Triple_Product(Dnu_hat(mi)(mj),Dmu_hat(ni)(nj),d_hat)+
            TV::Triple_Product(Dnu_hat(mi)(mj),mu_hat,Dd_hat(ni)(nj))+
            TV::Triple_Product(Dnu_hat(ni)(nj),Dmu_hat(mi)(mj),d_hat)+
            TV::Triple_Product(nu_hat,DDmu_hat(mi)(ni)(mj)(nj),d_hat)+
            TV::Triple_Product(nu_hat,Dmu_hat(mi)(mj),Dd_hat(ni)(nj))+
            TV::Triple_Product(Dnu_hat(ni)(nj),mu_hat,Dd_hat(mi)(mj))+
            TV::Triple_Product(nu_hat,Dmu_hat(ni)(nj),Dd_hat(mi)(mj))+
            TV::Triple_Product(nu_hat,mu_hat,DDd_hat(mi)(ni)(mj)(nj));

        DDc(mi)(ni)(mj,nj)=
            TV::Dot_Product(Dnu_hat(ni)(nj),Dmu_hat(mi)(mj))+
            TV::Dot_Product(nu_hat,DDmu_hat(mi)(ni)(mj)(nj))+
            TV::Dot_Product(Dnu_hat(mi)(mj),Dmu_hat(ni)(nj))+
            TV::Dot_Product(DDnu_hat(mi)(ni)(mj)(nj),mu_hat);

        DDtheta(mi)(ni)(mj,nj)=
            Dc(ni)(nj)*Ds(mi)(mj)+c*DDs(mi)(ni)(mj,nj)-
            Ds(ni)(nj)*Dc(mi)(mj)-s*DDc(mi)(ni)(mj,nj);

        DDE(mi)(ni)(mj,nj)=
            Dtheta(ni)(nj)*Dtheta(mi)(mj)+
            theta*DDtheta(mi)(ni)(mj,nj);
    }
    
    return E;
}

//#####################################################################
// Function Set_Matrix_Block_And_Rhs
//#####################################################################
template<class TV> typename TV::SCALAR MARCHING_CUBES_SYSTEM<TV>::
Set_Matrix_Block_And_Rhs(const VECTOR<int,TV::m+1> index,const VECTOR<TV,TV::m+1> particles,INDIRECT_ARRAY<ARRAY<TV>,VECTOR<int,TV::m+1>&> rhs)
{
    if(index.Sum()==-index.m) return 0;
    blocks.Add_End();
    BLOCK& block=blocks.Last();
    block.index=index;
    VECTOR<TV,TV::m+1> DE;
    T E=Bending_Energy(particles,DE,block.matrix);
    for(int i=0;i<TV::m+1;i++) rhs(i)+=DE(i);
    return E;
}

template class MARCHING_CUBES_SYSTEM<VECTOR<float,2> >;
template class MARCHING_CUBES_SYSTEM<VECTOR<float,3> >;
template class MARCHING_CUBES_SYSTEM<VECTOR<double,2> >;
template class MARCHING_CUBES_SYSTEM<VECTOR<double,3> >;
