//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,3> TV;
typedef MATRIX<T,3> TM;
typedef DIAGONAL_MATRIX<T,3> TDM;
typedef SYMMETRIC_MATRIX<T,3> TSM;

void Compute(const VECTOR<TV,4>& x,T& E,VECTOR<TV,4>& DE,VECTOR<VECTOR<TM,4>,4>& DDE)
{
    // ENERGY

    TV a=x(0)-x(1);
    TV b=x(3)-x(1);
    TV d=x(2)-x(1);

    TV mu=TV::Cross_Product(d,b);
    TV nu=TV::Cross_Product(a,d);

    T one_over_d_norm=(T)1/d.Magnitude();
    T one_over_mu_norm=(T)1/mu.Magnitude();
    T one_over_nu_norm=(T)1/nu.Magnitude();

    TV d_hat=d*one_over_d_norm;
    TV mu_hat=mu*one_over_mu_norm;
    TV nu_hat=nu*one_over_nu_norm;

    T s=TV::Triple_Product(nu_hat,mu_hat,d_hat);
    T c=TV::Dot_Product(nu_hat,mu_hat);

    T theta=atan2(s,c);
    E=(T).5*sqr(theta);

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
}


int main()
{
    const T e=1e-6;
    RANDOM_NUMBERS<T> rand;
    VECTOR<TV,4> x,dx,x_old,x_new;

    for(int i=0;i<4;i++){
        rand.Fill_Uniform(x(i),(T)(-1),(T)1);
        rand.Fill_Uniform(dx(i),-e,e);
        x_new(i)=x(i)+dx(i);
        x_old(i)=x(i)-dx(i);}

    T E,E_old,E_new;
    VECTOR<TV,4> DE,DE_old,DE_new;
    VECTOR<VECTOR<TM,4>,4> DDE,DDE_old,DDE_new;
    VECTOR<TV,4> delta_DE;
    Compute(x,E,DE,DDE);
    Compute(x_old,E_old,DE_old,DDE_old);
    Compute(x_new,E_new,DE_new,DDE_new);

    T delta_E=0;
    for(int i=0;i<4;i++) delta_E+=TV::Dot_Product(DE(i),dx(i));
    LOG::cout<<(E_new-E_old-delta_E*2)/e<<std::endl;

    for(int mi=0;mi<4;mi++)
    for(int ni=0;ni<4;ni++)
        delta_DE(mi)+=DDE(mi)(ni)*dx(ni);
        
    T norm=0;
    for(int i=0;i<4;i++)
        norm+=(DE_new(i)-DE_old(i)-delta_DE(i)*2).Magnitude_Squared();
    LOG::cout<<sqrt(norm)/e<<std::endl;

    return 0;
}
