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

void Compute(const VECTOR<TV,4>& x,T& E,VECTOR<VECTOR<T,3>,4>& DE)
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
    VECTOR<VECTOR<T,3>,4> Ds,Dc,Dtheta;

    for(int i=0;i<3;i++){
        Da(0)(i)(i)=1;Da(1)(i)(i)=-1;
        Db(3)(i)(i)=1;Db(1)(i)(i)=-1;
        Dd(2)(i)(i)=1;Dd(1)(i)(i)=-1;}

    for(int m=0;m<4;m++)
    for(int n=0;n<3;n++)
    {
        Dmu(m)(n)=TV::Cross_Product(Dd(m)(n),b)+TV::Cross_Product(d,Db(m)(n));
        Dnu(m)(n)=TV::Cross_Product(Da(m)(n),d)+TV::Cross_Product(a,Dd(m)(n));
        
        Dd_hat(m)(n)=(TDM::Identity_Matrix()-TSM::Outer_Product(d_hat))*Dd(m)(n)*one_over_d_norm;
        Dmu_hat(m)(n)=(TDM::Identity_Matrix()-TSM::Outer_Product(mu_hat))*Dmu(m)(n)*one_over_mu_norm;
        Dnu_hat(m)(n)=(TDM::Identity_Matrix()-TSM::Outer_Product(nu_hat))*Dnu(m)(n)*one_over_nu_norm;

        Ds(m)(n)=
            TV::Triple_Product(Dnu_hat(m)(n),mu_hat,d_hat)+
            TV::Triple_Product(nu_hat,Dmu_hat(m)(n),d_hat)+
            TV::Triple_Product(nu_hat,mu_hat,Dd_hat(m)(n));

        Dc(m)(n)=
            TV::Dot_Product(Dnu_hat(m)(n),mu_hat)+
            TV::Dot_Product(nu_hat,Dmu_hat(m)(n));

        Dtheta(m)(n)=c*Ds(m)(n)-s*Dc(m)(n);
        DE(m)(n)=theta*Dtheta(m)(n);
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
    VECTOR<VECTOR<T,3>,4> DE,DE_old,DE_new;
    Compute(x,E,DE);
    Compute(x_old,E_old,DE_old);
    Compute(x_new,E_new,DE_new);

    T delta_E=0;
    for(int i=0;i<4;i++)
        delta_E+=TV::Dot_Product(DE(i),dx(i));

    LOG::cout<<(E_new-E_old-2*delta_E)/e<<std::endl;
    return 0;
}
