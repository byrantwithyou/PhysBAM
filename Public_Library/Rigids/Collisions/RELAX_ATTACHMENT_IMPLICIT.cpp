//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Rigids/Collisions/RELAX_ATTACHMENT_IMPLICIT.h>
using namespace PhysBAM;
//#####################################################################
// Function Relax
//#####################################################################
template<class TV> void RELAX_ATTACHMENT_IMPLICIT<TV>::
Relax(const TV& Z,const TV& X,const TV& W,T mu)
{
    TV ZW=Z-W,XW=X-W;
    T zw=ZW.Normalize();
    T xw=XW.Normalize();
    T q=mu*zw;
    if(xw<=q){K=X;dKdX+=1;dynamic=false;/*LOG::puts("stick");*/return;}
    K=W+q*XW;
    dKdX=q/xw*((T)1-Outer_Product(XW,XW));
    dKdZ=Outer_Product(XW,mu*ZW);
    dKdW=(T)1-dKdX-dKdZ;
    dynamic=true;
    //LOG::puts("dynamic");
}
//#####################################################################
// Function Relax_Search
//#####################################################################
template<class TV> void RELAX_ATTACHMENT_IMPLICIT<TV>::
Relax_Search(const TV& Z,const TV& X,const TV& W,const IMPLICIT_OBJECT<TV>* io,T mu)
{
    T mu_bar=1/(sqr(mu)+1);
    auto criterion=[X,W,Z,io,mu_bar](T s)
        {
            TV K=W+s*(X-W);
            T phi=io->Extended_Phi(K);
            TV N=io->Extended_Normal(K);
            TV Y=K-phi*N;
            TV u=Y-Z;
            return sqr(u.Dot(N))-mu_bar*u.Magnitude_Squared();
        };

    if(criterion(1)>=0){K=X;dKdX+=1;dynamic=false;/*puts("stick");*/return;}
    if(criterion(0)<=0){K=W;dKdW+=1;dynamic=true;/*puts("frictionless");*/return;}
    T hi=1,lo=0;
    while(hi-lo>std::numeric_limits<T>::epsilon()*2){
        T md=lo+(hi-lo)/2;
        if(criterion(md)>=0) lo=md;
        else hi=md;}

    T s=lo;
    TV U=X-W;
    TV K=W+s*U;
    T phi=io->Extended_Phi(K);
    TV N=io->Extended_Normal(K);
    TV Y=K-phi*N;
    TV u=Y-Z;
    T a=u.Dot(N);
    TV r=a*N-mu_bar*u;
    TV v=phi*r-a*u;
    T j=1/(r.Dot(U));
    
    K=K;
    dKdphi=a*j*(1-mu_bar)*U;
    dKdN=Outer_Product(U,j*v);
    dKdZ=Outer_Product(U,j*r);
    dKdW=(1-s)-(1-s)*dKdZ;
    dKdX=s-s*dKdZ;
    
    //puts("dynamic");
    dynamic=true;
}
template class RELAX_ATTACHMENT_IMPLICIT<VECTOR<float,2> >;
template class RELAX_ATTACHMENT_IMPLICIT<VECTOR<float,3> >;
template class RELAX_ATTACHMENT_IMPLICIT<VECTOR<double,2> >;
template class RELAX_ATTACHMENT_IMPLICIT<VECTOR<double,3> >;
