//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_SYSTEM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

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
template<class TV,class T_ARRAY> typename TV::SCALAR
Bending_Energy(const VECTOR<TV,3>& x,ARRAY_BASE<TV,T_ARRAY>* DE,VECTOR<VECTOR<MATRIX<typename TV::SCALAR,2>,3>,3>* DDE)
{
    PHYSBAM_FATAL_ERROR();
}

//#####################################################################
// Function Bending_Energy
//#####################################################################
template<class TV,class T_ARRAY> typename TV::SCALAR
Bending_Energy(const VECTOR<TV,4>& x,ARRAY_BASE<TV,T_ARRAY>* DE,VECTOR<VECTOR<MATRIX<typename TV::SCALAR,3>,4>,4>* DDE)
{
    typedef typename TV::SCALAR T;
    typedef DIAGONAL_MATRIX<T,TV::m> TDM;
    typedef SYMMETRIC_MATRIX<T,TV::m> TSM;
    typedef MATRIX<T,TV::m> TM;

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

    if(!DE && !DDE) return E;
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
    for(int mj=0;mj<3;mj++){
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
        if(DE) (*DE)(mi)(mj)+=theta*Dtheta(mi)(mj);}

    if(!DDE) return E;
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

        (*DDE)(mi)(ni)(mj,nj)=
            Dtheta(ni)(nj)*Dtheta(mi)(mj)+
            theta*DDtheta(mi)(ni)(mj,nj);
    }
    
    return E;
}
//#####################################################################
// Function Set_Matrix_Block_And_Rhs
//#####################################################################
template<class TV> typename TV::SCALAR MARCHING_CUBES_SYSTEM<TV>::
Set_Matrix_Block_And_Rhs(const VECTOR<int,TV::m+1> index,const VECTOR<TV,TV::m+1> particles,INDIRECT_ARRAY<ARRAY<TV>,VECTOR<int,TV::m+1>&>* rhs)
{
    blocks.Add_End();
    BLOCK& block=blocks.Last();
    block.index=index;
    T E=Bending_Energy(particles,rhs,&block.matrix);
    return E;
}
//#####################################################################
// Function Set_Matrix_And_Rhs
//#####################################################################
template<class TV> typename TV::SCALAR MARCHING_CUBES_SYSTEM<TV>::
Set_Matrix_And_Rhs(MARCHING_CUBES_VECTOR<TV>& rhs,const ARRAY<VECTOR<int,TV::m+1> >& active_list,const ARRAY<int>& index_map,const ARRAY<int>& reverse_index_map,ARRAY_VIEW<const TV> X)
{
    blocks.Remove_All();
    rhs.x.Fill(TV());
    ARRAY<TV> rhs_full(X.m);
    T E=0;
    for(int i=0;i<active_list.m;i++){
        const VECTOR<int,TV::m+1>& nodes=active_list(i);
        INDIRECT_ARRAY<ARRAY<TV>,VECTOR<int,TV::m+1>&> ia(rhs_full.Subset(nodes));
        E+=Set_Matrix_Block_And_Rhs((VECTOR<int,TV::m+1>)reverse_index_map.Subset(nodes),
            (VECTOR<TV,TV::m+1>)X.Subset(nodes),&ia);}
    rhs.x=rhs_full.Subset(index_map);
    return E;
}
//#####################################################################
// Function Set_Matrix_And_Rhs
//#####################################################################
template<class TV> typename TV::SCALAR MARCHING_CUBES_SYSTEM<TV>::
Set_Rhs(MARCHING_CUBES_VECTOR<TV>& rhs,const ARRAY<VECTOR<int,TV::m+1> >& active_list,const ARRAY<int>& index_map,ARRAY_VIEW<const TV> X)
{
    rhs.x.Fill(TV());
    ARRAY<TV> rhs_full(X.m);
    T E=0;
    for(int i=0;i<active_list.m;i++){
        const VECTOR<int,TV::m+1>& nodes=active_list(i);
        INDIRECT_ARRAY<ARRAY<TV>,VECTOR<int,TV::m+1>&> ia(rhs_full.Subset(nodes));
        VECTOR<TV,TV::m+1> X(X.Subset(nodes));
        E+=Bending_Energy(X,&ia,(VECTOR<VECTOR<TM,TV::m+1>,TV::m+1>*)0);}
    rhs.x=rhs_full.Subset(index_map);
    return E;
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MARCHING_CUBES_SYSTEM<TV>::
Compute_Energy(const ARRAY<VECTOR<int,TV::m+1> >& active_list,ARRAY_VIEW<const TV> X)
{
    T E=0;
    for(int i=0;i<active_list.m;i++){
        const VECTOR<int,TV::m+1>& nodes=active_list(i);
        VECTOR<TV,TV::m+1> VX(X.Subset(nodes));
        E+=Bending_Energy(VX,(INDIRECT_ARRAY<ARRAY<TV>,VECTOR<int,TV::m+1>&>*)0,(VECTOR<VECTOR<TM,TV::m+1>,TV::m+1>*)0);}
    return E;
}
//#####################################################################
// Function Compute_Active_List
//#####################################################################
template<class TV> void MARCHING_CUBES_SYSTEM<TV>::
Compute_Active_List(ARRAY<VECTOR<int,TV::m+1> >& active_list,const HASHTABLE<VECTOR<int,2>,T_SURFACE*>& interface,const ARRAY<int>& reverse_index_map)
{
    for(typename HASHTABLE<VECTOR<int,2>,T_SURFACE*>::CONST_ITERATOR it(interface);it.Valid();it.Next()){
        T_SURFACE& surf=*it.Data();
        surf.Update_Number_Nodes();
        surf.mesh.Initialize_Adjacent_Elements();
        const ARRAY<ARRAY<int> >& adjacent_elements=*surf.mesh.adjacent_elements;
        for(int i=0;i<adjacent_elements.m;i++){
            for(int j=0;j<adjacent_elements(i).m;j++){
                int k=adjacent_elements(i)(j);
                if(i<k){
                    VECTOR<int,TV::m+1> nodes;
                    TV_INT ei=surf.mesh.elements(i);
                    TV_INT ek=surf.mesh.elements(k);
                    int u=-1;
                    for(int m=0;m<TV::m;m++)
                        if(!ek.Contains(ei(m))){
                            u=m;
                            break;}
                    for(int m=0;m<TV::m;m++){
                        nodes(m)=ei(u++);
                        if(u==TV::m) u=0;}
                    nodes(TV::m)=ek.Sum()-ei.Sum()+nodes(0);
                    if(reverse_index_map.Subset(nodes).Count_Matches(-1)!=TV::m+1)
                        active_list.Append(nodes);}}}}
}
//#####################################################################
// Function Set_Matrix_Block_And_Rhs
//#####################################################################
template<class TV> void MARCHING_CUBES_SYSTEM<TV>::
Test_System(const ARRAY<VECTOR<int,TV::m+1> >& active_list,const ARRAY<int>& index_map,const ARRAY<int>& reverse_index_map)
{
    T e=(T)1e-6;
    MARCHING_CUBES_SYSTEM<TV> system0,system1;
    MARCHING_CUBES_VECTOR<TV> rhs0,rhs1,vec,s,t;
    RANDOM_NUMBERS<T> random;
    ARRAY<TV> X0(reverse_index_map.m),dX(reverse_index_map.m);
    vec.x.Resize(index_map.m);
    s.x.Resize(index_map.m);
    t.x.Resize(index_map.m);
    random.Fill_Uniform(X0,-(T)1,(T)1);
    random.Fill_Uniform(vec.x,-(T)e,(T)e);
    dX.Subset(index_map)=vec.x;
    ARRAY<TV> X1(X0+dX);
    T E0=system0.Set_Matrix_And_Rhs(rhs0,active_list,index_map,reverse_index_map,X0);
    T E1=system1.Set_Matrix_And_Rhs(rhs1,active_list,index_map,reverse_index_map,X1);
    LOG::cout<<"E "<<((E1-E0)-(rhs1.x.Dot(vec.x)+rhs0.x.Dot(vec.x))/2)/e/maxabs(E0,E1,(T)1e-30)<<std::endl;

    system0.Multiply(vec,s);
    system1.Multiply(vec,t);
    T dif=(rhs1.x-rhs0.x-(T).5*s.x-(T).5*t.x).Magnitude();
    LOG::cout<<"dE "<<dif/e/max(rhs1.x.Magnitude(),rhs0.x.Magnitude(),(T)1e-30)<<std::endl;
}
template class MARCHING_CUBES_SYSTEM<VECTOR<float,2> >;
template class MARCHING_CUBES_SYSTEM<VECTOR<float,3> >;
template class MARCHING_CUBES_SYSTEM<VECTOR<double,2> >;
template class MARCHING_CUBES_SYSTEM<VECTOR<double,3> >;
