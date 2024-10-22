//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_ENO_RF  
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_PDE/Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Function Conservation_Solver
//#####################################################################
// psi is size (0,m) - U is size 3 by (-3,m+3) with 3 ghost cells - Fx is size 3 by (0,m)
template<class TV,int d> void CONSERVATION_ENO_RF<TV,d>::
Conservation_Solver(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,d>& eigensystem,EIGENSYSTEM<T,d>& eigensystem_explicit,
        const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux)
{
    switch(order){
        case 1:Conservation_Solver_Helper<1>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries);break;
        case 2:Conservation_Solver_Helper<2>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries);break;
        case 3:Conservation_Solver_Helper<3>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries);break;
        default:PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Conservation_Solver_Helper
//#####################################################################
template<class TV,int d> template<int eno_order> void CONSERVATION_ENO_RF<TV,d>::
Conservation_Solver_Helper(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,d>& eigensystem,EIGENSYSTEM<T,d>& eigensystem_explicit,
        const VECTOR<bool,2>& outflow_boundaries)
{
    int k,i,j;

    // divided differences
    ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > DU(0,d,-3,m+3),DF(0,d,-3,m+3);
    ARRAY<TV_DIMENSION,VECTOR<int,1> > F(-3,m+3);eigensystem.Flux(m,U,F); 
    for(i=-3;i<m+3;i++) for(k=0;k<d;k++){DU(k,i)(0)=U(i)(k);DF(k,i)(0)=F(i)(k);}
    for(j=2;j<eno_order;j++) for(k=0;k<d;k++) for(i=-3;i<m+4-j;i++){DU(k,i)(j)=(DU(k,i+1)(j-1)-DU(k,i)(j-1))/((j+1)*dx);DF(k,i)(j)=(DF(k,i+1)(j-1)-DF(k,i)(j-1))/((j+1)*dx);}

    // calculate the fluxes 
    ARRAY<bool,VECTOR<int,1> > psi_ghost(0,m+1);ARRAY<bool,VECTOR<int,1> >::Put(psi,psi_ghost); // ghost points for the if statement below  
    ARRAY<TV_DIMENSION,VECTOR<int,1> > flux(0,m); // fluxes to the right of each point
    VECTOR<T,d> lambda,lambda_left,lambda_right;
    MATRIX<T,d,d> L,R;
    ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > LDU(0,d,-3,m+3),LDF(0,d,-3,m+3);
    ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > *Dstate_ptr,*Dflux_ptr;
    for(i=0;i<m;i++) if(psi_ghost(i) || psi_ghost(i+1)){ // compute flux
        // eigensystem
        eigensystem.Eigenvalues(U,i,lambda,lambda_left,lambda_right);
        if(!eigensystem.All_Eigenvalues_Same()){
            eigensystem.Eigenvectors(U,i,L,R);
            // transfer the divided differences into the characteristic fields
            for(j=0;j<eno_order;j++) for(int ii=i-j;ii<i+1;ii++) if(ii >= -3 && ii < m+3-j) for(k=0;k<d;k++){
                LDU(k,ii)(j)=LDF(k,ii)(j)=0;
                for(int kk=0;kk<d;kk++){LDU(k,ii)(j)+=L(k,kk)*DU(kk,ii)(j);LDF(k,ii)(j)+=L(k,kk)*DF(kk,ii)(j);}}
            Dstate_ptr=&LDU;Dflux_ptr=&LDF;}
        else{Dstate_ptr=&DU;Dflux_ptr=&DF;}
        ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > &Dstate=*Dstate_ptr,&Dflux=*Dflux_ptr;
        // find a flux in each characteristic field
        T flux_total;
        if(eno_order == 1) for(k=0;k<d;k++){   
            if(lambda_left(k)*lambda_right(k) > 0)
                if(lambda(k) > 0) flux_total=Dflux(k,i)(0);
                else flux_total=Dflux(k,i+1)(0);
            else{
                T alpha=Alpha(lambda_left,lambda_right,k,d);
                T flux_left=Dflux(k,i)(0)+alpha*Dstate(k,i)(0);
                T flux_right=Dflux(k,i+1)(0)-alpha*Dstate(k,i+1)(0);
                flux_total=(T).5*(flux_left+flux_right);}
            if(!eigensystem.All_Eigenvalues_Same()) for(int kk=0;kk<d;kk++) flux(i)(kk)+=flux_total*R(k,kk);
            else flux(i)(k)=flux_total;}
        else if(eno_order == 2) for(k=0;k<d;k++){ 
            if(lambda_left(k)*lambda_right(k) > 0)
                if(lambda(k) > 0) flux_total=ADVECTION_SEPARABLE_UNIFORM<TV,T>::ENO(dx,Dflux(k,i)(0),Dflux(k,i-1)(1),Dflux(k,i)(1));
                else flux_total=ADVECTION_SEPARABLE_UNIFORM<TV,T>::ENO(dx,Dflux(k,i+1)(0),-Dflux(k,i+1)(1),-Dflux(k,i)(1));
            else{
                T alpha=Alpha(lambda_left,lambda_right,k,d);
                T flux_left=ADVECTION_SEPARABLE_UNIFORM<TV,T>::ENO(dx,Dflux(k,i)(0)+alpha*Dstate(k,i)(0),Dflux(k,i-1)(1)+alpha*Dstate(k,i-1)(1),Dflux(k,i)(1)+alpha*Dstate(k,i)(1));
                T flux_right=ADVECTION_SEPARABLE_UNIFORM<TV,T>::ENO(dx,Dflux(k,i+1)(0)-alpha*Dstate(k,i+1)(0),-(Dflux(k,i+1)(1)-alpha*Dstate(k,i+1)(1)),-(Dflux(k,i)(1)-alpha*Dstate(k,i)(1)));
                flux_total=(T).5*(flux_left+flux_right);}
            if(!eigensystem.All_Eigenvalues_Same()) for(int kk=0;kk<d;kk++) flux(i)(kk)+=flux_total*R(k,kk);
            else flux(i)(k)=flux_total;}
        else if(eno_order == 3) for(k=0;k<d;k++){
            if(lambda_left(k)*lambda_right(k) > 0)
                if(lambda(k) > 0) flux_total=ADVECTION_SEPARABLE_UNIFORM<TV,T>::ENO(dx,Dflux(k,i)(0),Dflux(k,i-1)(1),Dflux(k,i)(1),Dflux(k,i-3)(2),Dflux(k,i-1)(2),Dflux(k,i)(2));
                else flux_total=ADVECTION_SEPARABLE_UNIFORM<TV,T>::ENO(dx,Dflux(k,i+1)(0),-Dflux(k,i+1)(1),-Dflux(k,i)(1),Dflux(k,i+1)(2),Dflux(k,i)(2),Dflux(k,i-1)(2));
            else{
                T alpha=Alpha(lambda_left,lambda_right,k,d);
                T flux_left=ADVECTION_SEPARABLE_UNIFORM<TV,T>::ENO(dx,Dflux(k,i)(0)+alpha*Dstate(k,i)(0),Dflux(k,i-1)(1)+alpha*Dstate(k,i-1)(1),Dflux(k,i)(1)+alpha*Dstate(k,i)(1),
                    Dflux(k,i-3)(2)+alpha*Dstate(k,i-3)(2),Dflux(k,i-1)(2)+alpha*Dstate(k,i-1)(2),Dflux(k,i)(2)+alpha*Dstate(k,i)(2));
                T flux_right=ADVECTION_SEPARABLE_UNIFORM<TV,T>::ENO(dx,Dflux(k,i+1)(0)-alpha*Dstate(k,i+1)(0),-(Dflux(k,i+1)(1)-alpha*Dstate(k,i+1)(1)),-(Dflux(k,i)(1)-alpha*Dstate(k,i)(1)),
                    Dflux(k,i+1)(2)-alpha*Dstate(k,i+1)(2),Dflux(k,i)(2)-alpha*Dstate(k,i)(2),Dflux(k,i-1)(2)-alpha*Dstate(k,i-1)(2));
                flux_total=(T).5*(flux_left+flux_right);}
            if(!eigensystem.All_Eigenvalues_Same()) for(int kk=0;kk<d;kk++) flux(i)(kk)+=flux_total*R(k,kk);
            else flux(i)(k)=flux_total;}}

    // difference the fluxes
    T one_over_dx=1/dx;
    for(i=0;i<m;i++) if(psi_ghost(i)) Fx(i)=(flux(i)-flux(i-1))*one_over_dx;
    if(save_fluxes) for(i=0;i<m;i++) if(psi_ghost(i) || psi_ghost(i+1)) flux_temp(i)=flux(i);
}
//#####################################################################
#define INSTANTIATION_HELPER(TV) \
    template class CONSERVATION_ENO_RF<GRID<TV>,1>; \
    template class CONSERVATION_ENO_RF<GRID<TV>,2>; \
    template class CONSERVATION_ENO_RF<GRID<TV>,3>; \
    template class CONSERVATION_ENO_RF<GRID<TV>,4>; \
    template class CONSERVATION_ENO_RF<GRID<TV>,5>; \
    template class CONSERVATION_ENO_RF<GRID<TV>,6>;
#define P(...) __VA_ARGS__
#if 0 // broken
INSTANTIATION_HELPER(P(VECTOR<float,1>))
INSTANTIATION_HELPER(P(VECTOR<float,2>))
INSTANTIATION_HELPER(P(VECTOR<float,3>))
INSTANTIATION_HELPER(P(VECTOR<double,1>))
INSTANTIATION_HELPER(P(VECTOR<double,2>))
INSTANTIATION_HELPER(P(VECTOR<double,3>))
#endif
