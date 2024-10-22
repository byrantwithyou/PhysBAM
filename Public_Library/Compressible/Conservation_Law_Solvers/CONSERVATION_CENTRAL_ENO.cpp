//#####################################################################
// Copyright 2002-2004, Doug Enright, Ronald Fedkiw, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_CENTRAL_ENO  
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/Minmod.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_CENTRAL_ENO.h>
#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Function Conservation_Solver
//#####################################################################
// psi is size (0,m) - U is size 3 by (-3,m+3) with 3 ghost cells - Fx is size 3 by (0,m)
template<class TV,int d> void CONSERVATION_CENTRAL_ENO<TV,d>::
Conservation_Solver(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,d>& eigensystem,EIGENSYSTEM<T,d>& eigensystem_explicit,
        const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux)
{
    switch(order){
        case 1:Conservation_Solver_Helper<1>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries);break;
        case 2:Conservation_Solver_Helper<2>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries);break;
        default: PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Conservation_Solver_Helper
//#####################################################################
template<class TV,int d> template<int eno_order> void CONSERVATION_CENTRAL_ENO<TV,d>::
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
    for(i=0;i<m;i++) if( psi_ghost(i)|| psi_ghost(i+1) ){ // compute flux
        eigensystem.Eigenvalues(U,i,lambda,lambda_left,lambda_right);
        // find a flux in each component
        if(eno_order == 1) for(k=0;k<d;k++){
            T alpha=Alpha(lambda_left,lambda_right,k,d);
            T flux_left=DF(k,i)(0)+alpha*DU(k,i)(0);
            T flux_right=DF(k,i+1)(0)-alpha*DU(k,i+1)(0);
            flux(i)(k)=(T).5*(flux_left+flux_right);}
        else if(eno_order == 2) for(k=0;k<d;k++){
            T alpha=Alpha(lambda_left,lambda_right,k,d);
            T flux_left=Minmod(dx,DF(k,i)(0)+alpha*DU(k,i)(0),DF(k,i-1)(1)+alpha*DU(k,i-1)(1),DF(k,i)(1)+alpha*DU(k,i)(1));
            T flux_right=Minmod(dx,DF(k,i+1)(0)-alpha*DU(k,i+1)(0),-(DF(k,i+1)(1)-alpha*DU(k,i+1)(1)),-(DF(k,i)(1)-alpha*DU(k,i)(1)));
            flux(i)(k)=(T).5*(flux_left+flux_right);}}

    // difference the fluxes
    T one_over_dx=1/dx;
    for(i=0;i<m;i++) if(psi_ghost(i)) Fx(i)=(flux(i)-flux(i-1))*one_over_dx;
}
//#####################################################################
#define INSTANTIATION_HELPER(TV) \
    template class CONSERVATION_CENTRAL_ENO<TV,1>; \
    template class CONSERVATION_CENTRAL_ENO<TV,2>;
#define P(...) __VA_ARGS__
#if 0 // broken
INSTANTIATION_HELPER(P(VECTOR<float,1>))
INSTANTIATION_HELPER(P(VECTOR<float,2>))
INSTANTIATION_HELPER(P(VECTOR<float,3>))
INSTANTIATION_HELPER(P(VECTOR<double,1>))
INSTANTIATION_HELPER(P(VECTOR<double,2>))
INSTANTIATION_HELPER(P(VECTOR<double,3>))
#endif
