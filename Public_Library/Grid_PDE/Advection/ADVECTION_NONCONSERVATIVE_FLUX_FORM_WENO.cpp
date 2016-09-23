//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_PDE/Advection/ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO.h>
#include <Grid_PDE/Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{
//#####################################################################
// Function Advection_Solver
//#####################################################################
template<class T,class T2> void ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO<T,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    T one_over_dx=1/dx;
    ARRAY<T2,VECTOR<int,1> > D0(-2,m+3); // 1st divided difference
    for(int i=-3;i<m+3;i++) D0(i)=Z(i);

    if(compute_epsilon){
        epsilon=1e-6*sqr(maxabs(D0)); 
        if(epsilon == 0) epsilon=1e-6;} // epsilon=0 implies all v_i=0 and all u_Zx=0

    for(int i=0;i<m;i++){
        if(u(i) > 0){
            T2 flux_left=WENO(D0(i-3),D0(i-2),D0(i-1),D0(i),D0(i+1),epsilon),flux_right=WENO(D0(i-2),D0(i-1),D0(i),D0(i+1),D0(i+2),epsilon);
            u_Zx(i)=u(i)*(flux_right-flux_left)*one_over_dx;}
        else{ 
            T2 flux_left=WENO(D0(i+2),D0(i+1),D0(i),D0(i-1),D0(i-2),epsilon),flux_right=WENO(D0(i+3),D0(i+2),D0(i+1),D0(i),D0(i-1),epsilon);
            u_Zx(i)=u(i)*(flux_right-flux_left)*one_over_dx;}}
}
}

