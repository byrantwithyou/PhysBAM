//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_PDE/Advection/ADVECTION_CONSERVATIVE_WENO.h>
#include <Grid_PDE/Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{
//#####################################################################
// Function Advection_Solver
//#####################################################################
// finds (uZ)_x with Local Lax Friedrichs
template<class T,class T2> void ADVECTION_CONSERVATIVE_WENO<T,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    ARRAY<T2,VECTOR<int,1> > DZ0(-2,m+3); // 1st divided difference
    for(int i=-3;i<m+3;i++) DZ0(i)=Z(i);
    ARRAY<T2,VECTOR<int,1> > DUZ0(-2,m+3); // 1st divided difference
    for(int i=-3;i<m+3;i++) DUZ0(i)=u(i)*Z(i);

    if(compute_epsilon){
        epsilon=1e-6*sqr(maxabs(DUZ0)); // only DUZ used to find epsilon 
        if(epsilon == 0) epsilon=1e-6;} // epsilon=0 implies all v_i=0 and all u_Zx=0

    ARRAY<T2,VECTOR<int,1> > flux(0,m); // flux is to the right of each point 
    for(int i=0;i<m;i++){
        T2 alpha=maxabs(u(i),u(i+1));           
        T2 fluxleft=WENO(DUZ0(i-2)+alpha*DZ0(i-2),DUZ0(i-1)+alpha*DZ0(i-1),DUZ0(i)+alpha*DZ0(i),DUZ0(i+1)+alpha*DZ0(i+1),DUZ0(i+2)+alpha*DZ0(i+2),epsilon);
        T2 fluxright=WENO(DUZ0(i+3)-alpha*DZ0(i+3),DUZ0(i+2)-alpha*DZ0(i+2),DUZ0(i+1)-alpha*DZ0(i+1),DUZ0(i)-alpha*DZ0(i),DUZ0(i-1)-alpha*DZ0(i-1),epsilon);
        flux(i)=.5*(fluxleft+fluxright);}

    T one_over_dx=1/dx;
    for(int i=0;i<m;i++) u_Zx(i)=(flux(i)-flux(i-1))*one_over_dx;
}
}

