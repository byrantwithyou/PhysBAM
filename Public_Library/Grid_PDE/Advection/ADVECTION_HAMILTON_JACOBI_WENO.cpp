//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_HAMILTON_JACOBI_WENO  
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_PDE/Advection/ADVECTION_HAMILTON_JACOBI_WENO.h>
using namespace PhysBAM;
//##################################################################### 
// Function Advection_Solver
//#####################################################################
template<class TV,class T2> void ADVECTION_HAMILTON_JACOBI_WENO<TV,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{      
    int i;T one_over_dx=1/dx;
    RANGE<VECTOR<int,1> > R;
    R.min_corner.x=-2;
    R.max_corner.x=m+2;
    ARRAY<T2,VECTOR<int,1> > D0(R);
    for(i=-3;i<m+2;i++) D0(i)=(Z(i+1)-Z(i))*one_over_dx;

    if(compute_epsilon){
        epsilon=(T)1e-6*sqr(D0.Max_Abs());
        if(epsilon == 0) epsilon=(T)1e-6;} // epsilon=0 implies all v_i=0 and all u_Zx=0  

    for(i=0;i<m;i++){
        if(u(i) > 0) u_Zx(i)=u(i)*WENO(D0(i-3),D0(i-2),D0(i-1),D0(i),D0(i+1),epsilon);
        else u_Zx(i)=u(i)*WENO(D0(i+2),D0(i+1),D0(i),D0(i-1),D0(i-2),epsilon);}
}
//#####################################################################
namespace PhysBAM{
template class ADVECTION_HAMILTON_JACOBI_WENO<VECTOR<float,1>,float>;
template class ADVECTION_HAMILTON_JACOBI_WENO<VECTOR<float,2>,float>;
template class ADVECTION_HAMILTON_JACOBI_WENO<VECTOR<float,3>,float>;
template class ADVECTION_HAMILTON_JACOBI_WENO<VECTOR<double,1>,double>;
template class ADVECTION_HAMILTON_JACOBI_WENO<VECTOR<double,2>,double>;
template class ADVECTION_HAMILTON_JACOBI_WENO<VECTOR<double,3>,double>;
}
