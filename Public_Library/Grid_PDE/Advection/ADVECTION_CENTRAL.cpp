//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_PDE/Advection/ADVECTION_CENTRAL.h>
#include <Grid_PDE/Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{
//#####################################################################
// Function Advection_Solver
//#####################################################################
template<class TV,class T2> void ADVECTION_CENTRAL<TV,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    T one_over_two_dx=1/(2*dx);
    for(int i=0;i<m;i++) u_Zx(i)=u(i)*(Z(i+1)-Z(i-1))*one_over_two_dx;
}
}    
