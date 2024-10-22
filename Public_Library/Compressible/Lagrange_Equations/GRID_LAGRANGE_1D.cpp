//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Compressible/Lagrange_Equations/GRID_LAGRANGE_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void GRID_LAGRANGE_1D<T>::
Euler_Step(const ARRAY<T,VECTOR<int,1> >& u,const T dt)
{       
    for(int i=0;i<m;i++) x(i)+=dt*u(i); 
}
//#####################################################################
// Function Get_Lengths
//#####################################################################
// size (1,m-1)
template<class T> void GRID_LAGRANGE_1D<T>::
Get_Lengths(ARRAY<T,VECTOR<int,1> >& L)
{       
    for(int i=0;i<m-1;i++) L(i)=abs(x(i+1)-x(i)); 
}
//#####################################################################
// Function Get_Midpoints
//#####################################################################
// size (1,m-1)
template<class T> void GRID_LAGRANGE_1D<T>::
Get_Midpoints(ARRAY<T,VECTOR<int,1> >& M)
{       
    for(int i=0;i<m-1;i++) M(i)=(x(i)+x(i+1))/2; 
}
//#####################################################################
namespace PhysBAM{
template class GRID_LAGRANGE_1D<float>;
template class GRID_LAGRANGE_1D<double>;
}
