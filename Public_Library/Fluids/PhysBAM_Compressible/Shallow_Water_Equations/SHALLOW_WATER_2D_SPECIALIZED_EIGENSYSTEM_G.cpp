//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G  
//##################################################################### 
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// G(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G<T>::
Flux(const int m,const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& F,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-3;i<m+3;i++){
        F(i)(0)=U(i)(0)*U(i)(1); // h*v
        F(i)(1)=(T).5*sqr(U(i)(1))+gravity*eta_ghost(slice_index.x,i);} // .5*v^2+g*eta
}
//#####################################################################
namespace PhysBAM{
template class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G<float>;
template class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G<double>;
}
