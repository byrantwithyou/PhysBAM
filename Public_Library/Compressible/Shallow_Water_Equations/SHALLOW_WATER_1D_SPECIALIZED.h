//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_1D_SPECIALIZED  
//##################################################################### 
//
// Input U as 2 by (1,m) for h and h*u
//
//#####################################################################
#ifndef __SHALLOW_WATER_1D_SPECIALIZED__
#define __SHALLOW_WATER_1D_SPECIALIZED__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F.h>
namespace PhysBAM{

template<class T_input>
class SHALLOW_WATER_1D_SPECIALIZED:public SHALLOW_WATER<GRID<VECTOR<T_input,1> > >
{
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<T,2> TV_DIMENSION;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    enum {d=2};
public:
    using SHALLOW_WATER<TV>::boundary;using SHALLOW_WATER<TV>::conservation;

    T gravity;
    T min_height;
    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,VECTOR<int,1> >& U; // h and u
    ARRAY<T,VECTOR<int,1> >* ground;
    ARRAY<T,VECTOR<int,1> > eta_ghost;
    ARRAY<VECTOR<T,4> ,VECTOR<int,1> > postprocessed_flux; // for debugging
    SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F<T> eigensystem_F;

    SHALLOW_WATER_1D_SPECIALIZED(GRID<TV>& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_input,const T gravity_input=9.8,const T min_height_input=1e-3)
        :gravity(gravity_input),min_height(min_height_input),grid(grid_input),U(U_input),ground(0),eta_ghost(-2,grid.counts.x+3),
        postprocessed_flux(0,grid.counts.x),eigensystem_F(gravity_input,eta_ghost,min_height_input)
    {}
    
    void Set_Ground(ARRAY<T,VECTOR<int,1> >& ground_input)
    {ground=&ground_input;}

//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};    
}
#endif
