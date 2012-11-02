//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTON_JACOBI_2D  
//##################################################################### 
//
// Input a GRID_2D class.
// Input phi as (1,m) by (1,n).
//
//#####################################################################
#ifndef __HAMILTON_JACOBI_2D__
#define __HAMILTON_JACOBI_2D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTON_JACOBI.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_UNIFORM.h>
namespace PhysBAM{

template<class T> class HAMILTONIAN_2D;

template<class T_input>
class HAMILTON_JACOBI_2D:public HAMILTON_JACOBI,public LEVELSET<VECTOR<T_input,2> >,public LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef LEVELSET<TV> BASE;
    using BASE::grid;using BASE::phi;using BASE::boundary;using BASE::max_time_step;
    using BASE::curvature;using BASE::curvature_motion;using BASE::sigma;using BASE::Compute_Curvature;
    using LEVELSET_ADVECTION_UNIFORM<GRID<TV> >::HJ_WENO;using LEVELSET_ADVECTION_UNIFORM<GRID<TV> >::HJ_ENO;

    HAMILTONIAN_2D<T>& hamiltonian;

    HAMILTON_JACOBI_2D(HAMILTONIAN_2D<T>& hamiltonian_input,GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input) 
        :LEVELSET<TV>(grid_input,phi_input),LEVELSET_ADVECTION_UNIFORM<GRID<TV> >((LEVELSET<TV>*)this),hamiltonian(hamiltonian_input)
    {}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    void Calculate_Derivatives(ARRAY<T,TV_INT>& phi_ghost,ARRAY<T,TV_INT>& phix_minus,ARRAY<T,TV_INT>& phix_plus,ARRAY<T,TV_INT>& phiy_minus,ARRAY<T,TV_INT>& phiy_plus);
    T CFL(const T time=0);
//#####################################################################
};   
}
#endif
