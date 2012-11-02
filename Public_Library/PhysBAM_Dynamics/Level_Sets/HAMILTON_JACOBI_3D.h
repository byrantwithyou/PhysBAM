//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTON_JACOBI_3D  
//##################################################################### 
//
// Input a GRID_3D class.
// Input phi as (1,m) by (1,n) by (1,mn).
//
//#####################################################################
#ifndef __HAMILTON_JACOBI_3D__
#define __HAMILTON_JACOBI_3D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTON_JACOBI.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_UNIFORM.h>
namespace PhysBAM{

template<class T> class HAMILTONIAN_3D;

template<class T_input>
class HAMILTON_JACOBI_3D:public HAMILTON_JACOBI,public LEVELSET<VECTOR<T_input,3> >,LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef LEVELSET<TV> BASE;
    using BASE::grid;using BASE::phi;using BASE::boundary;using BASE::max_time_step;
    using BASE::curvature;using BASE::curvature_motion;using BASE::sigma;using BASE::Compute_Curvature;
    using LEVELSET_ADVECTION_UNIFORM<GRID<TV> >::HJ_WENO;using LEVELSET_ADVECTION_UNIFORM<GRID<TV> >::HJ_ENO;

    HAMILTONIAN_3D<T>& hamiltonian;

    HAMILTON_JACOBI_3D(HAMILTONIAN_3D<T>& hamiltonian_input,GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input) 
        :LEVELSET<TV>(grid_input,phi_input),LEVELSET_ADVECTION_UNIFORM<GRID<TV> >((LEVELSET<TV>*)this),hamiltonian(hamiltonian_input)
    {}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    void Calculate_Derivatives(ARRAY<T,TV_INT>& phi_ghost,ARRAY<T,TV_INT>& phix_minus,ARRAY<T,TV_INT>& phix_plus,ARRAY<T,TV_INT>& phiy_minus,ARRAY<T,TV_INT>& phiy_plus,ARRAY<T,TV_INT>& phiz_minus,ARRAY<T,TV_INT>& phiz_plus);
    T CFL(const T time=0);
//#####################################################################
};   
}
#endif
