//#####################################################################
// Copyright 2002-2005, Doug Enright, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_1D_HELPER.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> LEVELSET<VECTOR<T,1> >::
LEVELSET(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,const int number_of_ghost_cells_input)
    :LEVELSET_BASE<TV>(grid_input,phi_input,number_of_ghost_cells_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> LEVELSET<VECTOR<T,1> >::
~LEVELSET()
{}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,0> LEVELSET<VECTOR<T,1> >::
Principal_Curvatures(const VECTOR<T,1>& X) const
{
    return VECTOR<T,0>(); // not much curvature in 1D
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
// note that abs(phix)=1 if it's a distance function
template<class T> void LEVELSET<VECTOR<T,1> >::
Compute_Normals(const T time)
{
    T one_over_two_dx=1/(2*grid.dX.x);
    int ghost_cells=3;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
        
    if(!normals) normals=new ARRAY<VECTOR<T,1> ,TV_INT>(grid.Domain_Indices(ghost_cells-1));
    for(int i=normals->domain.min_corner.x;i<normals->domain.max_corner.x;i++)(*normals)(i)=VECTOR<T,1>((phi_ghost(i+1)-phi_ghost(i-1))*one_over_two_dx).Normalized();
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T> void LEVELSET<VECTOR<T,1> >::
Compute_Curvature(const T time)
{      
    if(!curvature) curvature=new ARRAY<T,TV_INT>(grid.counts.x,2);
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T> void LEVELSET<VECTOR<T,1> >::
Fast_Marching_Method(const T time,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells,int process_sign)
{
    Get_Signed_Distance_Using_FMM(phi,time,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells,process_sign);
}
//#####################################################################
// Function Get_Signed_Distance_Using_FMM
//#####################################################################
template<class T> void LEVELSET<VECTOR<T,1> >::
Get_Signed_Distance_Using_FMM(ARRAY<T,TV_INT>& signed_distance,const T time,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,
        const bool add_seed_indices_for_ghost_cells,int process_sign)
{       
    const int ghost_cells=2*number_of_ghost_cells+1;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(*this,ghost_cells);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells,process_sign);
    ARRAY<T,TV_INT>::Get(signed_distance,phi_ghost);
    boundary->Apply_Boundary_Condition(grid,signed_distance,time);
}
//#####################################################################
template class LEVELSET<VECTOR<float,1> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET<VECTOR<double,1> >;
#endif
