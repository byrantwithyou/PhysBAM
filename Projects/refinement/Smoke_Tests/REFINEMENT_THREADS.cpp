//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include "REFINEMENT_THREADS.h"
#include "SMOKE_TESTS.h"

using namespace PhysBAM;

template<class TV> void REFINEMENT_TASK<TV>::Run()
{        
    GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*smoke_tests->sub_scale,RANGE<TV>::Centered_Box(),true);
    ARRAY<T,FACE_INDEX<TV::dimension> > local_face_velocities(local_mac_grid);
    PROJECTION_UNIFORM<GRID<TV> > local_projection(local_mac_grid);
    smoke_tests->Map_Fine_To_Local_Boundary_For_Cell(local_mac_grid,local_face_velocities,cell_index);
    smoke_tests->Map_Fine_To_Local_Interior_For_Cell(local_mac_grid,local_face_velocities,cell_index,false);
    smoke_tests->Map_Fine_To_Local_Boundaries_For_Cell(local_mac_grid,local_projection.elliptic_solver->psi_N,cell_index);
    for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        local_projection.p(local_iterator.Cell_Index())=smoke_tests->projection.p(cell_index);}
    local_projection.elliptic_solver->Set_Neumann_Outer_Boundaries();
    local_projection.p*=dt;        
    local_projection.Make_Divergence_Free(local_face_velocities,dt,time);
    local_projection.p/=dt;
    smoke_tests->Map_Local_To_Fine_Interior_For_Cell(local_mac_grid,local_face_velocities,cell_index);
}
template class REFINEMENT_TASK<VECTOR<float,2> >;
template class REFINEMENT_TASK<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REFINEMENT_TASK<VECTOR<double,2> >;
template class REFINEMENT_TASK<VECTOR<double,3> >;
#endif
