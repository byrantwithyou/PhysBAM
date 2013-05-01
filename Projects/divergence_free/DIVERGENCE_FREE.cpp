//#####################################################################
// Copyright 2013
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include "DIVERGENCE_FREE.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DIVERGENCE_FREE<TV>::
DIVERGENCE_FREE(const GRID<TV>& grid_in,const GRID<TV>& mac_grid_in,ARRAY<TV,TV_INT>& v_in,const ARRAY<T,TV_INT>& m_in,const T m_eps_in)
    :grid(grid_in),mac_grid(mac_grid_in),v(v_in),m(m_in),m_eps(m_eps_in),projection(mac_grid,false,true)
{
    // initializations
    projection.Initialize_Grid(mac_grid);
    face_velocities.Resize(mac_grid);
    projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DIVERGENCE_FREE<TV>::
~DIVERGENCE_FREE()
{}
//#####################################################################
// Function Interpolate_Velocities_Nodes_To_Faces
//#####################################################################
template<class TV> void DIVERGENCE_FREE<TV>::
Interpolate_Velocities_Nodes_To_Faces()
{
//    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
}
//#####################################################################
// Function Identify_Dirichlet_Cells
//#####################################################################
template<class TV> void DIVERGENCE_FREE<TV>::
Identify_Dirichlet_Cells()
{
    // Dirichlet
    TV_INT cell; // pressure cell index
    projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;
}
//#####################################################################
// Function Identify_Neumann_Cells
//#####################################################################
template<class TV> void DIVERGENCE_FREE<TV>::
Identify_Neumann_Cells()
{
    // Naumann
    // for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
    //     if(source.Lazy_Inside(iterator.Location())){
    //         projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
    //         if(iterator.Axis()==1)face_velocities(iterator.Full_Index())=1;
    //         else face_velocities(iterator.Full_Index())=0;}}
}
//#####################################################################
// Function Set_Up_Solver
//#####################################################################
template<class TV> void DIVERGENCE_FREE<TV>::
Set_Up_Solver()
{
    // CG solver paramaters
    projection.elliptic_solver->Set_Relative_Tolerance(1e-9);
    projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    projection.elliptic_solver->pcg.cg_restart_iterations=40;
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class TV> void DIVERGENCE_FREE<TV>::
Make_Divergence_Free(T dt,T time)
{
    projection.p*=dt; // rescale pressure for guess
    projection.Make_Divergence_Free(face_velocities,dt,time);
    projection.p*=(1/dt); // unscale pressure
}
//#####################################################################
template class DIVERGENCE_FREE<VECTOR<float,2> >;
template class DIVERGENCE_FREE<VECTOR<float,3> >;
template class DIVERGENCE_FREE<VECTOR<double,2> >;
template class DIVERGENCE_FREE<VECTOR<double,3> >;
}
