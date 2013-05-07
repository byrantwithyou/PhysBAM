//#####################################################################
// Copyright 2013
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include "MPM_PROJECTION.h"
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PROJECTION<TV>::
MPM_PROJECTION(MPM_SIMULATION<TV>& sim_in)
    :sim(sim_in)
{
    mac_grid.Initialize(sim.grid.numbers_of_cells+1,RANGE<TV>(sim.grid.domain.min_corner-sim.grid.dX*0.5,sim.grid.domain.max_corner+sim.grid.dX*0.5),true);
    face_velocities.Resize(mac_grid);
    LOG::cout<<"mac_grid.counts = "<<mac_grid.counts<<std::endl; // count of cell centers
    cell_dirichlet.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    LOG::cout<<"number of cells via mac_grid = "<<mac_grid.numbers_of_cells<<std::endl;
    LOG::cout<<"number of cells via cell_dirichlet = "<<cell_dirichlet.Size()<<std::endl;
    cell_neumann.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PROJECTION<TV>::
~MPM_PROJECTION()
{}

//#####################################################################
// Function Reinitialize
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Reinitialize()
{
    face_velocities.Fill((T)0);
    cell_dirichlet.Fill(false);
    cell_neumann.Fill(false);
}

//#####################################################################
// Function Identify_Dirichlet_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Identify_Dirichlet_Cells()
{
    // Cells without any particles are dirichlet cells
    // (zero pressure, extrapolated face velocity).
    cell_dirichlet.Fill(false);
#pragma omp parallel for
    for(int p=0;p<sim.particles.number;p++)
        cell_dirichlet(mac_grid.Cell(sim.particles.X(p),0))=true;
}

//#####################################################################
// Function Identify_Neumann_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Identify_Neumann_Cells()
{
    // TODO
}

//#####################################################################
// Function Interpolate_Velocities_To_Faces
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Interpolate_Velocities_To_Faces()
{
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face_index=iterator.Full_Index();
        int axis=iterator.Axis(); // the axis of first_cell -> second_cell
        TV_INT first_cell=iterator.First_Cell_Index();
        TV_INT second_cell=iterator.Second_Cell_Index();        
        if(first_cell(axis)>=0&&second_cell(axis)<mac_grid.counts(axis)) // only deal with non-boundary faces
            face_velocities(face_index)=0.5*(sim.node_V(first_cell)(axis)+sim.node_V(second_cell)(axis));
        LOG::cout<<"axis: "<<axis<<"first_cell: "<<first_cell<<"second_cell: "<<second_cell<<"v: "<<face_velocities(face_index)<<std::endl;}
}

//#####################################################################
template class MPM_PROJECTION<VECTOR<float,2> >;
template class MPM_PROJECTION<VECTOR<float,3> >;
template class MPM_PROJECTION<VECTOR<double,2> >;
template class MPM_PROJECTION<VECTOR<double,3> >;
}
