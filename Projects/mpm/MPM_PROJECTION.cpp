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
    // mac_grid=sim.grid.Get_MAC_Grid();
    mac_grid.Initialize(sim.grid.numbers_of_cells+1,RANGE<TV>(sim.grid.domain.min_corner-sim.grid.dX*0.5,sim.grid.domain.max_corner+sim.grid.dX*0.5),true);
    face_velocities.Resize(mac_grid);
    LOG::cout<<"mac_grid.counts = "<<mac_grid.counts<<std::endl; // count of cell centers
    cell_dirichlet.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts-1));
    cell_neumann.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts-1));
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
        cell_dirichlet(sim.grid.Cell(sim.particles.X(p),0))=true;
}

//#####################################################################
// Function Identify_Neumann_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Identify_Neumann_Cells()
{
    // TODO
}

// //#####################################################################
// // Function Interpolate_Velocities_Nodes_To_Faces
// //#####################################################################
// template<class TV> void MPM_PROJECTION<TV>::
// Interpolate_Velocities_Nodes_To_Faces()
// {
//     for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
//         FACE_INDEX<TV::dimention> face_index=iterator.Full_Index();
//         face_velocities(face_index)=(T)0;
//         int axis=iterator.Axis();
//         TV_INT second_cell=iterator.Second_Cell_Index();
//         for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(second_cell,second_cell+1));it.Valid();it.Next())
//             face_velocities(face_index)+=sim.node_V(it.index)(axis);
//         face_velocities(face_index)/=TV::dimension*2-2;
//     }
//     // Some example code of using FACE_ITERATOR..
//     for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
//         face_velocities(iterator.Full_Index())=0;
//         RANGE<TV> dual_cell=iterator.Dual_Cell();
//         TV location=iterator.Location();
//         int axis=iterator.Axis();
//         TV_INT face_index=iterator.Face_Index();
//         for(int side=0;side<2;side++){
//             TV_INT cell_index=face_index+(1-side)*TV_INT::Axis_Vector(axis);
//         }
//     }
//     for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
//         int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
//         if(!psi_N.Component(axis)(face_index) && !(psi_D(first_cell) && psi_D(second_cell)))
//             face_velocities.Component(axis)(face_index)-=poisson->beta_face.Component(axis)(face_index)*(p(second_cell)-p(first_cell))*one_over_dx[axis];}
// }

//#####################################################################
template class MPM_PROJECTION<VECTOR<float,2> >;
template class MPM_PROJECTION<VECTOR<float,3> >;
template class MPM_PROJECTION<VECTOR<double,2> >;
template class MPM_PROJECTION<VECTOR<double,3> >;
}
