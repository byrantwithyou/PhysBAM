//#####################################################################
// Copyright 2013
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include "MPM_POISSON_SYSTEM.h"
#include "MPM_POISSON_VECTOR.h"
#include "MPMAC.h"
namespace PhysBAM{

//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MPMAC<TV>::
Initialize()
{
    mac_grid.Initialize(grid.numbers_of_cells,RANGE<TV>(grid.domain.min_corner,grid.domain.max_corner),true);
    for(int d=0;d<TV::m;d++) face_grid[d]=mac_grid.Get_Face_Grid(d);
    face_velocities.Resize(mac_grid);
    face_masses.Resize(mac_grid);
    face_momenta.Resize(mac_grid);
    cell_dirichlet.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    cell_neumann.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    div_u.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    pressure.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    for(int d=0;d<TV::m;d++){
        influence_corner[d].Resize(particles.number);
        weight[d].Resize(particles.number);
        for(int p=0;p<particles.number;p++)
            weight[d](p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));}
    min_mass=particles.mass.Min()*(T)1e-5;
    frame=0;
}

//#####################################################################
// Function Reinitialize
//#####################################################################
template<class TV> void MPMAC<TV>::
Reinitialize()
{
    face_velocities.Fill((T)0);
    face_masses.Fill((T)0);
    face_momenta.Fill((T)0);
    cell_dirichlet.Fill(false);
    cell_neumann.Fill(false);
    div_u.Fill((T)0);
    pressure.Fill((T)0);
}

//#####################################################################
// Function Weights
//#####################################################################
template<class TV> void MPMAC<TV>::
Weights()
{
#pragma omp parallel for
    for(int p=0;p<particles.number;p++)
        for(int axis=0;axis<TV::m;axis++)
            grid_basis_function.Build_Weights_Exact(particles.X(p),face_grid[axis],influence_corner[axis](p),weight[axis](p));
}

//#####################################################################
// Function Rasterize();
//#####################################################################
template<class TV> void MPMAC<TV>::
Rasterize()
{
    for(int p=0;p<particles.number;p++){
        for(int axis=0;axis<TV::m;axis++){
        }}
        
}

//#####################################################################
template class MPMAC<VECTOR<float,2> >;
template class MPMAC<VECTOR<float,3> >;
template class MPMAC<VECTOR<double,2> >;
template class MPMAC<VECTOR<double,3> >;
}
