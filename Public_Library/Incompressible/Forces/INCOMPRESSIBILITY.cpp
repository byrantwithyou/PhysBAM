//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <Tools/Log/LOG.h>
#include <Incompressible/Forces/INCOMPRESSIBILITY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INCOMPRESSIBILITY<TV>::
INCOMPRESSIBILITY(PROJECTION_UNIFORM<TV>& projection_input)
    :projection(projection_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INCOMPRESSIBILITY<TV>::
~INCOMPRESSIBILITY()
{
}
//#####################################################################
// Function Add_Implicit_Forces
//#####################################################################
template<class TV> void INCOMPRESSIBILITY<TV>::
Add_Implicit_Forces_Projection(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    projection.Make_Divergence_Free(face_velocities,dt,time);
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class TV> void INCOMPRESSIBILITY<TV>::
Initialize_Grids(const GRID<TV>& grid)
{
    projection.Initialize_Grid(grid);
}
//#####################################################################
namespace PhysBAM{
template class INCOMPRESSIBILITY<VECTOR<float,2> >;
template class INCOMPRESSIBILITY<VECTOR<float,3> >;
template class INCOMPRESSIBILITY<VECTOR<double,2> >;
template class INCOMPRESSIBILITY<VECTOR<double,3> >;
}
