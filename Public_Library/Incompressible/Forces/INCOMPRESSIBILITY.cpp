//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
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
Add_Implicit_Forces_Projection(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
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
