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
template<class T_GRID> INCOMPRESSIBILITY<T_GRID>::
INCOMPRESSIBILITY(PROJECTION_UNIFORM<T_GRID>& projection_input)
    :projection(projection_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBILITY<T_GRID>::
~INCOMPRESSIBILITY()
{
}
//#####################################################################
// Function Add_Implicit_Forces
//#####################################################################
template<class T_GRID> void INCOMPRESSIBILITY<T_GRID>::
Add_Implicit_Forces_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    projection.Make_Divergence_Free(face_velocities,dt,time);
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class T_GRID> void INCOMPRESSIBILITY<T_GRID>::
Initialize_Grids(const T_GRID& grid)
{
    projection.Initialize_Grid(grid);
}
//#####################################################################
namespace PhysBAM{
template class INCOMPRESSIBILITY<GRID<VECTOR<float,2> > >;
template class INCOMPRESSIBILITY<GRID<VECTOR<float,3> > >;
template class INCOMPRESSIBILITY<GRID<VECTOR<double,2> > >;
template class INCOMPRESSIBILITY<GRID<VECTOR<double,3> > >;
}
