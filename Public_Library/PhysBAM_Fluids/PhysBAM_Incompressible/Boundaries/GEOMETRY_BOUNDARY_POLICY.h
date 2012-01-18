//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GEOMETRY_BOUNDARY_POLICY__
#define __GEOMETRY_BOUNDARY_POLICY__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_FORWARD.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T_GRID>
struct GEOMETRY_BOUNDARY_POLICY
{
    typedef PhysBAM::BOUNDARY_REFLECTION_UNIFORM<T_GRID,typename T_GRID::SCALAR> BOUNDARY_REFLECTION;
    typedef PhysBAM::BOUNDARY_PHI_WATER<T_GRID> BOUNDARY_PHI_WATER; 
    typedef PhysBAM::BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP;
};
}

#endif
