//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Dynamics/Coupled_Evolution/SURFACE_TENSION_BOUNDARY_CONDITION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_TENSION_BOUNDARY_CONDITION<TV>::
SURFACE_TENSION_BOUNDARY_CONDITION(const LEVELSET<TV>& levelset_input,T surface_tension_coefficient_input)
    :levelset(levelset_input),surface_tension_coefficient(surface_tension_coefficient_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SURFACE_TENSION_BOUNDARY_CONDITION<TV>::
~SURFACE_TENSION_BOUNDARY_CONDITION()
{
}
//#####################################################################
// Function Update_Boundary_Conditions
//#####################################################################
template<class TV> void SURFACE_TENSION_BOUNDARY_CONDITION<TV>::
Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& p,
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time)
{
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
        if(levelset.phi(iterator.Cell_Index())>0){
            PHYSBAM_ASSERT(psi_D(iterator.Cell_Index()));
            p(iterator.Cell_Index())-=surface_tension_coefficient*levelset.Compute_Curvature(levelset.phi,iterator.index);}
}
namespace PhysBAM{
template class SURFACE_TENSION_BOUNDARY_CONDITION<VECTOR<float,1> >;
template class SURFACE_TENSION_BOUNDARY_CONDITION<VECTOR<float,2> >;
template class SURFACE_TENSION_BOUNDARY_CONDITION<VECTOR<float,3> >;
template class SURFACE_TENSION_BOUNDARY_CONDITION<VECTOR<double,1> >;
template class SURFACE_TENSION_BOUNDARY_CONDITION<VECTOR<double,2> >;
template class SURFACE_TENSION_BOUNDARY_CONDITION<VECTOR<double,3> >;
}
