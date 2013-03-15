//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<TV>::
INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS(const VECTOR<VECTOR<bool,2>,TV::dimension>& walls_input,const VECTOR<VECTOR<bool,2>,TV::dimension>& mpi_boundary_input)
    :walls(walls_input),mpi_boundary(mpi_boundary_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<TV>::
~INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS()
{}
//#####################################################################
// Function Update_Boundary_Conditions
//#####################################################################
template<class TV> void INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<TV>::
Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
    ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time)
{   
    for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side;
        if(mpi_boundary(axis)(axis_side)){
            for(FACE_ITERATOR<TV> iterator(grid,1,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next())
                if(axis!=iterator.Axis()) psi_N(iterator.Full_Index())=true;}
        else if(walls(axis)(axis_side)){
            for(CELL_ITERATOR<TV> iterator(grid,1,GRID<TV>::GHOST_REGION,side);iterator.Valid();iterator.Next()){
                psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=0;}
            for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                psi_N(iterator.Full_Index())=true;face_velocities(iterator.Axis(),iterator.Face_Index())=0;}}
        else
            for(CELL_ITERATOR<TV> iterator(grid,1,GRID<TV>::GHOST_REGION,side);iterator.Valid();iterator.Next()){
                psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=0;}}
}
namespace PhysBAM{
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<VECTOR<float,1> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<VECTOR<float,2> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<VECTOR<float,3> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<VECTOR<double,1> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<VECTOR<double,2> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<VECTOR<double,3> >;
}
