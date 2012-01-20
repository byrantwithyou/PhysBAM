//#####################################################################
// Copyright 2009, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Function Use_Maccormack_Advection
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class T_GRID> void LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Use_Maccormack_Advection(const T_ARRAYS_BOOL& cell_mask)
{
    advection_maccormack=new ADVECTION_MACCORMACK_UNIFORM<T_GRID,SCALAR,ADVECTION<T_GRID,SCALAR> >(*advection,0,&cell_mask,0);
    Set_Custom_Advection(*advection_maccormack);
}
//#####################################################################
// Function Approximate_Negative_Material
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class T_GRID> typename T_GRID::VECTOR_T::SCALAR LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Approximate_Negative_Material(const SCALAR interface_thickness,const SCALAR time) const
{
    T_GRID& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    T_GRID node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    SCALAR interface_half_width=interface_thickness*grid.dX.Max()/2,volume=0;
    for(NODE_ITERATOR iterator(node_grid);iterator.Valid();iterator.Next()) volume+=LEVELSET_UTILITIES<SCALAR>::Heaviside(-phi(iterator.Node_Index()),interface_half_width);
    return volume*grid.Cell_Size();
}
//#####################################################################
// Function Approximate_Positive_Material
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class T_GRID> typename T_GRID::VECTOR_T::SCALAR LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Approximate_Positive_Material(const SCALAR interface_thickness,const SCALAR time) const
{
    T_GRID& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    T_GRID node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    SCALAR interface_half_width=interface_thickness*grid.dX.Max()/2,volume=0;
    for(NODE_ITERATOR iterator(node_grid);iterator.Valid();iterator.Next()) volume+=LEVELSET_UTILITIES<SCALAR>::Heaviside(phi(iterator.Node_Index()),interface_half_width);
    return volume*grid.Cell_Size();
}

template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<float,1> > >;
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<float,2> > >;
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<double,1> > >;
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<double,2> > >;
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
