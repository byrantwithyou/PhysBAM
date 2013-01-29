//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MAC_GRID_PERIODIC
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
using namespace PhysBAM;
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_MAC_GRID_PERIODIC<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    ARRAYS_ND_BASE<T2,TV_INT>::Put(u,u_ghost); // interior
    TV_INT periods=grid.Domain_Indices().Maximum_Corner();
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=0;axis<GRID<TV>::dimension;axis++)for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side;
        TV_INT period=(axis_side==0?1:-1)*periods[axis]*TV_INT::Axis_Vector(axis);
        for(NODE_ITERATOR iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            u_ghost(node)=u_ghost(node+period);}}
}
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class TV,class T2> void BOUNDARY_MAC_GRID_PERIODIC<TV,T2>::
Fill_Ghost_Cells_Face(const GRID<TV>& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells)
{
    assert(grid.Is_MAC_Grid());
    for(int axis=0;axis<GRID<TV>::dimension;axis++){
        //Fill_Ghost_Cells(grid.Get_Face_Grid(axis,u.Component(axis),u_ghost.Component(axis),0,time,number_of_ghost_cells);
        const ARRAYS_ND_BASE<T2,TV_INT>& u_axis=u.Component(axis);
        ARRAYS_ND_BASE<T2,TV_INT>& u_ghost_axis=u_ghost.Component(axis);
        ARRAYS_ND_BASE<T2,TV_INT>::Put(u_axis,u_ghost_axis); // interior
        GRID<TV> face_grid=grid.Get_Face_Grid(axis);
        TV_INT periods=grid.Domain_Indices().Maximum_Corner();
        ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int face_axis=0;face_axis<GRID<TV>::dimension;face_axis++)for(int axis_side=0;axis_side<2;axis_side++){
            int side=2*face_axis+axis_side;
            TV_INT period=(axis_side==0?1:-1)*periods[face_axis]*TV_INT::Axis_Vector(face_axis);
            for(NODE_ITERATOR iterator(face_grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
                u_ghost_axis(node)=u_ghost_axis(node+period);}}}
}
//#####################################################################
namespace PhysBAM{
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,2>,VECTOR<float,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,2>,float>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,3>,float>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,2>,VECTOR<double,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,2>,double>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,3>,double>;
}
