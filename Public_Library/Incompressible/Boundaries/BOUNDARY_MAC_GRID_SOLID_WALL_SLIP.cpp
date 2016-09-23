//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
using namespace PhysBAM;
template<class TV> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>::
BOUNDARY_MAC_GRID_SOLID_WALL_SLIP(const TV_SIDES& constant_extrapolation)
    :phi(0)
{
    Set_Constant_Extrapolation(constant_extrapolation);
}
template<class TV> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>::
~BOUNDARY_MAC_GRID_SOLID_WALL_SLIP()
{
}
//#####################################################################
// Function Fill_Ghost_Faces
//#####################################################################
template<class TV> void BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>::
Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,ARRAY<T,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells) const
{
    assert(grid.Is_MAC_Grid());
    ARRAY<T,FACE_INDEX<TV::m> >::Put(u,u_ghost); // interior
    for(int face_axis=0;face_axis<TV::m;face_axis++){
        GRID<TV> face_grid=grid.Get_Face_Grid(face_axis);
        T_ARRAYS_BASE& u_ghost_component=u_ghost.Component(face_axis);
        VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int side=0;side<GRID<TV>::number_of_faces_per_cell;side++){
            if(Constant_Extrapolation(side)) Fill_Single_Ghost_Region(face_grid,u_ghost_component,side,regions(side));
            else Reflect_Single_Ghost_Region(face_axis,face_grid,u_ghost_component,side,regions(side));}}
}
//#####################################################################
// Function Reflect_Single_Ghost_Region
//#####################################################################
template<class TV> void BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>::
Reflect_Single_Ghost_Region(const int face_axis,const GRID<TV>& face_grid,T_ARRAYS_BASE& u_ghost_component,const int side,const RANGE<TV_INT>& region) const
{
    int axis=side/2,axis_side=side&1;
    int boundary=Boundary(side,region),reflection_times_two,flip;
    if(face_axis==axis){reflection_times_two=2*boundary;flip=-1;}
    else{reflection_times_two=2*boundary+(axis_side==0?-1:1);flip=1;}
    for(NODE_ITERATOR<TV> iterator(face_grid,region);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
        TV_INT reflected_node=node;reflected_node[axis]=reflection_times_two-node[axis];
        u_ghost_component(node)=flip*u_ghost_component(reflected_node);}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV> void BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>::
Apply_Boundary_Condition_Face(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& u,const T time)  const
{
    assert(grid.Is_MAC_Grid());
    for(int side=0;side<GRID<TV>::number_of_faces_per_cell;side++)
        if(!Constant_Extrapolation(side)) Zero_Single_Boundary_Side(grid,u,side);
}
//#####################################################################
// Function Zero_Single_Boundary_Side
//#####################################################################
template<class TV> void BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>::
Zero_Single_Boundary_Side(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& u,const int side) const
{
    int axis=side/2,axis_side=side&1;
    FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,side);
    if(phi){
        TV_INT interior_cell_offset=axis_side==0?TV_INT():-TV_INT::Axis_Vector(axis);
        for(;iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();
            if((*phi)(face+interior_cell_offset)) u.Component(axis)(face)=0;}}
    else for(;iterator.Valid();iterator.Next())u.Component(axis)(iterator.Face_Index())=0;
}
//#####################################################################
namespace PhysBAM{
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<VECTOR<float,1> >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<VECTOR<float,2> >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<VECTOR<float,3> >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<VECTOR<double,1> >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<VECTOR<double,2> >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<VECTOR<double,3> >;
}
