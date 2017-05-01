//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> BOUNDARY<TV,T2>::
BOUNDARY()
{
    Set_Constant_Extrapolation();
    Set_Fixed_Boundary(false);
    Limit_Minimum_Boundary_Value(false);
    Limit_Maximum_Boundary_Value(false);
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    ARRAYS_ND_BASE<T2,TV_INT>::Put(u,u_ghost);
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=0;side<GRID<TV>::number_of_faces_per_cell;side++)Fill_Single_Ghost_Region(grid,u_ghost,side,regions(side));
}
//#####################################################################
// Function Fill_Ghost_Faces
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells) const
{
    assert(grid.Is_MAC_Grid() && !clamp_below && !clamp_above);BOUNDARY<TV,T2> temp_boundary;
    if(use_fixed_boundary) temp_boundary.Set_Fixed_Boundary(true,fixed_boundary_value);
    for(int axis=0;axis<TV::m;axis++)temp_boundary.Fill_Ghost_Cells(grid.Get_Face_Grid(axis),u.Component(axis),u_ghost.Component(axis),0,time,number_of_ghost_cells);
}
//#####################################################################
// Function Find_Ghost_Regions
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Find_Ghost_Regions(const GRID<TV>& grid,VECTOR<RANGE<TV_INT>,2*TV::m>& regions,const int ghost_cells)
{
    RANGE<TV_INT> inner=grid.Domain_Indices(),ghost=grid.Domain_Indices(ghost_cells),current=inner;
    for(int axis=0;axis<TV::m;axis++){
        current.min_corner(axis)=ghost.min_corner(axis);
        current.max_corner(axis)=inner.min_corner(axis);
        regions(2*axis)=current;
        current.min_corner(axis)=inner.max_corner(axis);
        current.max_corner(axis)=ghost.max_corner(axis);
        regions(2*axis+1)=current;
        current.min_corner(axis)=ghost.min_corner(axis);}
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Fill_Single_Ghost_Region(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const int side,const RANGE<TV_INT>& region) const
{
    int axis=side/2,boundary=side&1?region.Minimum_Corner()[axis]-1:region.Maximum_Corner()[axis];
    NODE_ITERATOR<TV> iterator(grid,region);
    if(use_fixed_boundary) for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
        u_ghost(node)=fixed_boundary_value;}
    else if(clamp_below&&clamp_above) for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary;
        u_ghost(node)=clamp(u_ghost(boundary_node),lower_threshold,upper_threshold);}
    else if(clamp_below) for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary;
        u_ghost(node)=clamp_min(u_ghost(boundary_node),lower_threshold);}
    else if(clamp_above) for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary;
        u_ghost(node)=clamp_max(u_ghost(boundary_node),upper_threshold);}
    else for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary;
        u_ghost(node)=u_ghost(boundary_node);}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const
{
}
//#####################################################################
// Function Apply_Boundary_Condition_Face
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Apply_Boundary_Condition_Face(const GRID<TV>& grid,ARRAY<T2,FACE_INDEX<TV::m> >& u,const T time) const
{
}
//#####################################################################
// Function Apply_Boundary_Condition_Single_Side
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Apply_Boundary_Condition_Single_Side(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const int side,const T time) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Fill_Single_Ghost_Region(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
namespace PhysBAM{
template class BOUNDARY<VECTOR<float,1>,float>;
template class BOUNDARY<VECTOR<float,1>,VECTOR<float,1> >;
template class BOUNDARY<VECTOR<float,1>,VECTOR<float,2> >;
template class BOUNDARY<VECTOR<float,1>,VECTOR<float,3> >;
template class BOUNDARY<VECTOR<float,1>,VECTOR<float,4> >;
template class BOUNDARY<VECTOR<float,2>,float>;
template class BOUNDARY<VECTOR<float,2>,VECTOR<float,2> >;
template class BOUNDARY<VECTOR<float,2>,VECTOR<float,3> >;
template class BOUNDARY<VECTOR<float,2>,VECTOR<float,4> >;
template class BOUNDARY<VECTOR<float,2>,VECTOR<float,5> >;
template class BOUNDARY<VECTOR<float,3>,float>;
template class BOUNDARY<VECTOR<float,3>,VECTOR<float,3> >;
template class BOUNDARY<VECTOR<float,3>,VECTOR<float,5> >;
template class BOUNDARY<VECTOR<float,3>,VECTOR<float,6> >;
template class BOUNDARY<VECTOR<float,3>,bool>;
template class BOUNDARY<VECTOR<float,1>,SYMMETRIC_MATRIX<float,1> >;
template class BOUNDARY<VECTOR<float,1>,MATRIX<float,1> >;
template class BOUNDARY<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2> >;
template class BOUNDARY<VECTOR<float,2>,MATRIX<float,2> >;
template class BOUNDARY<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3> >;
template class BOUNDARY<VECTOR<float,3>,MATRIX<float,3> >;
template class BOUNDARY<VECTOR<double,1>,double>;
template class BOUNDARY<VECTOR<double,1>,VECTOR<double,1> >;
template class BOUNDARY<VECTOR<double,1>,VECTOR<double,2> >;
template class BOUNDARY<VECTOR<double,1>,VECTOR<double,3> >;
template class BOUNDARY<VECTOR<double,1>,VECTOR<double,4> >;
template class BOUNDARY<VECTOR<double,2>,double>;
template class BOUNDARY<VECTOR<double,2>,VECTOR<double,2> >;
template class BOUNDARY<VECTOR<double,2>,VECTOR<double,3> >;
template class BOUNDARY<VECTOR<double,2>,VECTOR<double,4> >;
template class BOUNDARY<VECTOR<double,2>,VECTOR<double,5> >;
template class BOUNDARY<VECTOR<double,3>,double>;
template class BOUNDARY<VECTOR<double,3>,VECTOR<double,3> >;
template class BOUNDARY<VECTOR<double,3>,VECTOR<double,5> >;
template class BOUNDARY<VECTOR<double,3>,VECTOR<double,6> >;
template class BOUNDARY<VECTOR<double,3>,bool>;
template class BOUNDARY<VECTOR<double,1>,SYMMETRIC_MATRIX<double,1> >;
template class BOUNDARY<VECTOR<double,1>,MATRIX<double,1> >;
template class BOUNDARY<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2> >;
template class BOUNDARY<VECTOR<double,2>,MATRIX<double,2> >;
template class BOUNDARY<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3> >;
template class BOUNDARY<VECTOR<double,3>,MATRIX<double,3> >;
template class BOUNDARY<VECTOR<float,1>,int>;
template class BOUNDARY<VECTOR<float,2>,int>;
template class BOUNDARY<VECTOR<float,3>,int>;
template class BOUNDARY<VECTOR<double,1>,int>;
template class BOUNDARY<VECTOR<double,2>,int>;
template class BOUNDARY<VECTOR<double,3>,int>;
}
