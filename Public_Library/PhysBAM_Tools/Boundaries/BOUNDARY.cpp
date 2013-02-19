//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
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
// Destructor
//#####################################################################
template<class TV,class T2> BOUNDARY<TV,T2>::
~BOUNDARY()
{
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    ARRAYS_ND_BASE<T2,TV_INT>::Put(u,u_ghost);
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=0;side<GRID<TV>::number_of_faces_per_cell;side++)Fill_Single_Ghost_Region(grid,u_ghost,side,regions(side));
}
//#####################################################################
// Function Fill_Ghost_Faces
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells)
{
    assert(grid.Is_MAC_Grid() && !clamp_below && !clamp_above);BOUNDARY<TV,T2> temp_boundary;
    if(use_fixed_boundary) temp_boundary.Set_Fixed_Boundary(true,fixed_boundary_value);
    for(int axis=0;axis<GRID<TV>::dimension;axis++)temp_boundary.Fill_Ghost_Cells(grid.Get_Face_Grid(axis),u.Component(axis),u_ghost.Component(axis),0,time,number_of_ghost_cells);
}
//#####################################################################
// Function Find_Ghost_Regions
//#####################################################################
static inline void Find_Ghost_Regions_Helper(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& domain,const RANGE<VECTOR<int,1> >& ghost)
{
    regions(0)=RANGE<VECTOR<int,1> >(VECTOR<int,1>(ghost.min_corner.x),VECTOR<int,1>(domain.min_corner.x)); // left
    regions(1)=RANGE<VECTOR<int,1> >(VECTOR<int,1>(domain.max_corner.x),VECTOR<int,1>(ghost.max_corner.x)); // right
}
static inline void Find_Ghost_Regions_Helper(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& domain,const RANGE<VECTOR<int,2> >& ghost)
{
    regions(0)=RANGE<VECTOR<int,2> >(VECTOR<int,2>(ghost.min_corner.x,domain.min_corner.y),VECTOR<int,2>(domain.min_corner.x,domain.max_corner.y)); // left
    regions(1)=RANGE<VECTOR<int,2> >(VECTOR<int,2>(domain.max_corner.x,domain.min_corner.y),VECTOR<int,2>(ghost.max_corner.x,domain.max_corner.y)); // right
    regions(2)=RANGE<VECTOR<int,2> >(VECTOR<int,2>(ghost.min_corner.x,ghost.min_corner.y),VECTOR<int,2>(ghost.max_corner.x,domain.min_corner.y)); // bottom
    regions(3)=RANGE<VECTOR<int,2> >(VECTOR<int,2>(ghost.min_corner.x,domain.max_corner.y),VECTOR<int,2>(ghost.max_corner.x,ghost.max_corner.y)); // top
}
static inline void Find_Ghost_Regions_Helper(ARRAY<RANGE<VECTOR<int,3> > >& regions,const RANGE<VECTOR<int,3> >& domain,const RANGE<VECTOR<int,3> >& ghost)
{
    regions(0)=RANGE<VECTOR<int,3> >(VECTOR<int,3>(ghost.min_corner.x,domain.min_corner.y,domain.min_corner.z),VECTOR<int,3>(domain.min_corner.x,domain.max_corner.y,domain.max_corner.z)); // left
    regions(1)=RANGE<VECTOR<int,3> >(VECTOR<int,3>(domain.max_corner.x,domain.min_corner.y,domain.min_corner.z),VECTOR<int,3>(ghost.max_corner.x,domain.max_corner.y,domain.max_corner.z)); // right
    regions(2)=RANGE<VECTOR<int,3> >(VECTOR<int,3>(ghost.min_corner.x,ghost.min_corner.y,domain.min_corner.z),VECTOR<int,3>(ghost.max_corner.x,domain.min_corner.y,domain.max_corner.z)); // bottom
    regions(3)=RANGE<VECTOR<int,3> >(VECTOR<int,3>(ghost.min_corner.x,domain.max_corner.y,domain.min_corner.z),VECTOR<int,3>(ghost.max_corner.x,ghost.max_corner.y,domain.max_corner.z)); // top
    regions(4)=RANGE<VECTOR<int,3> >(VECTOR<int,3>(ghost.min_corner.x,ghost.min_corner.y,ghost.min_corner.z),VECTOR<int,3>(ghost.max_corner.x,ghost.max_corner.y,domain.min_corner.z)); // front
    regions(5)=RANGE<VECTOR<int,3> >(VECTOR<int,3>(ghost.min_corner.x,ghost.min_corner.y,domain.max_corner.z),VECTOR<int,3>(ghost.max_corner.x,ghost.max_corner.y,ghost.max_corner.z)); // back
}
template<class TV,class T2> void BOUNDARY<TV,T2>::
Find_Ghost_Regions(const GRID<TV>& grid,ARRAY<RANGE<TV_INT> >& regions,const int ghost_cells) const
{
    RANGE<TV_INT> domain=grid.Domain_Indices(),ghost=domain.Thickened(ghost_cells);
    regions.Resize(GRID<TV>::number_of_faces_per_cell);
    Find_Ghost_Regions_Helper(regions,domain,ghost);
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Fill_Single_Ghost_Region_Threaded(RANGE<TV_INT>& region,const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const int side)
{
    Fill_Single_Ghost_Region(grid,u_ghost,side,region);
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Fill_Single_Ghost_Region(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const int side,const RANGE<TV_INT>& region) const
{
    int axis=side/2,boundary=side&1?region.Minimum_Corner()[axis]-1:region.Maximum_Corner()[axis];
    UNIFORM_GRID_ITERATOR_NODE<TV> iterator(grid,region);
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
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time)
{
}
//#####################################################################
// Function Apply_Boundary_Condition_Face
//#####################################################################
template<class TV,class T2> void BOUNDARY<TV,T2>::
Apply_Boundary_Condition_Face(const GRID<TV>& grid,ARRAY<T2,FACE_INDEX<TV::m> >& u,const T time)
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
template class BOUNDARY<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2> >;
template class BOUNDARY<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3> >;
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
template class BOUNDARY<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2> >;
template class BOUNDARY<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3> >;
template class BOUNDARY<VECTOR<float,1>,int>;
template class BOUNDARY<VECTOR<float,2>,int>;
template class BOUNDARY<VECTOR<float,3>,int>;
template class BOUNDARY<VECTOR<double,1>,int>;
template class BOUNDARY<VECTOR<double,2>,int>;
template class BOUNDARY<VECTOR<double,3>,int>;
}
