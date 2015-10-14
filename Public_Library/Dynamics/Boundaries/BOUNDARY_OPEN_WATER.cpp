//#####################################################################
// Copyright 2002-2006, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Incompressible/Boundaries/BOUNDARY_FORWARD.h>
#include <Dynamics/Boundaries/BOUNDARY_OPEN_WATER.h>
namespace PhysBAM{
//#####################################################################
// Function Fill_Ghost_Faces
//#####################################################################
template<class TV> void BOUNDARY_OPEN_WATER<TV>::
Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,ARRAY<T,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells) const
{
    assert(grid.Is_MAC_Grid());
    ARRAY<T,FACE_INDEX<TV::m> >::Put(u,u_ghost); // interior
    for(int face_axis=0;face_axis<TV::m;face_axis++){
        GRID<TV> face_grid=grid.Get_Face_Grid(face_axis);
        T_ARRAYS_BASE& u_ghost_component=u_ghost.Component(face_axis);
        ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int side=0;side<GRID<TV>::number_of_faces_per_cell;side++){
            RANGE<TV_INT>& region=regions(side);
            int axis=side/2,boundary=side&1?region.Minimum_Corner()[axis]-1:region.Maximum_Corner()[axis]+1;
            if(open_boundary(side)) for(NODE_ITERATOR<TV> iterator(face_grid,region);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary; 
                if(side&1) if(u_ghost_component(boundary_node)<0) u_ghost_component(node)=attenuate_inflow*u_ghost_component(boundary_node);else u_ghost_component(node)=u_ghost_component(boundary_node);
                else if(u_ghost_component(boundary_node)>0) u_ghost_component(node)=attenuate_inflow*u_ghost_component(boundary_node);else u_ghost_component(node)=u_ghost_component(boundary_node);}
            else if(Constant_Extrapolation(side)) Fill_Single_Ghost_Region(face_grid,u_ghost_component,side,regions(side));
            else boundary_mac_grid_solid_wall_slip.Reflect_Single_Ghost_Region(face_axis,face_grid,u_ghost_component,side,regions(side));}}
}
}
