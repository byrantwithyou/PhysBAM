//#####################################################################
// Copyright 2004-2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
#include <Incompressible/Boundaries/BOUNDARY_SOLID_WALL_SLIP_OUTFLOW.h>
namespace PhysBAM{
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV> void BOUNDARY_SOLID_WALL_SLIP_OUTFLOW<TV>::
Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,ARRAY<T,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells) const
{
    assert(grid.Is_MAC_Grid());
    ARRAY<T,FACE_INDEX<TV::m> >::Put(u,u_ghost); // interior
    lower_threshold=upper_threshold=0;
    for(int face_axis=0;face_axis<TV::m;face_axis++){
        GRID<TV> face_grid=grid.Get_Face_Grid(face_axis);
        ARRAY<T,TV_INT>& u_ghost_component=u_ghost.Component(face_axis);
        ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int side=0;side<GRID<TV>::number_of_faces_per_cell;side++){
            if(Constant_Extrapolation(side)){
                clamp_above=side&1;clamp_below=!clamp_above;
                Fill_Single_Ghost_Region(face_grid,u_ghost_component,side,regions(side));}
            else Reflect_Single_Ghost_Region(face_axis,face_grid,u_ghost_component,side,regions(side));}}
}
}
