//#####################################################################
// Copyright 2004-2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_SOLID_WALL_SLIP_OUTFLOW
//#####################################################################
#ifndef __BOUNDARY_SOLID_WALL_SLIP_OUTFLOW__
#define __BOUNDARY_SOLID_WALL_SLIP_OUTFLOW__

#include <Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
namespace PhysBAM{

template<class TV>
class BOUNDARY_SOLID_WALL_SLIP_OUTFLOW:public BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV> BASE;
    using BASE::Constant_Extrapolation;using BASE::lower_threshold;using BASE::upper_threshold;using BASE::clamp_below;using BASE::clamp_above;

    BOUNDARY_SOLID_WALL_SLIP_OUTFLOW(const bool left_constant_extrapolation_input=false,const bool right_constant_extrapolation_input=false,
        const bool bottom_constant_extrapolation_input=false,const bool top_constant_extrapolation_input=false,
        const bool front_constant_extrapolation_input=false,const bool back_constant_extrapolation_input=false)
        :BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>(left_constant_extrapolation_input,right_constant_extrapolation_input,bottom_constant_extrapolation_input,top_constant_extrapolation_input,
            front_constant_extrapolation_input,back_constant_extrapolation_input)
    {}

//#####################################################################
    void Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,ARRAY<T,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells=3) const override;
//#####################################################################
};
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
//#####################################################################
}
#endif
