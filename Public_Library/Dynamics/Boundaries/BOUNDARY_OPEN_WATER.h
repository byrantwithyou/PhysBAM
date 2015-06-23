//#####################################################################
// Copyright 2002-2006, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OPEN_WATER
//#####################################################################
#ifndef __BOUNDARY_OPEN_WATER__
#define __BOUNDARY_OPEN_WATER__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Incompressible/Boundaries/BOUNDARY_FORWARD.h>
namespace PhysBAM{

template<class TV>
class BOUNDARY_OPEN_WATER:public BOUNDARY<TV,typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
public:
    typedef BOUNDARY<TV,T> BASE;
    using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;using BASE::Find_Ghost_Regions;

    ARRAY<bool> open_boundary;
    T attenuate_inflow;
    BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV> boundary_mac_grid_solid_wall_slip;

    BOUNDARY_OPEN_WATER(const T attenuate_inflow_input=T(1),const bool left_open_boundary_input=false,const bool right_open_boundary_input=false,const bool bottom_open_boundary_input=false,const bool top_open_boundary_input=false,
        const bool front_open_boundary_input=false,const bool back_open_boundary_input=false)
        :open_boundary(6,true),boundary_mac_grid_solid_wall_slip()
    {
        attenuate_inflow=attenuate_inflow_input;
        open_boundary(0)=left_open_boundary_input;
        open_boundary(1)=right_open_boundary_input;
        open_boundary(2)=bottom_open_boundary_input;
        open_boundary(3)=top_open_boundary_input;
        open_boundary(4)=front_open_boundary_input;
        open_boundary(5)=back_open_boundary_input;
    }

    ~BOUNDARY_OPEN_WATER(){}

public:

//#####################################################################
    void Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,ARRAY<T,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells=3) const override;
//#####################################################################
};
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
#endif
