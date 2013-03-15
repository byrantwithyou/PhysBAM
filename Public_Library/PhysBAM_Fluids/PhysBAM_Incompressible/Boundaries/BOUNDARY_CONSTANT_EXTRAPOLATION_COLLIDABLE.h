//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE
//#####################################################################
#ifndef __BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE__
#define __BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE__

#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename GRID<TV>::RIGID_BODY_LIST T_RIGID_BODY_LIST;
public:
    T_RIGID_BODY_LIST& body_list;

    BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE(T_RIGID_BODY_LIST& body_list_input)
        :body_list(body_list_input)
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<T2,TV_INT>& u,ARRAY<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<T2,TV_INT>& u,const T time) const PHYSBAM_OVERRIDE {} // do nothing
    void Collision_Aware_Extrapolate(const GRID<TV>& grid,ARRAY<T2,TV_INT>& u_ghost,const TV_INT& source_index,const TV_INT& ghost_index,const TV& direction,const T ray_length);
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<T2,TV_INT>& u,ARRAY<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    ARRAY<T2,TV_INT>::Put(u,u_ghost); // interior
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=0;axis<GRID<TV>::dimension;axis++)for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side,outward_sign=axis_side?-1:1;
        TV direction=outward_sign*TV::Axis_Vector(axis);
        T signed_dx=outward_sign*grid.DX()[axis];
        int boundary=Boundary(side,regions(side));
        for(NODE_ITERATOR<TV> iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            TV_INT boundary_node=node;boundary_node[axis]=boundary;T ray_length=signed_dx*(node[axis]-boundary);
            Collision_Aware_Extrapolate(grid,u_ghost,boundary_node,node,direction,ray_length);}}
}
//#####################################################################
// Function Collision_Aware_Extrapolate
//#####################################################################
template<class TV,class T2> void BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE<TV,T2>::
Collision_Aware_Extrapolate(const GRID<TV>& grid,ARRAY<T2,TV_INT>& u_ghost,const TV_INT& source_index,const TV_INT& ghost_index,const TV& direction,const T ray_length)
{
    int body_id;
    RAY<TV> ray=RAY<TV>(grid.X(source_index),direction,true);ray.semi_infinite=false;ray.t_max=ray_length;
    PHYSBAM_NOT_IMPLEMENTED(); // update to use fluid collision body list
    if(body_list.Intersection_With_Any_Simplicial_Object(ray,body_id)) u_ghost(ghost_index)=0; // ZERO MIGHT NOT BE GREAT FOR EVERYTHING!
    else u_ghost(ghost_index)=u_ghost(source_index);
}
//#####################################################################
}
#endif
