//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace RASTERIZATION
//##################################################################### 
#include <Tools/Grids_Uniform/BLOCK_UNIFORM.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Incompressible/Collisions_And_Interactions/OBJECTS_IN_CELL.h>
#include <Incompressible/Collisions_And_Interactions/RIGID_BODY_RASTERIZATION_HELPER.h>
namespace PhysBAM{

template<class TV> class GRID;

namespace RASTERIZATION{
//#####################################################################
// Function Rasterize_Object
//#####################################################################
template<class TV> void Rasterize_Object(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID& id)
{
    Rasterize_Object_Generic(collision_geometry,grid,objects,id);
}
//#####################################################################
// Function Compute_Occupied_Blocks
//#####################################################################
template<class T,class TV> void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,
    ARRAY<bool,VECTOR<int,TV::m> >& occupied,const bool with_body_motion,const T& extra_thickness,const T& body_thickness_factor)
{
    typedef VECTOR<int,TV::m> TV_INT;
    if(collision_geometry.Number_Of_Simplices()) Compute_Occupied_Blocks_Generic(collision_geometry,grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
    else{
        for(NODE_ITERATOR<TV> node_iterator(grid,2);node_iterator.Valid();node_iterator.Next()){
            TV_INT block_index=node_iterator.Node_Index();BLOCK_UNIFORM<TV> block(grid,block_index);
            for(int cell_index=0;cell_index<GRID<TV>::number_of_cells_per_block;cell_index++){TV_INT cell=block.Cell(cell_index);
                if(collision_geometry.Implicit_Geometry_Extended_Value(grid.X(cell))<=extra_thickness){occupied(block_index)=true;break;}}}}
}
//#####################################################################
// Function Rasterize_Box_Onto_Blocks
//#####################################################################
template<class TV> void Rasterize_Box_Onto_Blocks(const GRID<TV>& grid,ARRAY<bool,VECTOR<int,TV::m> >& occupied,const RANGE<TV>& box)
{
    TV DX_over_two=(typename TV::SCALAR).5*grid.dX;
    for(CELL_ITERATOR<TV> iterator(grid,grid.Clamp_To_Cell(box.Translated(-DX_over_two),3));iterator.Valid();iterator.Next())
        occupied(iterator.Cell_Index()+1)=true;
}
//#####################################################################
// Function Rasterize_Box
//#####################################################################
template<class TV> void Rasterize_Box(const GRID<TV>& grid,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects_in_cell,const RANGE<TV>& box,const COLLISION_GEOMETRY_ID id)
{
    for(CELL_ITERATOR<TV> iterator(grid,grid.Clamp_To_Cell(box,3));iterator.Valid();iterator.Next())
        objects_in_cell.Add_Object_To_Cell(iterator.Cell_Index(),id);
}
//#####################################################################
template void Compute_Occupied_Blocks<double,VECTOR<double,1> >(COLLISION_GEOMETRY<VECTOR<double,1> > const&,GRID<VECTOR<double,1> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<double,1>::m> >&,bool,double const&,double const&);
template void Compute_Occupied_Blocks<double,VECTOR<double,2> >(COLLISION_GEOMETRY<VECTOR<double,2> > const&,GRID<VECTOR<double,2> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<double,2>::m> >&,bool,double const&,double const&);
template void Compute_Occupied_Blocks<double,VECTOR<double,3> >(COLLISION_GEOMETRY<VECTOR<double,3> > const&,GRID<VECTOR<double,3> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<double,3>::m> >&,bool,double const&,double const&);
template void Compute_Occupied_Blocks<float,VECTOR<float,1> >(COLLISION_GEOMETRY<VECTOR<float,1> > const&,GRID<VECTOR<float,1> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<float,1>::m> >&,bool,float const&,float const&);
template void Compute_Occupied_Blocks<float,VECTOR<float,2> >(COLLISION_GEOMETRY<VECTOR<float,2> > const&,GRID<VECTOR<float,2> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<float,2>::m> >&,bool,float const&,float const&);
template void Compute_Occupied_Blocks<float,VECTOR<float,3> >(COLLISION_GEOMETRY<VECTOR<float,3> > const&,GRID<VECTOR<float,3> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<float,3>::m> >&,bool,float const&,float const&);
template void Rasterize_Object<VECTOR<double,1> >(COLLISION_GEOMETRY<VECTOR<double,1> > const&,GRID<VECTOR<double,1> > const&,
    OBJECTS_IN_CELL<VECTOR<double,1>,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<double,2> >(COLLISION_GEOMETRY<VECTOR<double,2> > const&,GRID<VECTOR<double,2> > const&,
    OBJECTS_IN_CELL<VECTOR<double,2>,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<double,3> >(COLLISION_GEOMETRY<VECTOR<double,3> > const&,GRID<VECTOR<double,3> > const&,
    OBJECTS_IN_CELL<VECTOR<double,3>,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<float,1> >(COLLISION_GEOMETRY<VECTOR<float,1> > const&,GRID<VECTOR<float,1> > const&,
    OBJECTS_IN_CELL<VECTOR<float,1>,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<float,2> >(COLLISION_GEOMETRY<VECTOR<float,2> > const&,GRID<VECTOR<float,2> > const&,
    OBJECTS_IN_CELL<VECTOR<float,2>,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<float,3> >(COLLISION_GEOMETRY<VECTOR<float,3> > const&,GRID<VECTOR<float,3> > const&,
    OBJECTS_IN_CELL<VECTOR<float,3>,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
}
}
