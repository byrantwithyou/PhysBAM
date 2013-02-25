//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace RASTERIZATION
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/BLOCK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Rasterization/RIGID_GEOMETRY_RASTERIZATION_HELPER.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class TV> class GRID;

namespace RASTERIZATION{
//#####################################################################
// Function Rasterize_Object
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Object(const COLLISION_GEOMETRY<TV>& collision_geometry,const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID& id)
{
    Rasterize_Object_Generic(collision_geometry,grid,objects,id);
}
//#####################################################################
// Function Compute_Occupied_Blocks
//#####################################################################
template<class T,class TV,class T_GRID> void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<TV>& collision_geometry,const T_GRID& grid,
    ARRAY<bool,VECTOR<int,TV::m> >& occupied,const bool with_body_motion,const T& extra_thickness,const T& body_thickness_factor)
{
    typedef typename REBIND<TV,int>::TYPE TV_INT;
    if(collision_geometry.Number_Of_Simplices()) Compute_Occupied_Blocks_Generic(collision_geometry,grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
    else{
        for(UNIFORM_GRID_ITERATOR_NODE<TV> node_iterator(grid,2);node_iterator.Valid();node_iterator.Next()){
            TV_INT block_index=node_iterator.Node_Index();BLOCK_UNIFORM<T_GRID> block(grid,block_index);
            for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++){TV_INT cell=block.Cell(cell_index);
                if(collision_geometry.Implicit_Geometry_Extended_Value(grid.X(cell))<=extra_thickness){occupied(block_index)=true;break;}}}}
}
//#####################################################################
// Function Rasterize_Box_Onto_Blocks
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Box_Onto_Blocks(const T_GRID& grid,ARRAY<bool,VECTOR<int,TV::m> >& occupied,const RANGE<TV>& box)
{
    TV DX_over_two=(typename TV::SCALAR).5*grid.dX;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,grid.Clamp_To_Cell(box.Translated(-DX_over_two),3));iterator.Valid();iterator.Next())
        occupied(iterator.Cell_Index()+1)=true;
}
//#####################################################################
// Function Rasterize_Box
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Box(const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const RANGE<TV>& box,const COLLISION_GEOMETRY_ID id)
{
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,grid.Clamp_To_Cell(box,3));iterator.Valid();iterator.Next())
        objects_in_cell.Add_Object_To_Cell(iterator.Cell_Index(),id);
}
//#####################################################################
template void Compute_Occupied_Blocks<float,VECTOR<float,1>,GRID<VECTOR<float,1> > >(COLLISION_GEOMETRY<VECTOR<float,1> > const&,GRID<VECTOR<float,1> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<float,1>::m> >&,bool,float const&,float const&);
template void Compute_Occupied_Blocks<float,VECTOR<float,2>,GRID<VECTOR<float,2> > >(COLLISION_GEOMETRY<VECTOR<float,2> > const&,GRID<VECTOR<float,2> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<float,2>::m> >&,bool,float const&,float const&);
template void Compute_Occupied_Blocks<float,VECTOR<float,3>,GRID<VECTOR<float,3> > >(COLLISION_GEOMETRY<VECTOR<float,3> > const&,GRID<VECTOR<float,3> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<float,3>::m> >&,bool,float const&,float const&);
template void Rasterize_Object<VECTOR<double,1>,GRID<VECTOR<double,1> > >(COLLISION_GEOMETRY<VECTOR<double,1> > const&,GRID<VECTOR<double,1> > const&,
    OBJECTS_IN_CELL<GRID<VECTOR<double,1> >,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<double,2>,GRID<VECTOR<double,2> > >(COLLISION_GEOMETRY<VECTOR<double,2> > const&,GRID<VECTOR<double,2> > const&,
    OBJECTS_IN_CELL<GRID<VECTOR<double,2> >,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<double,3>,GRID<VECTOR<double,3> > >(COLLISION_GEOMETRY<VECTOR<double,3> > const&,GRID<VECTOR<double,3> > const&,
    OBJECTS_IN_CELL<GRID<VECTOR<double,3> >,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Compute_Occupied_Blocks<double,VECTOR<double,1>,GRID<VECTOR<double,1> > >(COLLISION_GEOMETRY<VECTOR<double,1> > const&,GRID<VECTOR<double,1> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<double,1>::m> >&,bool,double const&,double const&);
template void Compute_Occupied_Blocks<double,VECTOR<double,2>,GRID<VECTOR<double,2> > >(COLLISION_GEOMETRY<VECTOR<double,2> > const&,GRID<VECTOR<double,2> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<double,2>::m> >&,bool,double const&,double const&);
template void Compute_Occupied_Blocks<double,VECTOR<double,3>,GRID<VECTOR<double,3> > >(COLLISION_GEOMETRY<VECTOR<double,3> > const&,GRID<VECTOR<double,3> > const&,
    ARRAY<bool,VECTOR<int,VECTOR<double,3>::m> >&,bool,double const&,double const&);
template void Rasterize_Object<VECTOR<float,1>,GRID<VECTOR<float,1> > >(COLLISION_GEOMETRY<VECTOR<float,1> > const&,GRID<VECTOR<float,1> > const&,
    OBJECTS_IN_CELL<GRID<VECTOR<float,1> >,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<float,2>,GRID<VECTOR<float,2> > >(COLLISION_GEOMETRY<VECTOR<float,2> > const&,GRID<VECTOR<float,2> > const&,
    OBJECTS_IN_CELL<GRID<VECTOR<float,2> >,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
template void Rasterize_Object<VECTOR<float,3>,GRID<VECTOR<float,3> > >(COLLISION_GEOMETRY<VECTOR<float,3> > const&,GRID<VECTOR<float,3> > const&,
    OBJECTS_IN_CELL<GRID<VECTOR<float,3> >,COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID const&);
};
};
