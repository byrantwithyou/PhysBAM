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
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Incompressible/Collisions_And_Interactions/OBJECTS_IN_CELL.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV,class COLLISION_GEOMETRY_ID> class OBJECTS_IN_CELL;

namespace RASTERIZATION{
template<class TV> void Rasterize_Object_Generic(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID id);
template<class T,class TV> void
Compute_Occupied_Blocks_Generic(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,ARRAY<bool,VECTOR<int,TV::m> >& occupied,
    const bool with_body_motion,const T extra_thickness,const T body_thickness_factor);
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
// Function Prepare_For_World_Space_Simplex_Bounding_Box
//#####################################################################
template<class T> void Prepare_For_World_Space_Simplex_Bounding_Box(const RIGID_BODY<VECTOR<T,1> >& rigid_body)
{}
//#####################################################################
// Function Prepare_For_World_Space_Simplex_Bounding_Box
//#####################################################################
template<class T> void Prepare_For_World_Space_Simplex_Bounding_Box(const RIGID_BODY<VECTOR<T,2> >& rigid_body)
{
    if(!rigid_body.simplicial_object->segment_list) const_cast<RIGID_BODY<VECTOR<T,2> >&>(rigid_body).simplicial_object->Update_Segment_List();
}
//#####################################################################
// Function Prepare_For_World_Space_Simplex_Bounding_Box
//#####################################################################
template<class T> void Prepare_For_World_Space_Simplex_Bounding_Box(const RIGID_BODY<VECTOR<T,3> >& rigid_body)
{
    if(!rigid_body.simplicial_object->triangle_list) const_cast<RIGID_BODY<VECTOR<T,3> >&>(rigid_body).simplicial_object->Update_Triangle_List();
}
//#####################################################################
// Function Rasterize_Object_Generic
//#####################################################################
template<class TV> void Rasterize_Object_Generic(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID id)
{
    const RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry;
    rigid_collision_geometry=dynamic_cast<const RIGID_COLLISION_GEOMETRY<TV>*>(&collision_geometry);
    if(rigid_collision_geometry) Prepare_For_World_Space_Simplex_Bounding_Box(rigid_collision_geometry->rigid_body);
    if(collision_geometry.Number_Of_Simplices())
        for(int t=0;t<collision_geometry.Number_Of_Simplices();t++){
            typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SIMPLEX;
            T_SIMPLEX simplex=collision_geometry.World_Space_Simplex(t);
            VECTOR<TV,TV::dimension> pts;
            for(int d=0;d<TV::dimension;d++) pts(d)=simplex.X(d);
            RANGE<TV> box=RANGE<TV>::Bounding_Box(pts).Thickened(collision_geometry.collision_thickness);
            Rasterize_Box(grid,objects,box,id);}
    else
        for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry.collision_geometries_for_rasterization->m;i++)
            if((*collision_geometry.collision_geometries_for_rasterization)(i) && (*collision_geometry.collision_geometries_for_rasterization)(i)->active)
                Rasterize_Object(*(*collision_geometry.collision_geometries_for_rasterization)(i),grid,objects,id);
}
//#####################################################################
// Macro SPECIALIZE_RASTERIZE_OBJECT
//#####################################################################
#define SPECIALIZE_RASTERIZE_OBJECT(TV) \
    void Rasterize_Object_Generic(const COLLISION_GEOMETRY<TV>&,const GRID<TV>&,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>&,const COLLISION_GEOMETRY_ID&) \
    {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Compute_Occupied_Blocks_Generic
//#####################################################################
// NOTE: Since Compute_Occupied_Blocks is a virtual function we can't directly make it a template member
template<class T,class TV> void
Compute_Occupied_Blocks_Generic(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,ARRAY<bool,VECTOR<int,TV::m> >& occupied,
    const bool with_body_motion,const T extra_thickness,const T body_thickness_factor)
{
    const RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry;
    rigid_collision_geometry=dynamic_cast<const RIGID_COLLISION_GEOMETRY<TV>*>(&collision_geometry);
    if(rigid_collision_geometry) Prepare_For_World_Space_Simplex_Bounding_Box(rigid_collision_geometry->rigid_body);
    for(int t=0;t<collision_geometry.Number_Of_Simplices();t++){
        typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SIMPLEX;
        T_SIMPLEX simplex=collision_geometry.World_Space_Simplex(t);
        VECTOR<TV,TV::dimension> pts;
        for(int d=0;d<TV::dimension;d++) pts(d)=simplex.X(d);
        RANGE<TV> box=RANGE<TV>::Bounding_Box(pts);
        if(with_body_motion){
            T_SIMPLEX simplex_saved=collision_geometry.World_Space_Simplex(t,1);
            VECTOR<TV,TV::dimension> pts_saved;
            for(int d=0;d<TV::dimension;d++) pts_saved(d)=simplex_saved.X(d);
            box.Enlarge_To_Include_Box(RANGE<TV>::Bounding_Box(pts_saved));}
        box.Change_Size(extra_thickness+body_thickness_factor*collision_geometry.collision_thickness);
        Rasterize_Box_Onto_Blocks(grid,occupied,box);}
}
//#####################################################################
// Macro SPECIALIZE_COMPUTE_OCCUPIED
//#####################################################################
#define SPECIALIZE_COMPUTE_OCCUPIED(TV) \
    void Compute_Occupied_Blocks_Generic(const COLLISION_GEOMETRY<TV>&,const GRID<TV>&,ARRAY<bool>&,const bool,const TV::SCALAR,const TV::SCALAR) \
    {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
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
