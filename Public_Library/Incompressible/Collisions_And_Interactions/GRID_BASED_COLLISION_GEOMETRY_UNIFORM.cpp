//#####################################################################
// Copyright 2005-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/PAIR.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Intersections/BOX_POINT_SIMPLEX_1D_INTERSECTION.h>
#include <Geometry/Intersections/BOX_SEGMENT_2D_INTERSECTION.h>
#include <Geometry/Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Intersections/RAY_POINT_SIMPLEX_1D_INTERSECTION.h>
#include <Geometry/Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
#include <Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Collisions_And_Interactions/RIGID_BODY_RASTERIZATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
GRID_BASED_COLLISION_GEOMETRY_UNIFORM(GRID<TV>& grid_input)
    :GRID_BASED_COLLISION_GEOMETRY<TV>(grid_input),use_collision_face_neighbors(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
~GRID_BASED_COLLISION_GEOMETRY_UNIFORM()
{}
//##################################################################### 
// Function Initialize_Grids
//##################################################################### 
template<class TV> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Initialize_Grids()
{
    collision_thickness=(T)1e-3*grid.dX.Min();
    collision_geometry_collection.Set_Collision_Body_Thickness(collision_thickness);

    assert(grid.Is_MAC_Grid());
    VECTOR<bool,TV::m> all_true;all_true.Fill(true);
    cell_neighbors_visible.Resize(grid.Domain_Indices(3),false);cell_neighbors_visible.Fill(all_true); // initialize here so collision aware redistancing works in Initialize
    face_neighbors_visible.Resize(grid,1,false);face_neighbors_visible.Fill(all_true);
}
//##################################################################### 
// Function Compute_Occupied_Blocks
//##################################################################### 
template<class TV> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Compute_Occupied_Blocks(const bool with_body_motion,const T extra_thickness,const T body_thickness_factor)
{
    ARRAY<bool,TV_INT>& occupied=with_body_motion?swept_occupied_blocks:occupied_blocks;
    occupied.Resize(grid.Block_Indices(3),false,false);occupied.Fill(false);
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++)
        if(collision_geometry_collection.bodies(i) && collision_geometry_collection.bodies(i)->active)
            RASTERIZATION::Compute_Occupied_Blocks(*collision_geometry_collection.bodies(i),grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
}
//#####################################################################
// Function Compute_Grid_Visibility
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Compute_Grid_Visibility()
{
    VECTOR<bool,TV::m> all_true;all_true.Fill(true);
    cell_neighbors_visible.Fill(all_true);face_neighbors_visible.Fill(all_true);
    if(!collision_geometry_collection.bodies.m) return;

    // cell neighbors
    for(FACE_ITERATOR<TV> iterator(grid,3,GRID<TV>::INTERIOR_REGION);iterator.Valid();iterator.Next()){
        TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1,cell2,collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
        if(collision_geometry_collection.Intersection_Between_Points(grid.Center(cell1),grid.Center(cell2),&objects)) cell_neighbors_visible(cell1)(iterator.Axis())=false;}

    // face neighbors
    for(FACE_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
        for(int direction=0;direction<TV::m;direction++){TV_INT direction_offset=TV_INT::Axis_Vector(direction);
            ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(iterator.Second_Cell_Index(),iterator.Second_Cell_Index()+direction_offset,collision_geometry_collection.bodies.m,objects);
            if(!objects.m) continue;
            if(collision_geometry_collection.Intersection_Between_Points(iterator.Location(),grid.Face(FACE_INDEX<TV::m>(axis,iterator.Face_Index()+direction_offset)),&objects))
                face_neighbors_visible.Component(axis)(iterator.Face_Index())(direction)=false;}}
}
//#####################################################################
// Function Compute_Psi_N
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Compute_Psi_N(ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,FACE_INDEX<TV::m> >* face_velocities) const
{
    if(!collision_geometry_collection.bodies.m) return;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1,cell2,collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
        int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();COLLISION_GEOMETRY_ID body_id;int count=0;T velocity=0;
        RAY<TV> ray(iterator.First_Cell_Center(),TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)){count++;
            velocity=collision_geometry_collection(body_id).Pointwise_Object_Pseudo_Velocity(ray.aggregate_id,ray.Point(ray.t_max),COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE)[axis];}
        ray.Initialize(iterator.Second_Cell_Center(),-TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)){count++;
            velocity+=collision_geometry_collection(body_id).Pointwise_Object_Pseudo_Velocity(ray.aggregate_id,ray.Point(ray.t_max),COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE)[axis];}
        if(count){psi_N.Component(axis)(face_index)=true;if(face_velocities) (*face_velocities).Component(axis)(face_index)=velocity/count;}}
}
//#####################################################################
// Function Compute_Psi_N
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Compute_Psi_N_Zero_Velocity(ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,FACE_INDEX<TV::m> >* face_velocities) const
{
    if(!collision_geometry_collection.bodies.m) return;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1,cell2,collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
        int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();COLLISION_GEOMETRY_ID body_id;int count=0;
        RAY<TV> ray(iterator.First_Cell_Center(),TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)) count++;
        ray.Initialize(iterator.Second_Cell_Center(),-TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)) count++;
        if(count){psi_N.Component(axis)(face_index)=true;if(face_velocities) (*face_velocities)(axis,face_index)=(T)0;}}
}
//#####################################################################
// Function Compute_Simplices_In_Cell
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV_INT>& simplices_in_cell,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies,
    int ghost_cells,T thickness,bool assume_active) const
{
    simplices_in_cell.Resize(grid.Domain_Indices(ghost_cells+1));
    for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++) if(assume_active || Is_Active(i)){ // Is_Active is not safe to call if different bodies list is used
        COLLISION_GEOMETRY<TV>* body=bodies(i);
        int n=body->Number_Of_Simplices();
        for(int e=0;e<n;e++){
            RANGE<TV> box=body->World_Space_Simplex(e).Bounding_Box().Thickened(thickness);
            for(CELL_ITERATOR<TV> iterator(grid,grid.Clamp_To_Cell(box,ghost_cells+1));iterator.Valid();iterator.Next())
                simplices_in_cell(iterator.Cell_Index()).Append(PAIR<COLLISION_GEOMETRY_ID,int>(i,e));}}
}
//#####################################################################
// Function Inside_Any_Simplex_Of_Any_Body
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Inside_Any_Simplex_Of_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id) const
{
    ARRAY<COLLISION_GEOMETRY_ID> objects;
    TV_INT cell_index=grid.Clamp_To_Cell(location);
    objects_in_cell.Get_Objects_For_Cell(cell_index,objects);
    if(!objects.m) return false;
    return collision_geometry_collection.Inside_Any_Simplex_Of_Any_Body(location,body_id,aggregate_id,&objects);
}
//#####################################################################
// Function Inside_Any_Body
//#####################################################################
// TODO(jontg): Slow, but it works...
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Inside_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id) const
{
    return collision_geometry_collection.Inside_Any_Body(location,collision_thickness*(T).5,body_id);
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside_Any_Body
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Implicit_Geometry_Lazy_Inside_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id) const
{
    return collision_geometry_collection.Implicit_Geometry_Lazy_Inside_Any_Body(location,body_id);
}
//#####################################################################
// Function Closest_Non_Intersecting_Point_Of_Any_Body
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Closest_Non_Intersecting_Point_Of_Any_Body(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id) const
{
    // since the ray is a general ray, we don't get specific id's.
    return collision_geometry_collection.Closest_Non_Intersecting_Point_Of_Any_Simplicial_Object(ray,body_id);
}
//#####################################################################
// Function Cell_Center_Visible_From_Face
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Cell_Center_Visible_From_Face(const TV_INT& cell,const int axis,const TV_INT& face_index) const
{
    ARRAY<COLLISION_GEOMETRY_ID> objects;
    objects_in_cell.Get_Objects_For_Cell(cell,objects);
    if(!objects.m) return true;
    return !collision_geometry_collection.Intersection_Between_Points(grid.Center(cell),grid.Face(FACE_INDEX<TV::m>(axis,face_index)),&objects);
}
//#####################################################################
// Function Face_Velocity
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Face_Velocity(const int axis,const TV_INT& face_index,const TV_INT* cells,const int number_of_cells,const TV& X,T& face_velocity) const
{
    ARRAY<COLLISION_GEOMETRY_ID> objects;
    objects_in_cell.Get_Objects_For_Cells(cells,number_of_cells,collision_geometry_collection.bodies.m,objects);
    if(!objects.m) return false;
    COLLISION_GEOMETRY_ID body_id;
    int triangle_id;
    TV intersection_point;
    if(collision_geometry_collection.Intersection_Between_Points(X,grid.Face(FACE_INDEX<TV::m>(axis,face_index)),body_id,triangle_id,intersection_point,&objects)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(triangle_id,intersection_point)[axis];
        return true;}
    return false;
}
//#####################################################################
// Function Face_Velocity
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Face_Velocity(const int side,const int axis,const TV_INT& face_index,const TV_INT* cells,const int number_of_cells,const TV& X,T& face_velocity) const
{
    ARRAY<COLLISION_GEOMETRY_ID> objects;
    objects_in_cell.Get_Objects_For_Cells(cells,number_of_cells,collision_geometry_collection.bodies.m,objects);
    if(!objects.m) return false;
    COLLISION_GEOMETRY_ID body_id;
    int triangle_id;
    TV intersection_point;
    if(collision_geometry_collection.Intersection_Between_Points(X,grid.X(side==0?face_index-TV_INT::Axis_Vector(axis):face_index),body_id,triangle_id,intersection_point,&objects)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(triangle_id,intersection_point)[axis];
        return true;}
    return false;
}
//#####################################################################
// Function Latest_Velocity_Crossover
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Latest_Velocity_Crossover(const int side,const int axis,const TV_INT& face_index,const T dt,T& face_velocity) const
{
    COLLISION_GEOMETRY_ID body_id;
    int aggregate_id;
    TV initial_hit_point,X=grid.Face(FACE_INDEX<TV::m>(axis,face_index));
    if(Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point)[axis];
        return true;}
    ARRAY<COLLISION_GEOMETRY_ID> objects;
    objects_in_cell.Get_Objects_For_Cells(face_index,face_index-TV_INT::Axis_Vector(axis),collision_geometry_collection.bodies.m,objects);
    if(!objects.m) return false;
    int triangle_id;
    if(collision_geometry_collection.Intersection_Between_Points(grid.Center((side==0)?face_index:(face_index-TV_INT::Axis_Vector(axis))),X,body_id,triangle_id,initial_hit_point,&objects)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(triangle_id,initial_hit_point)[axis];
        return true;}
    return false;
}
//#####################################################################
// Function Cell_Center_Intersection
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Cell_Center_Intersection(const TV_INT& cell_index,const TV_INT* cell_indices_for_body_id,const int number_of_cells_for_body_id,const TV& X,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id,
    TV& intersection_point) const
{
    ARRAY<COLLISION_GEOMETRY_ID> objects;
    objects_in_cell.Get_Objects_For_Cells(cell_indices_for_body_id,number_of_cells_for_body_id,collision_geometry_collection.bodies.m,objects);
    if(!objects.m) return false;
    return collision_geometry_collection.Intersection_Between_Points(grid.Center(cell_index),X,body_id,aggregate_id,intersection_point,&objects);}
//#####################################################################
// Function Latest_Velocity_Crossover
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Latest_Velocity_Crossover(const int axis,const TV_INT& face_index,const T dt,T& face_velocity) const
{
    COLLISION_GEOMETRY_ID body_id;
    int aggregate_id;
    TV initial_hit_point,X=grid.Face(FACE_INDEX<TV::m>(axis,face_index));
    if(Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point)[axis];
        return true;}
    return false;
}
//#####################################################################
// Function Compute_Simplices_In_Cell
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV_INT>& simplices_in_cell,int ghost_cells,T thickness) const
{
    Compute_Simplices_In_Cell(simplices_in_cell,collision_geometry_collection.bodies,ghost_cells,thickness,false);
}
//#####################################################################
// Function Latest_Cell_Crossover_And_Velocity
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Latest_Cell_Crossover_And_Velocity(const TV_INT& cell_index,const T dt,TV& velocity) const
{
    COLLISION_GEOMETRY_ID body_id;
    int aggregate_id;
    TV initial_hit_point,X=grid.Center(cell_index);
    if(Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point)){
        velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point);
        return true;}
    return false;
}
//#####################################################################
// Function Object_Velocity
//#####################################################################
template<class TV> TV GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>::
Object_Velocity(const COLLISION_GEOMETRY_ID body_id,const int aggregate_id,const TV& X) const
{
    return collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,X);
}
//#####################################################################
namespace PhysBAM{
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<float,1> >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<float,2> >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<float,3> >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<double,1> >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<double,2> >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<double,3> >;
}
