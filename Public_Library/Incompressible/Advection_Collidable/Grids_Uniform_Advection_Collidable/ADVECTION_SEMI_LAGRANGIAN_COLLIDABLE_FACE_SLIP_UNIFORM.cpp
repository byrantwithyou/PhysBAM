//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM  
//##################################################################### 
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Advection/ADVECTION.h>
#include <Grid_PDE/Interpolation/AVERAGING_UNIFORM.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T_FACE_LOOKUP> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<TV,T_FACE_LOOKUP>::
ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list_input)
    :body_list(body_list_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T_FACE_LOOKUP> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<TV,T_FACE_LOOKUP>::
~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM()
{
}
//#####################################################################
// Function Update_Advection_Equation_Face_Lookup
//#####################################################################
template<class TV,class T_FACE_LOOKUP> void ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<TV,T_FACE_LOOKUP>::
Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& Z,const T_FACE_LOOKUP& Z_ghost,const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,ARRAY<T,FACE_INDEX<TV::m> >* Z_min,ARRAY<T,FACE_INDEX<TV::m> >* Z_max)
{
    PHYSBAM_ASSERT(!Z_min_ghost && !Z_max_ghost && !Z_min && !Z_max);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
        if(!body_list.Occupied_Face_Center(iterator.Full_Index())){
            TV grid_point_location=iterator.Location(),length_and_direction=-dt*averaging.Face_To_Face_Vector(grid,axis,face,face_velocities.Nested()),
                interpolation_point=grid_point_location+length_and_direction;
            Z(axis,face)=linear_interpolation.Clamped_To_Array_Face_Component(axis,grid,Z_ghost.Nested().Starting_Point_Face(axis,face),interpolation_point);}
        else{
            const typename T_FACE_LOOKUP::LOOKUP& lookup=Z_ghost.Starting_Point_Face(axis,face);
            TV velocity=AVERAGING_UNIFORM<TV,T_FACE_LOOKUP>::Average_Face_To_Face_Vector_Helper(grid,iterator.Full_Index(),lookup);
            TV length_and_direction=-dt*velocity;
            TV_INT adjacent_cell_center=iterator.First_Cell_Index();
            if((*body_list.outside_fluid)(iterator.First_Cell_Index()))
                adjacent_cell_center=iterator.Second_Cell_Index();

            TV interpolation_point=iterator.Location()+length_and_direction;
            TV cell_center_location=grid.Center(adjacent_cell_center);
            length_and_direction=interpolation_point-cell_center_location;
            RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
            if(RAY<TV>::Create_Non_Degenerate_Ray(cell_center_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id)){
                int aggregate_id=0;
                body_list.collision_geometry_collection.Intersection_Between_Points(cell_center_location,cell_center_location+length_and_direction,body_id,aggregate_id,interpolation_point);
                Z(axis,face)=body_list.Object_Velocity(body_id,aggregate_id,interpolation_point)[axis];}
            else if(body_list.Inside_Any_Body(cell_center_location,body_id)){
                COLLISION_GEOMETRY<TV>& collision_geometry=body_list.collision_geometry_collection(body_id);
                int simplex_id;
                if(!collision_geometry.Has_Volumetric_Geometry())
                    Z(axis,face)=collision_geometry.Pointwise_Object_Velocity(cell_center_location)[axis];
                else if(collision_geometry.Inside_Any_Simplex(cell_center_location,simplex_id))
                    Z(axis,face)=collision_geometry.Pointwise_Object_Velocity(simplex_id,cell_center_location)[axis];
                else
                    PHYSBAM_FATAL_ERROR("Inconsistent inside checks");}
            else{
                const typename T_FACE_LOOKUP::LOOKUP lookup(Z_ghost,Z_ghost.nested_face_lookup);
                Z(axis,face)=linear_interpolation_collidable.From_Block_Face_Component(axis,grid,BLOCK_UNIFORM<TV>(grid,interpolation_point,lookup.Number_Of_Ghost_Cells()),lookup,interpolation_point);}}}
}
//#####################################################################
// Function Average_To_Invalidated_Face
//#####################################################################
template<class TV,class T_FACE_LOOKUP> void ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<TV,T_FACE_LOOKUP>::
Average_To_Invalidated_Face(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_values,ARRAY<bool,FACE_INDEX<TV::m> >* faces_not_to_revalidate)
{
    return;
    PHYSBAM_FATAL_ERROR("Not doing this");
}
namespace PhysBAM{
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<VECTOR<float,1>,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<float,1>,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<VECTOR<float,2>,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<VECTOR<float,3>,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<VECTOR<double,1>,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<double,1>,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<VECTOR<double,2>,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<VECTOR<double,3>,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > > >;
}
