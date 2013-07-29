//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM  
//##################################################################### 
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM__

#include <Tools/Advection/ADVECTION.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <Incompressible/Interpolation_Collidable/AVERAGING_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<TV>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM:public ADVECTION<TV,typename TV::SCALAR,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
public:
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM<TV,T,T_FACE_LOOKUP> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_UNIFORM<TV,T,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_UNIFORM<TV,typename T_FACE_LOOKUP::NESTED_LOOKUP> averaging;
    AVERAGING_COLLIDABLE_UNIFORM<TV,T_FACE_LOOKUP> averaging_collidable;
    ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid_mask;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list_input,ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid_mask_input)
        :body_list(body_list_input),averaging_collidable(body_list,0),
         face_velocities_valid_mask(face_velocities_valid_mask_input)
    {}

    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM()
    {}

    void Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,ARRAY<T,FACE_INDEX<TV::m> >* Z_min,ARRAY<T,FACE_INDEX<TV::m> >* Z_max)
    {ARRAY<bool,FACE_INDEX<TV::m> > face_velocities_valid_mask_next(grid,3,false);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
        if(!body_list.Swept_Occupied_Face_Center(iterator)){
            TV interpolation_point=iterator.Location()-dt*averaging.Face_To_Face_Vector(grid,axis,face,face_velocities.Nested());
            Z.Component(axis)(face)=linear_interpolation.Clamped_To_Array_Face_Component(axis,grid,Z_ghost.Nested().Starting_Point_Face(axis,face),interpolation_point);
            face_velocities_valid_mask_next.Component(axis)(face)=true;
            if(Z_min && Z_max){
                VECTOR<T,2> extrema=linear_interpolation.Extrema_Clamped_To_Array_Face_Component(axis,grid,Z_min_ghost->Nested().Starting_Point_Face(axis,face),
                    Z_max_ghost->Nested().Starting_Point_Face(axis,face),interpolation_point);
                Z_min->Component(axis)(face)=extrema.x;Z_max->Component(axis)(face)=extrema.y;}}
        else{
            T face_velocity=0; // avoid uninitialization warning
            TV grid_point_location=iterator.Location(),length_and_direction=-dt*averaging_collidable.Face_To_Face_Vector(grid,axis,face,face_velocities),
                interpolation_point=grid_point_location+length_and_direction;
            if(body_list.Latest_Velocity_Crossover(axis,face,dt,face_velocity)){
                face_velocities_valid_mask_next.Component(axis)(face)=false;Z.Component(axis)(face)=face_velocity;
                if(Z_min && Z_max){Z_min->Component(axis)(face)=Z_max->Component(axis)(face)=face_velocity;}}
            else{
                RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                if(RAY<TV>::Create_Non_Degenerate_Ray(grid_point_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id)) 
                    interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                const typename T_FACE_LOOKUP::LOOKUP& lookup=Z_ghost.Starting_Point_Face(axis,face);
                Z.Component(axis)(face)=linear_interpolation_collidable.Clamped_To_Array_Face_Component(axis,grid,lookup,interpolation_point);
                face_velocities_valid_mask_next.Component(axis)(face)=lookup.found_valid_point;
                if(Z_min && Z_max){
                    const typename T_FACE_LOOKUP::LOOKUP &Z_min_lookup=Z_min_ghost->Starting_Point_Face(axis,face),&Z_max_lookup=Z_max_ghost->Starting_Point_Face(axis,face);
                    VECTOR<T,2> extrema=linear_interpolation_collidable.Extrema_Clamped_To_Array_Face_Component(axis,grid,Z_min_lookup,Z_max_lookup,interpolation_point);
                    Z_min->Component(axis)(face)=extrema.x;Z_max->Component(axis)(face)=extrema.y;}}}}
    ARRAY<bool,FACE_INDEX<TV::m> >::Exchange(face_velocities_valid_mask,face_velocities_valid_mask_next);
    // ghost values should always be valid
    for(FACE_ITERATOR<TV> iterator(grid,3,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next())
        face_velocities_valid_mask(iterator.Full_Index())=true;}

    T Compute_Revalidation_Value(const int axis,const TV& from,const TV& to,const T& current_invalid_value,const T& default_value)
    {TV point;COLLISION_GEOMETRY_ID body_id;int aggregate_id;
    if(body_list.collision_geometry_collection.Intersection_Between_Points(from,to,body_id,aggregate_id,point)) return body_list.Object_Velocity(body_id,aggregate_id,point)[axis];
    else return default_value;}

    void Average_To_Invalidated_Face(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_values)
    {// average values collision aware in Gauss-Jacobi fashion
    VECTOR<ARRAY<PAIR<TV_INT,bool> >,TV::m> face_invalid_indices; // index and bool true if entry has been validated on iteration
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(!face_velocities_valid_mask.Component(iterator.Axis())(iterator.Face_Index())) 
        face_invalid_indices[iterator.Axis()].Append(PAIR<TV_INT,bool>(iterator.Face_Index(),false));
    
    for(int arrays_axis=0;arrays_axis<TV::m;arrays_axis++){
        ARRAY<PAIR<TV_INT,bool> >& invalid_indices=face_invalid_indices[arrays_axis];
        ARRAYS_ND_BASE<VECTOR<bool,TV::m>,TV_INT>& neighbors_visible=body_list.face_neighbors_visible.Component(arrays_axis);
        ARRAYS_ND_BASE<bool,TV_INT>& valid_points=face_velocities_valid_mask.Component(arrays_axis);T_ARRAYS_BASE& values=face_values.Component(arrays_axis);
        
        bool done=false;
        for(CELL_ITERATOR<TV> iterator(grid,3,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next())
            valid_points(iterator.index)=false; // don't average from boundaries

        while(!done){
            done=true;
            for(int k=0;k<invalid_indices.m;k++){ 
                T sum=0;int count=0;
                for(int axis=0;axis<TV::m;axis++){
                    TV_INT min_face=invalid_indices(k).x-TV_INT::Axis_Vector(axis),max_face=invalid_indices(k).x+TV_INT::Axis_Vector(axis);
                    if(neighbors_visible(min_face)(axis) && valid_points(min_face)){sum+=values(min_face);count++;}
                    if(neighbors_visible(invalid_indices(k).x)(axis) && valid_points(max_face)){sum+=values(max_face);count++;}}
                if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}}
            if(!done) for(int k=invalid_indices.m-1;k>=0;k--) if(invalid_indices(k).y){valid_points(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}

        // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
        done=false;
        while(!done){
            done=true;
            for(int k=0;k<invalid_indices.m;k++){ 
                T sum=0;int count=0;
                for(int axis=0;axis<TV::m;axis++){
                    TV_INT min_face=invalid_indices(k).x-TV_INT::Axis_Vector(axis),max_face=invalid_indices(k).x+TV_INT::Axis_Vector(axis);
                    if(neighbors_visible(min_face)(axis)){if(valid_points(min_face)){sum+=values(min_face);count++;}}
                    else{sum+=Compute_Revalidation_Value(arrays_axis,grid.X(invalid_indices(k).x),grid.X(min_face),values(invalid_indices(k).x),T());count++;}
                    if(neighbors_visible(invalid_indices(k).x)(axis)){if(valid_points(max_face)){sum+=values(max_face);count++;}}
                    else{sum+=Compute_Revalidation_Value(arrays_axis,grid.X(invalid_indices(k).x),grid.X(max_face),values(invalid_indices(k).x),T());count++;}}
                if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}
                else values(invalid_indices(k).x)=T();}
            if(!done) for(int k=invalid_indices.m-1;k>=0;k--) if(invalid_indices(k).y){valid_points(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
        for(CELL_ITERATOR<TV> iterator(grid,3,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next())
            valid_points(iterator.index)=true;}} // set valid for future advection

//#####################################################################
    };
}
#endif
