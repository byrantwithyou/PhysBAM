//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Level_Sets/FAST_MARCHING.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Level_Sets/FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<T_GRID>::
FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM(const LEVELSET_COLLIDABLE<TV>& levelset_collidable,const int ghost_cells,THREAD_QUEUE* thread_queue_input)
    :FAST_MARCHING_METHOD_UNIFORM<T_GRID>(levelset_collidable,ghost_cells,thread_queue_input),levelset_collidable(levelset_collidable)
{
    Neighbor_Visible=[&](const int neighbor_number,const TV_INT& current_index) -> bool
        {return !levelset_collidable.collision_body_list->cell_neighbors_visible.Valid_Index(current_index) || levelset_collidable.collision_body_list->cell_neighbors_visible(current_index)(neighbor_number);};
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<T_GRID>::
~FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM()
{
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class T_GRID> void FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<T_GRID>::
Initialize_Interface_Threaded(RANGE<TV_INT>& domain,ARRAY<T,TV_INT>& phi_ghost,ARRAY<T,TV_INT>& phi_new,ARRAY<bool,TV_INT>& done)
{ 
    LEVELSET<TV> levelset_ghost(cell_grid,phi_ghost);

    T fmm_initialization_iterative_tolerance_absolute=levelset_collidable.fmm_initialization_iterative_tolerance*cell_grid.dX.Min();    
 
    COLLISION_GEOMETRY_ID body_id;
    for(CELL_ITERATOR<TV> iterator(cell_grid,domain);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Cell_Index();
        if(!done(index)) continue;
        T value[3]={0}; // the phi value to use in the given direction
        int number_of_axis=0; // the number of axis that we want to use later
        bool really_clamp_phi_with_collision_bodies=levelset_collidable.clamp_phi_with_collision_bodies&&phi_ghost(index)<=0;
        T abs_phi=abs(phi_ghost(index));
        TV location=iterator.Location();
        for(int axis=0;axis<T_GRID::dimension;axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;
            T dx=cell_grid.dX[axis];
            bool use_low=false,use_high=false;T value_low=0,value_high=0;
            if(index[axis]>dimension_start[axis]){
                if(!Neighbor_Visible(axis,low)){
                    RAY<TV> ray(location,-TV::Axis_Vector(axis),true);ray.t_max=dx;ray.semi_infinite=false;
                    levelset_collidable.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
                    T phi_object=max(ray.t_max,levelset_collidable.small_number);
                    use_low=true;if(really_clamp_phi_with_collision_bodies) value_low=min(abs_phi,phi_object);else value_low=phi_object;}
                else if(done(low)&&LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(low))){
                    use_low=true;value_low=LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(low))*dx;}}
            if(index[axis]<dimension_end[axis]){
                if(!Neighbor_Visible(axis,index)){
                    RAY<TV> ray(location,TV::Axis_Vector(axis),true);ray.t_max=dx;ray.semi_infinite=false;
                    levelset_collidable.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
                    T phi_object=max(ray.t_max,levelset_collidable.small_number);
                    use_high=true;if(really_clamp_phi_with_collision_bodies) value_high=min(abs_phi,phi_object);else value_high=phi_object;}
                else if(done(high)&&LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(high))){
                    use_high=true;value_high=LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(high))*dx;}}
            if(!use_low){if(use_high) value[number_of_axis]=value_high;else number_of_axis--;}
            else if(!use_high) value[number_of_axis]=value_low;
            else value[number_of_axis]=min(value_low,value_high);
            number_of_axis++;}
        if(number_of_axis==1) phi_new(index)=value[0];
        else if(number_of_axis==2){
            if(T d2=sqr(value[0])+sqr(value[1])) phi_new(index)=value[0]*value[1]/sqrt(d2);
            else phi_new(index)=0;}
        else{PHYSBAM_ASSERT(T_GRID::dimension==3); // 2d should never get to this point
            T value_xy=value[0]*value[1],value_xz=value[0]*value[2],value_yz=value[1]*value[2],d2=sqr(value_xy)+sqr(value_xz)+sqr(value_yz);
            if(d2) phi_new(index)=value_xy*value[2]/sqrt(d2);
            else phi_new(index)=min(value[0],value[1],value[2]);}
        if(levelset_collidable.refine_fmm_initialization_with_iterative_solver){
            TV vertex=location;T phi_value=phi_ghost(index);
            for(int iterations=0;iterations<levelset_collidable.fmm_initialization_iterations;iterations++){
                if(abs(phi_value)>10*cell_grid.dX.Min())break; // stop if it looks like it's blowing up
                vertex-=phi_value*levelset_ghost.Normal(vertex);phi_value=levelset_ghost.Phi(vertex);
                if(abs(phi_value)<=fmm_initialization_iterative_tolerance_absolute){
                    T new_phi_value=(vertex-location).Magnitude();
                    if((new_phi_value-phi_new(index))/max(fmm_initialization_iterative_tolerance_absolute,phi_new(index))<levelset_collidable.fmm_initialization_iterative_drift_fraction && new_phi_value>0)
                        phi_new(index)=new_phi_value;
                    break;}}}
        phi_new(index)*=LEVELSET_UTILITIES<T>::Sign(phi_ghost(index));}
}
//#####################################################################
namespace PhysBAM{
template class FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<GRID<VECTOR<float,1> > >;
template class FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> > >;
template class FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> > >;
template class FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<GRID<VECTOR<double,1> > >;
template class FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> > >;
template class FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> > >;
}
