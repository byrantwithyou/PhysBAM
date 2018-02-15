//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_BASED_COLLISION_GEOMETRY
//#####################################################################
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/BLOCK_UNIFORM.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY.h>
#include <Incompressible/Collisions_And_Interactions/OBJECTS_IN_CELL.h>
#include <Incompressible/Collisions_And_Interactions/RIGID_BODY_RASTERIZATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_BASED_COLLISION_GEOMETRY<TV>::
GRID_BASED_COLLISION_GEOMETRY(GRID<TV>& grid_input)
    :collision_geometry_collection(*new COLLISION_BODY_COLLECTION<TV>),grid(grid_input),number_of_ghost_cells(3),
    objects_in_cell(*new OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GRID_BASED_COLLISION_GEOMETRY<TV>::
~GRID_BASED_COLLISION_GEOMETRY()
{
    delete &collision_geometry_collection;
    delete &objects_in_cell;
}
//##################################################################### 
// Function Add_Bodies
//##################################################################### 
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Add_Bodies(RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
{
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i))
        collision_geometry_collection.Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body_collection.Rigid_Body(i)),i,true);
}
//##################################################################### 
// Function Add_Bodies
//##################################################################### 
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Add_Body(RIGID_BODY<TV>& rigid_body)
{
    collision_geometry_collection.Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
}
//##################################################################### 
// Function Rasterize_Objects
//##################################################################### 
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Rasterize_Objects()
{
    objects_in_cell.Reset(grid,number_of_ghost_cells);
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++)
        if(collision_geometry_collection.bodies(i) && collision_geometry_collection.Is_Active(i) && collision_geometry_collection.bodies(i)->active)
            RASTERIZATION::Rasterize_Object(*collision_geometry_collection.bodies(i),grid,objects_in_cell,i);
}
//##################################################################### 
// Function Earliest_Simplex_Crossover
//##################################################################### 
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,TV& weights,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    T min_time=FLT_MAX;bool collision=false;T current_hit_time;TV current_weights;int current_simplex_id;
    if(!objects){for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++) if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && 
        collision_geometry_collection.bodies(i)->Earliest_Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,current_simplex_id) && current_hit_time < min_time){
            min_time=hit_time=current_hit_time;weights=current_weights;body_id=i;simplex_id=current_simplex_id;collision=true;}}
    else for(int k=0;k<objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && 
        collision_geometry_collection.bodies(i)->Earliest_Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,current_simplex_id) && current_hit_time < min_time){
            min_time=hit_time=current_hit_time;weights=current_weights;body_id=i;simplex_id=current_simplex_id;collision=true;}}
    return collision;
}
//##################################################################### 
// Function Latest_Crossover
//##################################################################### 
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Latest_Crossover(const TV& start_X,const TV& end_X,const T dt,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,TV& initial_hit_point,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    T hit_time;TV weights;POINT_SIMPLEX_COLLISION_TYPE returned_collision_type;
    bool crossover=Latest_Simplex_Crossover(start_X,end_X,dt,hit_time,weights,body_id,simplex_id,returned_collision_type,objects);
    if(crossover){
        T_SIMPLEX initial_simplex=collision_geometry_collection(body_id).World_Space_Simplex(simplex_id);
        initial_hit_point=initial_simplex.Point_From_Barycentric_Coordinates(weights);}
    return crossover;
}
//##################################################################### 
// Function Latest_Simplex_Crossover
//##################################################################### 
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,TV& weights,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type,
    const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    returned_collision_type=POINT_SIMPLEX_NO_COLLISION;POINT_SIMPLEX_COLLISION_TYPE collision_type;
    T max_time=-FLT_MAX;bool collision=false;T current_hit_time;TV current_weights;int current_simplex_id;
    if(!objects){for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++) if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && 
        collision_geometry_collection.bodies(i)->Latest_Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,current_simplex_id,collision_type) && current_hit_time > max_time){
            max_time=hit_time=current_hit_time;weights=current_weights;body_id=i;simplex_id=current_simplex_id;collision=true;
            returned_collision_type=collision_type;}}
    else for(int k=0;k<objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && 
        collision_geometry_collection.bodies(i)->Latest_Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,current_simplex_id,collision_type) && current_hit_time > max_time){
            max_time=hit_time=current_hit_time;weights=current_weights;body_id=i;simplex_id=current_simplex_id;collision=true;
            returned_collision_type=collision_type;}}
    return collision;
}
//##################################################################### 
// Function Any_Simplex_Crossover
//##################################################################### 
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    if(!objects){for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++) if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && collision_geometry_collection.bodies(i)->Any_Simplex_Crossover(start_X,end_X,dt)) return true;}
    else for(int k=0;k<objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && collision_geometry_collection.bodies(i)->Any_Simplex_Crossover(start_X,end_X,dt)) return true;}
    return false;
}
//##################################################################### 
// Function Update_Intersection_Acceleration_Structures
//##################################################################### 
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1,const int state2)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Update_Intersection_Acceleration_Structures(use_swept_simplex_hierarchy,state1,state2);
}
//##################################################################### 
// Function Get_Body_Penetration
//##################################################################### 
// normal is flipped to ensure that start_phi is positive
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Get_Body_Penetration(const TV& start_X,const TV& end_X,const T contour_value,const T dt,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,T& start_phi,T& end_phi,TV& end_body_normal,TV& body_velocity,
    const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    T hit_time=FLT_MAX;T current_hit_time;int current_simplex_id;T current_start_phi,current_end_phi;TV current_end_body_normal,current_body_velocity;
    if(objects) for(int k=0;k<objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);
        if(collision_geometry_collection.bodies(i)->Get_Body_Penetration(start_X,end_X,contour_value,dt,current_hit_time,current_simplex_id,current_start_phi,current_end_phi,current_end_body_normal,
             current_body_velocity) && current_hit_time<hit_time){
             body_id=i;hit_time=current_hit_time;
             simplex_id=current_simplex_id;start_phi=current_start_phi;end_phi=current_end_phi;end_body_normal=current_end_body_normal;body_velocity=current_body_velocity;}}
    else for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.Size();i++) 
        if(Is_Active(i) && collision_geometry_collection.bodies(i)->Get_Body_Penetration(start_X,end_X,contour_value,dt,current_hit_time,current_simplex_id,current_start_phi,current_end_phi,current_end_body_normal,
                current_body_velocity) && current_hit_time<hit_time){
            body_id=i;hit_time=current_hit_time;
            simplex_id=current_simplex_id;start_phi=current_start_phi;end_phi=current_end_phi;end_body_normal=current_end_body_normal;body_velocity=current_body_velocity;}
    return hit_time<FLT_MAX;
}
//##################################################################### 
// Function Push_Out_Point
//##################################################################### 
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Push_Out_Point(TV& X,const T collision_distance,const bool check_particle_crossover,bool& particle_crossover,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    T distance=FLT_MAX,current_distance;TV X_old=X;
    if(objects) for(int k=0;k<objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);TV current_X=X_old;
        if(collision_geometry_collection.bodies(i)->Push_Out_Point(current_X,collision_distance,current_distance) && current_distance<distance){X=current_X;distance=current_distance;}}
    else for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.Size();i++){TV current_X=X_old;
        if(Is_Active(i) && collision_geometry_collection.bodies(i)->Push_Out_Point(current_X,collision_distance,current_distance) && current_distance<distance){X=current_X;distance=current_distance;}}
    if(distance<FLT_MAX){COLLISION_GEOMETRY_ID body_id;int simplex_id;TV intersection_point;
        if(check_particle_crossover) particle_crossover=collision_geometry_collection.Intersection_Between_Points(X_old,X,body_id,simplex_id,intersection_point,objects);
        return true;}
    return false;
}
//##################################################################### 
// Function Occupied_Block
//##################################################################### 
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Occupied_Block(const T_BLOCK& block) const
{
    return occupied_blocks(block.Block());
}
//##################################################################### 
// Function Swept_Occupied_Block
//##################################################################### 
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Swept_Occupied_Block(const T_BLOCK& block) const
{
    return swept_occupied_blocks(block.Block());
}
//#####################################################################
// Function Save_State
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Save_State(const int state_index,const T time)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++)
        if(Is_Active(i))
            collision_geometry_collection.bodies(i)->Save_State(state_index,time);
}
//#####################################################################
// Function Restore_State
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Restore_State(const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++)
        if(Is_Active(i))
            collision_geometry_collection.bodies(i)->Restore_State(state_index);
}
//#####################################################################
// Function Average_States
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Average_States(const int state1, const int state2,const int result_state,const T interpolation_distance)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++)
        if(Is_Active(i))
            collision_geometry_collection.bodies(i)->Average_States(state1,state2,result_state,interpolation_distance);
}
//#####################################################################
// Function Delete_State
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Delete_State(const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++)
        if(Is_Active(i))
            collision_geometry_collection.bodies(i)->Delete_State(state_index);
}
//#####################################################################
// Function Intersection_With_Any_Simplicial_Object
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Intersection_With_Any_Simplicial_Object(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    return collision_geometry_collection.Intersection_With_Any_Simplicial_Object(ray,body_id,objects);
}
//#####################################################################
// Function Add_Bodies
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Add_Bodies(COLLISION_BODY_COLLECTION<TV>& collision_geometry_list)
{
    collision_geometry_collection.Add_Bodies(collision_geometry_list);
}
//#####################################################################
// Function Remove_Body
//#####################################################################
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Remove_Body(COLLISION_GEOMETRY_ID id)
{
    collision_geometry_collection.Remove_Body(id);
}
//#####################################################################
// Function Is_Active
//#####################################################################
template<class TV> bool GRID_BASED_COLLISION_GEOMETRY<TV>::
Is_Active(COLLISION_GEOMETRY_ID id) const
{
    return collision_geometry_collection.Is_Active(id);
}
//##################################################################### 
// Function Read_State
//##################################################################### 
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Read_State(TYPED_ISTREAM input,const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Read_State(input,state_index);
}
//##################################################################### 
// Function Write_State
//##################################################################### 
template<class TV> void GRID_BASED_COLLISION_GEOMETRY<TV>::
Write_State(TYPED_OSTREAM output,const int state_index) const
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Write_State(output,state_index);
}
//#####################################################################
namespace PhysBAM{
template class GRID_BASED_COLLISION_GEOMETRY<VECTOR<float,1> >;
template class GRID_BASED_COLLISION_GEOMETRY<VECTOR<float,2> >;
template class GRID_BASED_COLLISION_GEOMETRY<VECTOR<float,3> >;
template class GRID_BASED_COLLISION_GEOMETRY<VECTOR<double,1> >;
template class GRID_BASED_COLLISION_GEOMETRY<VECTOR<double,2> >;
template class GRID_BASED_COLLISION_GEOMETRY<VECTOR<double,3> >;
}
