//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_COLLISION_BODY_INACCURATE_UNION
//##################################################################### 
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
FLUID_COLLISION_BODY_INACCURATE_UNION(GRID<TV>& grid_input,T contour_value_input)
    :collision_bodies(grid_input),contour_value(contour_value_input),grid(grid_input),levelset(grid_input,phi)
{
    collision_geometries_for_rasterization=&collision_bodies.collision_geometry_collection.bodies;
}    
//#####################################################################
// Function Destructor
//#####################################################################
template<class TV> FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
~FLUID_COLLISION_BODY_INACCURATE_UNION()
{
}
//#####################################################################
// Function Implicit_Geometry_Extended_Value
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Implicit_Geometry_Extended_Value_Helper(const TV& location,UNIFORM_TAG<TV>) const
{
    return interpolation.Clamped_To_Array(grid,phi,location);
}
//#####################################################################
// Function Implicit_Geometry_Extended_Value
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Implicit_Geometry_Extended_Value(const TV& location) const
{
    return Implicit_Geometry_Extended_Value_Helper(location,typename GRID<TV>::GRID_TAG());
}
//#####################################################################
// Function Update_Intersection_Acceleration_Structures
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1,const int state2)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Update_Intersection_Acceleration_Structures(use_swept_simplex_hierarchy,state1,state2);
}
//#####################################################################
// Function Restore_State
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Restore_State(const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Restore_State(state_index);
}
//#####################################################################
// Function Save_State
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Save_State(const int state_index,const T time)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Save_State(state_index,time);
}
//#####################################################################
// Function Read_State
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Read_State(TYPED_ISTREAM input,const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Read_State(input,state_index);
}
//#####################################################################
// Function Write_State
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Write_State(TYPED_OSTREAM output,const int state_index) const
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Write_State(output,state_index);
}
//#####################################################################
// Function Initialize_Grid_Structures
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Initialize_Grid_Structures()
{
    collision_bodies.collision_geometry_collection.Update_Bounding_Boxes();
    // phi and velocity
    phi.Resize(grid.Cell_Indices(3),no_init);
    phi.Fill(10*grid.dX.Min());
    face_velocities.Resize(grid,3,init_all,0);
    face_velocities_set.Resize(grid,3,init_all,false);
    T_FACE_ARRAYS_INT face_velocities_count(grid,3);
    T_FACE_ARRAYS_COLLISION_GEOMETRY_ID face_operations(grid,3);
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++)
        if(collision_bodies.Is_Active(i) && collision_bodies.collision_geometry_collection.bodies(i)->active)
            Initialize_Grid_Structures_Subobject(face_velocities_count,face_operations,i,typename GRID<TV>::GRID_TAG());
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
        if(face_velocities_count.Component(iterator.Axis())(iterator.Face_Index())){
            face_velocities.Component(iterator.Axis())(iterator.Face_Index())/=face_velocities_count.Component(iterator.Axis())(iterator.Face_Index());
            face_velocities_set.Component(iterator.Axis())(iterator.Face_Index())=true;}
}
//#####################################################################
// Function Initialize_Grid_Structures_Subobject
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,T_FACE_ARRAYS_COLLISION_GEOMETRY_ID& face_operations,const COLLISION_GEOMETRY_ID subobject,UNIFORM_TAG<TV>)
{
    COLLISION_GEOMETRY<TV>& collision_body=*collision_bodies.collision_geometry_collection.bodies(subobject);
    RANGE<TV> bounding_box=collision_body.Axis_Aligned_Bounding_Box();
    RANGE<TV_INT> box=grid.Clamp_To_Cell(bounding_box,2).Thickened(1);
    for(CELL_ITERATOR<TV> iterator(grid,box);iterator.Valid();iterator.Next()){
        T phi_value=collision_body.Implicit_Geometry_Extended_Value(iterator.Location());
        phi(iterator.Cell_Index())=min(phi_value,phi(iterator.Cell_Index()));
        if(phi_value<0) for(int axis=0;axis<TV::m;axis++){
                TV_INT face1=iterator.First_Face_Index(axis),face2=iterator.Second_Face_Index(axis);
                if(face_operations.Component(axis)(face1)!=subobject){face_operations.Component(axis)(face1)=subobject;
                    face_velocities.Component(axis)(face1)=collision_body.Pointwise_Object_Velocity(grid.Face(FACE_INDEX<TV::m>(axis,face1)))[axis];
                    face_velocities_count.Component(axis)(face1)++;}
                if(face_operations.Component(axis)(face2)!=subobject){face_operations.Component(axis)(face2)=subobject;
                    face_velocities.Component(axis)(face2)=collision_body.Pointwise_Object_Velocity(grid.Face(FACE_INDEX<TV::m>(axis,face2)))[axis];
                    face_velocities_count.Component(axis)(face2)++;}}}
}
//#####################################################################
// Function Pointwise_Object_Velocity
//#####################################################################
template<class TV> TV FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Pointwise_Object_Velocity(const int aggregate_id,const TV& X) const
{
    T_BLOCK block(grid,X);
    if(!Block_Valid(block,typename GRID<TV>::GRID_TAG())) return TV();
    return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Normalized(block,face_velocities,face_velocities_set,X);
}
//#####################################################################
// Function Pointwise_Object_Velocity
//#####################################################################
template<class TV> TV FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Pointwise_Object_Velocity(const TV& X) const
{
    return Pointwise_Object_Velocity(0,X);
}
//#####################################################################
// Function Pointwise_Object_Pseudo_Velocity
//#####################################################################
template<class TV> TV FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Pointwise_Object_Pseudo_Velocity(const int aggregate_id,const TV& X,const int state1,const int state2) const
{
    return Pointwise_Object_Velocity(0,X);
}
//#####################################################################
// Function Latest_Simplex_Crossover
//#####################################################################
template<class TV> bool FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,
    TV& weights,int& simplex_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const
{
    return false;
}
//#####################################################################
// Function Any_Simplex_Crossover
//#####################################################################
template<class TV> bool FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const
{
    return levelset(end_X)<=0;
}
//#####################################################################
// Function Get_Body_Penetration
//#####################################################################
template<class TV> bool FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Get_Body_Penetration(const TV& start_X,const TV& end_X,const T contour_value,const T dt,T& hit_time,int& simplex_id,T& start_phi,T& end_phi,TV& end_body_normal,
    TV& body_velocity) const
{
    end_phi=levelset(end_X);
    if(end_phi>contour_value) return false;
    start_phi=levelset(start_X);
    end_body_normal=levelset.Normal(end_X);
    body_velocity=Pointwise_Object_Velocity(0,end_X);
    hit_time=dt;
    return true;
}
//#####################################################################
// Function Number_Of_Simplices
//#####################################################################
template<class TV> int FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Number_Of_Simplices() const
{
    return 0;
}
//#####################################################################
// Function Push_Out_Point
//#####################################################################
template<class TV> bool FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Push_Out_Point(TV& X,const T collision_distance,T& distance) const
{
    T current_distance=levelset(X);
    if(current_distance<collision_distance){
        X+=(collision_distance-current_distance)*levelset.Normal(X);
        distance=current_distance;
        return true;}
    return false;
}
//#####################################################################
// Function Inside_Any_Simplex
//#####################################################################
template<class TV> bool FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Inside_Any_Simplex(const TV& location,int& simplex_id) const
{
    return levelset(location)<=contour_value;
}
//#####################################################################
// Function Simplex_Closest_Non_Intersecting_Point
//#####################################################################
template<class TV> bool FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const
{
    auto implicit_object_on_a_ray=[this,&ray](T x){return levelset(ray.Point(x));};

    if(implicit_object_on_a_ray(0)<=0){
        ray.t_max=0;
        return true;}
    if(implicit_object_on_a_ray(ray.t_max)>0) return false;
    ray.t_max=Bisection_Secant_Root(implicit_object_on_a_ray,(T)0,ray.t_max,(T).001*grid.dX.Min());
    ray.t_max=max((T)0,ray.t_max-(T).01*grid.dX.Min()); // TODO: probably make this shift a parameter, or find an entirely different cleaner way
    return true;
}
//#####################################################################
// Function Simplex_Intersection
//#####################################################################
template<class TV> bool FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Simplex_Intersection(RAY<TV>& ray) const
{
    auto implicit_object_on_a_ray=[this,&ray](T x){return levelset(ray.Point(x));};

    if(implicit_object_on_a_ray(0)<=0){
        ray.t_max=0;
        return true;}
    if(implicit_object_on_a_ray(ray.t_max)>0) return false;
    ray.t_max=Bisection_Secant_Root(implicit_object_on_a_ray,(T)0,ray.t_max,(T).001*grid.dX.Min());
    return true;
}
//#####################################################################
// Function Simplex_Closest_Point_On_Boundary
//#####################################################################
template<class TV> TV FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,
    const T thickness_over_2,int* simplex_id,T* returned_distance) const
{
    T distance=levelset(location);
    if(returned_distance) *returned_distance=distance;
    return location-distance*levelset.Normal(location);
}
//##################################################################### 
namespace PhysBAM
{
template class FLUID_COLLISION_BODY_INACCURATE_UNION<VECTOR<float,2> >;
template class FLUID_COLLISION_BODY_INACCURATE_UNION<VECTOR<float,3> >;
template class FLUID_COLLISION_BODY_INACCURATE_UNION<VECTOR<double,2> >;
template class FLUID_COLLISION_BODY_INACCURATE_UNION<VECTOR<double,3> >;
}
