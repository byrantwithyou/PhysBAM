//#####################################################################
// Copyright 2006, Geoffrey Irving, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_COLLISION_BODY_INACCURATE_UNION
//#####################################################################
#ifndef __FLUID_COLLISION_BODY_INACCURATE_UNION__
#define __FLUID_COLLISION_BODY_INACCURATE_UNION__

#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY.h>
namespace PhysBAM{

template<class TV>
class FLUID_COLLISION_BODY_INACCURATE_UNION:public COLLISION_GEOMETRY<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_T;
    typedef ARRAY<int,FACE_INDEX<TV::m> > T_FACE_ARRAYS_INT;
    typedef ARRAY<COLLISION_GEOMETRY_ID,FACE_INDEX<TV::m> > T_FACE_ARRAYS_COLLISION_GEOMETRY_ID;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SIMPLEX;
    typedef LINEAR_INTERPOLATION_MAC_HELPER<TV> T_LINEAR_INTERPOLATION_MAC_HELPER;
    typedef typename GRID<TV>::BLOCK T_BLOCK;
    typedef COLLISION_GEOMETRY<TV> BASE;
public:
    GRID_BASED_COLLISION_GEOMETRY<TV> collision_bodies;
    using BASE::collision_geometries_for_rasterization;
    T contour_value;
    GRID<TV>& grid;
    T_ARRAYS_T phi;
    LEVELSET_IMPLICIT_OBJECT<TV> levelset;
private:
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities;
    ARRAY<bool,FACE_INDEX<TV::m> > face_velocities_set;
    LINEAR_INTERPOLATION_UNIFORM<TV,T> interpolation;
public:

    FLUID_COLLISION_BODY_INACCURATE_UNION(GRID<TV>& grid_input);
    FLUID_COLLISION_BODY_INACCURATE_UNION(GRID<TV>& grid_input,T contour_value_input);
    
    ~FLUID_COLLISION_BODY_INACCURATE_UNION();

    TV Pointwise_Object_Pseudo_Velocity(const int aggregate_id,const TV& X,const int state1,const int state2) const override
    {return Pointwise_Object_Velocity(0,X);}

    // TODO: check whether this is still valid
    //      in the low accuracy case, the points inside the object are set to invalid in the driver, so this does not need to invalidate anything
    //      (i.e., it ignored crossovers)
    // If it isn't the case, we need to generalize COLLISION_BODY_COLLECTION::Latest_Crossover to work for this.
    bool Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,TV& weights,int& simplex_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const override
    {return false;}

    bool Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const override
    {return levelset(end_X)<=0;}

    bool Get_Body_Penetration(const TV& start_X,const TV& end_X,const T contour_value,const T dt,T& hit_time,int& simplex_id,T& start_phi,T& end_phi,TV& end_body_normal,
        TV& body_velocity) const override
    {end_phi=levelset(end_X);if(end_phi>contour_value) return false;
    start_phi=levelset(start_X);end_body_normal=levelset.Normal(end_X);
    body_velocity=Pointwise_Object_Velocity(0,end_X);
    hit_time=dt;return true;}

    int Number_Of_Simplices() const override
    {return 0;}

    bool Push_Out_Point(TV& X,const T collision_distance,T& distance) const override
    {T current_distance=levelset(X);
    if(current_distance<collision_distance){
        X+=(collision_distance-current_distance)*levelset.Normal(X);
        distance=current_distance;return true;}
    return false;}

    bool Inside_Any_Simplex(const TV& location,int& simplex_id) const override
    {return levelset(location)<=contour_value;}

    bool Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const override
    {IMPLICIT_OBJECT_ON_A_RAY<LEVELSET_IMPLICIT_OBJECT<TV> > implicit_object_on_a_ray(levelset,ray);
    if(implicit_object_on_a_ray(0)<=0){ray.t_max=0;return true;}
    if(implicit_object_on_a_ray(ray.t_max)>0) return false;
    ITERATIVE_SOLVER<T> solver;solver.tolerance=(T).001*grid.dX.Min();
    ray.t_max=solver.Bisection_Secant_Root(implicit_object_on_a_ray,0,ray.t_max);
    ray.t_max=max((T)0,ray.t_max-(T).01*grid.dX.Min()); // TODO: probably make this shift a parameter, or find an entirely different cleaner way
    return true;}

    bool Simplex_Intersection(RAY<TV>& ray) const override
    {IMPLICIT_OBJECT_ON_A_RAY<LEVELSET_IMPLICIT_OBJECT<TV> > implicit_object_on_a_ray(levelset,ray);
    if(implicit_object_on_a_ray(0)<=0){ray.t_max=0;return true;}
    if(implicit_object_on_a_ray(ray.t_max)>0) return false;
    ITERATIVE_SOLVER<T> solver;solver.tolerance=(T).001*grid.dX.Min();
    ray.t_max=solver.Bisection_Secant_Root(implicit_object_on_a_ray,0,ray.t_max);
    return true;}

    TV Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2=0,int* simplex_id=0,T* returned_distance=0) const override
    {T distance=levelset(location);if(returned_distance) *returned_distance=distance;return location-distance*levelset.Normal(location);}

private:
    template<class T_TAG> static bool Block_Valid(const T_BLOCK& block,T_TAG)
    {return true;}

public:

    TV Pointwise_Object_Velocity(const int aggregate_id,const TV& X) const override
    {T_BLOCK block(grid,X);if(!Block_Valid(block,typename GRID<TV>::GRID_TAG())) return TV();
    return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Normalized(block,face_velocities,face_velocities_set,X);}

    TV Pointwise_Object_Velocity(const TV& X) const override
    {return Pointwise_Object_Velocity(0,X);}

//#####################################################################
    void Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1=0,const int state2=0) override;
    void Save_State(const int state_index,const T time=0) override;
    void Restore_State(const int state_index) override;
    void Read_State(TYPED_ISTREAM& input,const int state_index) override;
    void Write_State(TYPED_OSTREAM& output,const int state_index) const override;
    typename TV::SCALAR Implicit_Geometry_Extended_Value(const TV& location) const override;
    void Initialize_Grid_Structures(const GRID<TV>& grid,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id) const;
private:
    typename TV::SCALAR Implicit_Geometry_Extended_Value_Helper(const TV& location,UNIFORM_TAG<TV>) const;
    void Initialize_Grid_Structures_Helper(OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id,UNIFORM_TAG<TV>);
    void Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,T_FACE_ARRAYS_COLLISION_GEOMETRY_ID& face_operations,const COLLISION_GEOMETRY_ID subobject,UNIFORM_TAG<TV>);
//#####################################################################
};
}
#endif
