//#####################################################################
// Copyright 2006-2008, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY
//#####################################################################
#ifndef __RIGID_BODY__
#define __RIGID_BODY__

#include <Tools/Data_Structures/ELEMENT_ID.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/FRAME.h>
#include <Tools/Utilities/Find_Type.h>
#include <Tools/Utilities/PHYSBAM_ATTRIBUTE.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>

namespace PhysBAM{

template<class TV>
class RIGID_BODY:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    typedef typename conditional<TV::m==2,SEGMENT_HIERARCHY<TV>,TRIANGLE_HIERARCHY<T> >::type T_SIMPLEX_HIERARCHY;

public:
    IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* implicit_object; // implicit representation of geometry
    T_SIMPLICIAL_OBJECT* simplicial_object; // discrete representation of geometry
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    int particle_index;
    bool is_static; // not saved to file - indicates whether this object is static in the scene
    T surface_roughness; // small number indicating errors in the geometry
    T coefficient_of_friction;
    T coefficient_of_restitution; // not saved to file
    T coefficient_of_rolling_friction; // not saved to file
    bool is_temporarily_static; // not saved to file
    T fracture_threshold;
    bool thin_shell;
    bool CFL_initialized; // true if precalculation is done, false if out of date
    T bounding_box_radius; // needed for CFL calculation

    T_ORIENTED_BOX oriented_box;
    RANGE<TV> axis_aligned_bounding_box;
    bool bounding_box_up_to_date;
    std::string name; // not saved to file - for convenience and debugging.
    T_SIMPLEX_HIERARCHY* moving_simplex_hierarchy; // bounds moving simplices in world space, to accelerate crossover check

    ARRAY<STRUCTURE<TV>*> structures;

    RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,bool create_collision_geometry=false,int index=-1);
    virtual ~RIGID_BODY();

    template<class T_STRUCTURE>
    T_STRUCTURE Find_Structure(const int index=0) const
    {return Find_Type<T_STRUCTURE>(structures,index);}

    virtual std::string Name() const {return Static_Name();}
    static std::string Static_Name()
    {return LOG::sprintf("RIGID_BODY<VECTOR<T,%d> >",TV::m);}

    static T Coefficient_Of_Friction(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1)
    {return min(body0.coefficient_of_friction,body1.coefficient_of_friction);}

    RAY<TV> Object_Space_Ray(const RAY<TV>& world_space_ray) const
    {RAY<TV> transformed_ray(Object_Space_Point(world_space_ray.endpoint),Object_Space_Vector(world_space_ray.direction));
    transformed_ray.semi_infinite=world_space_ray.semi_infinite;transformed_ray.t_max=world_space_ray.t_max;
    transformed_ray.aggregate_id=world_space_ray.aggregate_id;
    return transformed_ray;}

    TV Object_Space_Point(const TV& world_space_point) const
    {return Frame().Inverse_Times(world_space_point);}

    TV Object_Space_Vector(const TV& world_space_vector) const
    {return Frame().r.Inverse_Rotate(world_space_vector);}

    TV World_Space_Point(const TV& object_space_point) const
    {return Frame()*object_space_point;}

    TV World_Space_Vector(const TV& object_space_vector) const
    {return Frame().r.Rotate(object_space_vector);}

    FRAME<TV>& Frame() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.frame(particle_index);}

    const FRAME<TV>& Frame() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.frame(particle_index);}

    TWIST<TV>& Twist() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.twist(particle_index);}

    const TWIST<TV>& Twist() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.twist(particle_index);}

    static TV Pointwise_Object_Velocity(const TWIST<TV>& twist,const TV& t,const TV& X)
    {return twist.linear+TV::Cross_Product(twist.angular,X-t);}

    static TV Relative_Velocity(const RIGID_BODY<TV>& geometry0,const RIGID_BODY<TV>& geometry1,const TV& world_point)
    {return geometry0.Pointwise_Object_Velocity(world_point)-geometry1.Pointwise_Object_Velocity(world_point);}

    static TV Relative_Velocity_At_Geometry1_Particle(const RIGID_BODY<TV>& geometry0,const RIGID_BODY<TV>& geometry1,const TV& world_point,const int particle_index)
    {return geometry0.Pointwise_Object_Velocity_At_Particle(world_point,particle_index)-geometry1.Pointwise_Object_Velocity(world_point);}

    static T_SPIN Relative_Angular_Velocity(const RIGID_BODY<TV>& geometry0,const RIGID_BODY<TV>& geometry1) // make sure the angular velocities are updated before calling this!
    {return geometry0.Twist().angular-geometry1.Twist().angular;}

    static TWIST<TV> Relative_Twist(const RIGID_BODY<TV>& geometry0,const RIGID_BODY<TV>& geometry1,const TV& world_point)
    {return TWIST<TV>(geometry0.Pointwise_Object_Velocity(world_point)-geometry1.Pointwise_Object_Velocity(world_point),geometry0.Twist().angular-geometry1.Twist().angular);}

    const T_ORIENTED_BOX& Oriented_Bounding_Box() const
    {assert(bounding_box_up_to_date);return oriented_box;}

    const RANGE<TV>& Axis_Aligned_Bounding_Box() const
    {assert(bounding_box_up_to_date);return axis_aligned_bounding_box;}

    bool Simplex_Inside(const TV& location,const T collision_thickness) const
    {return simplicial_object->Inside(Object_Space_Point(location),collision_thickness);}

    bool Simplex_Outside(const TV& location,const T collision_thickness) const
    {return simplicial_object->Outside(Object_Space_Point(location),collision_thickness);}

    bool Bounding_Boxes_Intersect(const RIGID_BODY<TV>& rigid_body,const T thickness=0) const
    {if(!thickness) return Axis_Aligned_Bounding_Box().Lazy_Intersection(rigid_body.Axis_Aligned_Bounding_Box()) // check axis aligned first for speed
        && Oriented_Bounding_Box().Intersection(rigid_body.Oriented_Bounding_Box());
    return Axis_Aligned_Bounding_Box().Intersection(rigid_body.Axis_Aligned_Bounding_Box(),thickness)
        && Oriented_Bounding_Box().Thickened(thickness).Intersection(rigid_body.Oriented_Bounding_Box());} // Thickened is rather expensive; avoid it where possible

    TV Simplex_Normal(const TV& location,const int aggregate=0) const
    {return World_Space_Vector(simplicial_object->Normal(Object_Space_Point(location),aggregate));}

    bool Simplex_Boundary(const TV& location,const T collision_thickness) const
    {return simplicial_object->Boundary(Object_Space_Point(location),collision_thickness);}    

    void Set_Mass(const T mass_input) // rescales moments of inertia
    {T& mass=Mass();Inertia_Tensor()*=mass_input/mass;mass=mass_input;}

    void Set_Rigid_Mass(const RIGID_BODY_MASS<TV>& rigid_mass_input)
    {Inertia_Tensor()=rigid_mass_input.inertia_tensor;Mass()=rigid_mass_input.mass;}

    void Set_Inertia_Tensor(const DIAGONAL_MATRIX<T,TV::SPIN::m>& inertia_tensor_input) // already scaled by the mass of the object
    {Inertia_Tensor()=inertia_tensor_input;}

    void Rescale(const T scaling_factor,const bool rescale_mass=true)
    {Inertia_Tensor()*=sqr(scaling_factor);if(rescale_mass) Set_Mass(Mass()*pow(scaling_factor,TV::dimension-thin_shell));}

    T Length_Scale_Squared() const
    {return Inertia_Tensor().Max()/Mass();}

    void Set_Coefficient_Of_Restitution(const T coefficient_input=.5)
    {coefficient_of_restitution=coefficient_input;}

    static T Coefficient_Of_Restitution(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1)
    {return min(body0.coefficient_of_restitution,body1.coefficient_of_restitution);}

    void Set_Coefficient_Of_Rolling_Friction(const T coefficient_input=.5)
    {coefficient_of_rolling_friction=coefficient_input;}

    static T Coefficient_Of_Rolling_Friction(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1)
    {return min(body0.coefficient_of_rolling_friction,body1.coefficient_of_rolling_friction);}

    bool& Is_Kinematic() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.kinematic(particle_index);}
    
    const bool& Is_Kinematic() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.kinematic(particle_index);}

    bool Has_Infinite_Inertia() const
    {return is_static || is_temporarily_static || rigid_body_collection.rigid_body_particles.kinematic(particle_index);}

    bool Is_Simulated() const
    {return !is_static && !rigid_body_collection.rigid_body_particles.kinematic(particle_index);}

    SYMMETRIC_MATRIX<T,TV::SPIN::m> World_Space_Inertia_Tensor() const // relative to the center of mass
    {return Rigid_Mass().World_Space_Inertia_Tensor(Frame().r);}

    SYMMETRIC_MATRIX<T,TV::SPIN::m> World_Space_Inertia_Tensor_Inverse() const // relative to the center of mass
    {return Rigid_Mass().World_Space_Inertia_Tensor_Inverse(Frame().r);};

    SYMMETRIC_MATRIX<T,TV::SPIN::m> World_Space_Inertia_Tensor(const TV& reference_point) const // relative to a reference point
    {return Rigid_Mass().World_Space_Inertia_Tensor(Frame(),reference_point);}

    T_SPIN World_Space_Inertia_Tensor_Times(const T_SPIN& angular_velocity) const
    {return Rigid_Mass().World_Space_Inertia_Tensor_Times(Frame().r,angular_velocity);}

    T_SPIN World_Space_Inertia_Tensor_Inverse_Times(const T_SPIN& angular_momentum) const
    {return Rigid_Mass().World_Space_Inertia_Tensor_Inverse_Times(Frame().r,angular_momentum);}

    void Update_Angular_Velocity() // needs to be called to keep the angular velocity valid
    {Twist().angular=World_Space_Inertia_Tensor_Inverse_Times(Angular_Momentum());}

    void Update_Angular_Velocity(RIGID_BODY_STATE<TV>& state) // needs to be called to keep the angular velocity valid
    {state.Update_Angular_Velocity(Inertia_Tensor());}

    void Update_Angular_Momentum() // assumes a valid angular_velocity
    {Angular_Momentum()=World_Space_Inertia_Tensor_Times(Twist().angular);}

    void Update_Angular_Momentum(RIGID_BODY_STATE<TV>& state) // assumes a valid angular_velocity
    {state.Update_Angular_Momentum(Inertia_Tensor());}

    T_SPIN& Angular_Momentum() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.angular_momentum(particle_index);}

    const T_SPIN& Angular_Momentum() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.angular_momentum(particle_index);}

    TV Momentum() const
    {return Mass()*Twist().linear;}

    T& Mass() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.mass(particle_index);}

    const T& Mass() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.mass(particle_index);}

    DIAGONAL_MATRIX<T,TV::SPIN::m>& Inertia_Tensor() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.inertia_tensor(particle_index);}

    const DIAGONAL_MATRIX<T,TV::SPIN::m>& Inertia_Tensor() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particles.inertia_tensor(particle_index);}

    const RIGID_BODY_MASS<TV> Rigid_Mass() const PHYSBAM_ALWAYS_INLINE
    {return RIGID_BODY_MASS<TV>(Mass(),Inertia_Tensor());}

    T Kinetic_Energy() const
    {return Translational_Kinetic_Energy()+Rotational_Kinetic_Energy();}

    T Translational_Kinetic_Energy() const
    {return (T).5*Mass()*Twist().linear.Magnitude_Squared();}

    T Rotational_Kinetic_Energy() const
    {return Rigid_Mass().Rotational_Kinetic_Energy(Frame().r,Angular_Momentum());}

    T Kinetic_Energy(const TWIST<TV>& twist) const
    {if(Has_Infinite_Inertia()) return 0;return (T).5*Mass()*twist.linear.Magnitude_Squared()+(T).5*TV::SPIN::Dot_Product(World_Space_Inertia_Tensor_Times(twist.angular),twist.angular);}

private:
    static MATRIX<T,1> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,1>& vector,const MATRIX<T,0>& tensor)
    {return MATRIX<T,1>();}

    template<class T_TENSOR,int d>
    static SYMMETRIC_MATRIX<T,d> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,d>& vector,const T_TENSOR& tensor)
    {return tensor.Conjugate_With_Cross_Product_Matrix(vector);}
public:

    SYMMETRIC_MATRIX<T,TV::m> Impulse_Factor(const TV& location) const
    {if(Has_Infinite_Inertia()) return SYMMETRIC_MATRIX<T,TV::m>(); // return zero matrix
    return Conjugate_With_Cross_Product_Matrix(location-Frame().t,World_Space_Inertia_Tensor_Inverse())+1/Mass();}

    SYMMETRIC_MATRIX<T,TV::m> Object_Space_Impulse_Factor(const TV& object_space_location) const
    {if(Has_Infinite_Inertia()) return SYMMETRIC_MATRIX<T,TV::m>(); // return zero matrix
    return Conjugate_With_Cross_Product_Matrix(object_space_location,Inertia_Tensor().Inverse())+1/Mass();}

    static SYMMETRIC_MATRIX<T,TV::m> Impulse_Factor(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1,const TV& location)
    {assert(!body0.Has_Infinite_Inertia() || !body1.Has_Infinite_Inertia());
    return body0.Impulse_Factor(location)+body1.Impulse_Factor(location);}

    TV Simplex_Normal(const int id,const RIGID_BODY_STATE<TV>& state) const
    {return state.World_Space_Vector(simplicial_object->Normal(id));}

    void Save_State(RIGID_BODY_STATE<TV>& state,const T time=0) const
    {state.time=time;state.frame=Frame();state.twist=Twist();state.angular_momentum=Angular_Momentum();}

    void Restore_State(const RIGID_BODY_STATE<TV>& state)
    {Frame()=state.frame;Twist()=state.twist;Angular_Momentum()=state.angular_momentum;}

    TWIST<TV> Gather(const TWIST<TV>& wrench,const TV& location) const
    {return TWIST<TV>(wrench.linear,wrench.angular+TV::Cross_Product(location-Frame().t,wrench.linear));}

    TWIST<TV> Scatter(const TWIST<TV>& twist,const TV& location) const
    {return TWIST<TV>(twist.linear+TV::Cross_Product(twist.angular,location-Frame().t),twist.angular);}

    TWIST<TV> Inertia_Inverse_Times(const TWIST<TV>& wrench) const
    {if(Has_Infinite_Inertia()) return TWIST<TV>();
    return TWIST<TV>(wrench.linear/Mass(),World_Space_Inertia_Tensor_Inverse_Times(wrench.angular));}

    TWIST<TV> Inertia_Times(const TWIST<TV>& twist) const
    {PHYSBAM_ASSERT(!Has_Infinite_Inertia());
    return TWIST<TV>(twist.linear*Mass(),World_Space_Inertia_Tensor_Times(twist.angular));}

//#####################################################################
private:
    static const int mm=TV::m+T_SPIN::m;
public:
    void Effective_Inertia_Inverse(MATRIX<T,mm>& extended_mass_inverse,const TV& location) const;
    void Effective_Inertia(MATRIX<T,mm>& extended_mass_inverse,const TV& location) const;
    void Gather_Matrix(MATRIX<T,mm>& gather,const TV& location) const;
    void Scatter_Matrix(MATRIX<T,mm>& scatter,const TV& location) const;
    TWIST<TV> Effective_Inertia_Inverse_Times(const TWIST<TV>& wrench,const TV& location) const;
    TWIST<TV> Effective_Inertia_Times(const TWIST<TV>& twist,const TV& location) const;
    T Volume() const;
    static void Print_Pairs(const RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const ARRAY<VECTOR<int,2> >& pairs);
    void Initialize_CFL();
    virtual T CFL(const T max_distance_per_time_step,const T max_rotation_per_time_step,const bool verbose=false);
    void Interpolate_Between_States(const RIGID_BODY_STATE<TV>& state1,const RIGID_BODY_STATE<TV>& state2,const T time,RIGID_BODY_STATE<TV>& interpolated_state);
    void Compute_Velocity_Between_States(const RIGID_BODY_STATE<TV>& state1,const RIGID_BODY_STATE<TV>& state2,RIGID_BODY_STATE<TV>& result_state);
    void Apply_Impulse_To_Body(const TV& location,const TV& impulse,const T_SPIN& angular_impulse=T_SPIN(),const bool half_impulse_for_accumulator=false);
    static void Apply_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& impulse,const T_SPIN& angular_impulse=T_SPIN(),
        const bool half_impulse_for_accumulator=false);
    static void Compute_Clamped_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,TWIST<TV>& impulse,const ROTATION<TV>& saved_rotation_1,
        const ROTATION<TV>& saved_rotation_2);
    static TWIST<TV> Compute_Collision_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2,const TV& location,
        const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction=0,const bool clamp_friction_magnitude=true,
        const bool rolling_friction=false,const bool clamp_energy=false);
    static void Apply_Collision_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2,const TV& location,
        const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction=0,const bool clamp_friction_magnitude=true,
        const bool rolling_friction=false,const bool clamp_energy=false,const bool half_impulse_for_accumulator=false);
    static void Apply_Sticking_And_Angular_Sticking_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TWIST<TV>& delta_relative_twist,
        const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix);
    static TWIST<TV> Apply_Rolling_Friction(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& normal,const T normal_impulse);
    static TWIST<TV> Find_Impulse_And_Angular_Impulse(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1,const TV& location,const TWIST<TV>& delta_rel_twist_at_location);
    static TWIST<TV> Find_Impulse_And_Angular_Impulse(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1,const TV& location,const TWIST<TV>& delta_rel_twist_at_location,
        const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix);
    static void Apply_Push(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& normal,const T distance);
    void Apply_Push_To_Body(const TV& location,const TV& impulse,const T_SPIN& angular_impulse=T_SPIN());
    T Volumetric_Density() const;
    void Diagonalize_Inertia_Tensor(const SYMMETRIC_MATRIX<T,TV::SPIN::m>& inertia_tensor_at_center_of_mass);
    template<class T2> void Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface(TETRAHEDRALIZED_VOLUME<T2>& tetrahedralized_volume,
        TRIANGULATED_SURFACE<T>& triangulated_surface,const T cell_size,const int subdivision_loops=0,const bool (*create_levelset_test)(TETRAHEDRALIZED_VOLUME<T>&)=0,
        const bool use_implicit_surface_maker=true,const T shrink_levelset_amount=0);
    void Update_Bounding_Box();
    void Update_Bounding_Box_From_Implicit_Geometry();
    void Add_Structure(STRUCTURE<TV>& structure); // set up acceleration stuctures for certain types of structures
    static void Print_Names(const RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const ARRAY<int>& indices);
    virtual TV Pointwise_Object_Velocity(const TV& X) const;
    virtual TV Pointwise_Object_Velocity_At_Particle(const TV& X,const int particle_index) const;
    TV Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const;
    TV Implicit_Geometry_Normal(const TV& location,const int aggregate=-1) const;
    TV Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const;
    bool Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value=0) const;
    bool Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value=0) const;
    T Implicit_Geometry_Extended_Value(const TV& location) const;
    const RANGE<TV>& Object_Space_Bounding_Box() const;
    bool Simplex_Intersection(RAY<TV>& ray,const T collision_thickness) const;
    RANGE<TV> World_Space_Simplex_Bounding_Box(const int id) const;
    typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX World_Space_Simplex(const int id) const;
    RANGE<TV> World_Space_Simplex_Bounding_Box(const int id,const FRAME<TV>& frame) const;
protected:
    void Remove_Structure(STRUCTURE<TV>* structure);
    template<class T> friend void Add_Structure_Helper(RIGID_BODY<VECTOR<T,1> >& self,STRUCTURE<VECTOR<T,1> >& structure);
    template<class T> friend void Add_Structure_Helper(RIGID_BODY<VECTOR<T,2> >& self,STRUCTURE<VECTOR<T,2> >& structure);
    template<class T> friend void Add_Structure_Helper(RIGID_BODY<VECTOR<T,3> >& self,STRUCTURE<VECTOR<T,3> >& structure);
};
template<class TV,class SCALAR> struct REPLACE_FLOATING_POINT<RIGID_BODY<TV>,SCALAR>{typedef RIGID_BODY<typename REPLACE_FLOATING_POINT<TV,SCALAR>::TYPE> TYPE;};
}
#endif
