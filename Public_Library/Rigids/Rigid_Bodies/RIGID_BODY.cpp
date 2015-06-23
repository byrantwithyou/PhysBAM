//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/PROJECTED_ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Vectors/Dot_Product.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Intersections/RAY_POINT_SIMPLICES_1D_INTERSECTION.h>
#include <Geometry/Intersections/RAY_SEGMENTED_CURVE_2D_INTERSECTION.h>
#include <Geometry/Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
#include <Geometry/Level_Sets/LEVELSET_MAKER.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_HELPER.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Particles/RIGIDS_PARTICLES_FORWARD.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY<TV>::
RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,bool create_collision_geometry,int index)
    :implicit_object(0),simplicial_object(0),rigid_body_collection(rigid_body_collection_input),is_static(false),
    surface_roughness((T)1e-6),coefficient_of_friction((T).5),is_temporarily_static(false),fracture_threshold(FLT_MAX),
    thin_shell(false),CFL_initialized(false),bounding_box_up_to_date(false),moving_simplex_hierarchy(0)
{
    if(index>=0) particle_index=index;
    else{
        FRAME<TV>* old_frames=rigid_body_collection.rigid_body_particles.frame.base_pointer;
        particle_index=rigid_body_collection.rigid_body_particles.Add_Element();
        if(old_frames!=rigid_body_collection.rigid_body_particles.frame.base_pointer)
            rigid_body_collection.Update_Level_Set_Transforms();
        rigid_body_collection.rigid_body_particles.structure_ids(particle_index)=VECTOR<int,3>(-1,-1,-1);}
    if(particle_index>=rigid_body_collection.rigid_body_particles.rigid_body.m)
        rigid_body_collection.rigid_body_particles.rigid_body.Resize(particle_index+1);
    else assert(!rigid_body_collection.rigid_body_particles.rigid_body(particle_index));
    rigid_body_collection.rigid_body_particles.rigid_body(particle_index)=this;

    if(create_collision_geometry && rigid_body_collection.collision_body_list)
        rigid_body_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(*this),particle_index,true);

    Set_Rigid_Mass(RIGID_BODY_MASS<TV>::Identity_Mass());
    Angular_Momentum()=T_SPIN();
    rigid_body_collection.rigid_body_particles.kinematic(particle_index)=false;
    Set_Coefficient_Of_Restitution();
    Set_Coefficient_Of_Rolling_Friction();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY<TV>::
~RIGID_BODY()
{
    rigid_body_collection.collision_body_list->Remove_Body(rigid_body_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(particle_index));
    rigid_body_collection.rigid_body_particles.structure_ids(particle_index)=VECTOR<int,3>(-1,-1,-1);
    rigid_body_collection.rigid_body_particles.rigid_body(particle_index)=0;
    delete moving_simplex_hierarchy;
    delete implicit_object;
}
//#####################################################################
// Function Pointwise_Object_Velocity
//#####################################################################
template<class TV> TV RIGID_BODY<TV>::
Pointwise_Object_Velocity(const TV& X) const
{
    return Pointwise_Object_Velocity(Twist(),Frame().t,X);
}
//#####################################################################
// Function Pointwise_Object_Velocity_At_Particle
//#####################################################################
template<class TV> TV RIGID_BODY<TV>::
Pointwise_Object_Velocity_At_Particle(const TV& X,const int particle_index) const
{
    return Pointwise_Object_Velocity(X);
}
//#####################################################################
// Function Implicit_Geometry_Extended_Value
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY<TV>::
Implicit_Geometry_Extended_Value(const TV& location) const
{
    return implicit_object->Extended_Phi(location);
}
//#####################################################################
// Function Implicit_Geometry_Normal
//#####################################################################
template<class TV> TV RIGID_BODY<TV>::
Implicit_Geometry_Normal(const TV& location,const int aggregate) const
{
    return implicit_object->Normal(location,aggregate);
}
//#####################################################################
// Function Implicit_Geometry_Normal
//#####################################################################
template<class TV> TV RIGID_BODY<TV>::
Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
{
    phi_value=(*implicit_object)(location);
    return implicit_object->Normal(location,aggregate);
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside
//#####################################################################
template<class TV> bool RIGID_BODY<TV>::
Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value) const
{
    return implicit_object->Lazy_Inside(location,contour_value);
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool RIGID_BODY<TV>::
Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value) const
{
    return implicit_object->Lazy_Inside_And_Value(location,phi,contour_value);
}
//#####################################################################
// Function Implicit_Geometry_Extended_Normal
//#####################################################################
template<class TV> TV RIGID_BODY<TV>::
Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
{
    phi_value=implicit_object->Extended_Phi(location);
    return implicit_object->Extended_Normal(location,aggregate);
}
//#####################################################################
// Function Simplex_Intersection
//#####################################################################
template<class TV> bool RIGID_BODY<TV>::
Simplex_Intersection(RAY<TV>& ray,const T collision_thickness) const
{
    RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,*simplicial_object,collision_thickness)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    return false;
}
//#####################################################################
// Function Add_Structure_Helper
//#####################################################################
template<class T> void
Add_Structure_Helper(RIGID_BODY<VECTOR<T,1> >& self,STRUCTURE<VECTOR<T,1> >& structure)
{
    typedef VECTOR<T,1> TV;
    if(POINT_SIMPLICES_1D<T>* point_simplices=dynamic_cast<POINT_SIMPLICES_1D<T>*>(&structure)){ // set up acceleration stuctures too
        if(self.simplicial_object) self.Remove_Structure(self.simplicial_object);
        self.simplicial_object=point_simplices;
        if(!self.simplicial_object->bounding_box) self.simplicial_object->Update_Bounding_Box();}
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(self.implicit_object) self.Remove_Structure(self.implicit_object->object_space_implicit_object);
        self.implicit_object=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(implicit_object_input,false,&self.Frame());
        self.implicit_object->Update_Box();self.implicit_object->Compute_Cell_Minimum_And_Maximum(false);} // don't recompute cell min/max if already computed
    self.structures.Append(&structure);
}
template<class T> void
Add_Structure_Helper(RIGID_BODY<VECTOR<T,2> >& self,STRUCTURE<VECTOR<T,2> >& structure)
{
    typedef VECTOR<T,2> TV;
    if(SEGMENTED_CURVE_2D<T>* segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>*>(&structure)){ // set up acceleration stuctures too
        if(self.simplicial_object) self.Remove_Structure(self.simplicial_object);
        self.simplicial_object=segmented_curve;
        if(!self.simplicial_object->bounding_box) self.simplicial_object->Update_Bounding_Box();}
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(self.implicit_object) self.Remove_Structure(self.implicit_object->object_space_implicit_object);
        self.implicit_object=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(implicit_object_input,false,&self.Frame());
        self.implicit_object->Update_Box();self.implicit_object->Compute_Cell_Minimum_And_Maximum(false);} // don't recompute cell min/max if already computed
    else if(TRIANGULATED_AREA<T>* triangulated_area=dynamic_cast<TRIANGULATED_AREA<T>*>(&structure)){
        if(TRIANGULATED_AREA<T>* old_triangulated_area=self.template Find_Structure<TRIANGULATED_AREA<T>*>()) self.Remove_Structure(old_triangulated_area);
        if(!triangulated_area->bounding_box) triangulated_area->Update_Bounding_Box();}
    self.structures.Append(&structure);
}
template<class T> void
Add_Structure_Helper(RIGID_BODY<VECTOR<T,3> >& self,STRUCTURE<VECTOR<T,3> >& structure)
{
    typedef VECTOR<T,3> TV;
    if(TRIANGULATED_SURFACE<T>* triangulated_surface_input=dynamic_cast<TRIANGULATED_SURFACE<T>*>(&structure)){ // set up acceleration stuctures too
        if(self.simplicial_object) self.Remove_Structure(self.simplicial_object);
        self.simplicial_object=triangulated_surface_input;
        if(!self.simplicial_object->hierarchy){self.simplicial_object->Initialize_Hierarchy();self.simplicial_object->hierarchy->Update_Box_Radii();}
        if(!self.simplicial_object->bounding_box) self.simplicial_object->Update_Bounding_Box();}
    else if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(&structure)){
        if(self.implicit_object) self.Remove_Structure(self.implicit_object->object_space_implicit_object);
        self.implicit_object=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(implicit_object_input,false,&self.Frame());
        self.implicit_object->Update_Box();self.implicit_object->Compute_Cell_Minimum_And_Maximum(false);} // don't recompute cell min/max if already computed
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(self.implicit_object) self.Remove_Structure(self.implicit_object->object_space_implicit_object);
        self.implicit_object=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(implicit_object_input,false,&self.Frame());
        self.implicit_object->Update_Box();self.implicit_object->Compute_Cell_Minimum_And_Maximum(false);} // don't recompute cell min/max if already computed
    else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(&structure)){
        if(TETRAHEDRALIZED_VOLUME<T>* old_tetrahedralized_volume=self.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>()) self.Remove_Structure(old_tetrahedralized_volume);
        if(!tetrahedralized_volume->bounding_box) tetrahedralized_volume->Update_Bounding_Box();}
    self.structures.Append(&structure);
}
//#####################################################################
// Function Add_Structure
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Add_Structure(STRUCTURE<TV>& structure)
{
    Add_Structure_Helper(dynamic_cast<RIGID_BODY<TV>&>(*this),structure);
}
//#####################################################################
// Function Print_Names
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Print_Names(const RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const ARRAY<int>& ids)
{
    LOG::cout<<"{";
    for(int i=0;i<ids.m;i++){LOG::cout<<"\""<<rigid_body_collection.rigid_body_particles.rigid_body(ids(i))->name<<"\"";if(i<ids.m) LOG::cout<<", ";}
    LOG::cout<<"}";
}
//#####################################################################
// Function Object_Space_Bounding_Box
//#####################################################################
template<class TV> const RANGE<TV>& RIGID_BODY<TV>::
Object_Space_Bounding_Box() const // at least one of triangulated surface or implicit surface must exist for this to work
{
    if(simplicial_object){PHYSBAM_ASSERT(simplicial_object->bounding_box);return *simplicial_object->bounding_box;}
    else{PHYSBAM_ASSERT(implicit_object);return implicit_object->object_space_implicit_object->Box();}
}
//#####################################################################
// Function Remove_Structure
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Remove_Structure(STRUCTURE<TV>* structure)
{
    structures.Remove_Index_Lazy(structures.Find(structure));
}
//#####################################################################
// Function World_Space_Simplex_Bounding_Box
//#####################################################################
template<class TV> RANGE<TV> RIGID_BODY<TV>::
World_Space_Simplex_Bounding_Box(const int id) const
{
    const VECTOR<int,TV::dimension>& elements=simplicial_object->mesh.elements(id);
    VECTOR<TV,TV::dimension> pts;
    for(int i=0;i<TV::dimension;i++) pts(i)=World_Space_Point(simplicial_object->particles.X(elements(i)));
    return RANGE<TV>::Bounding_Box(pts);
}
//#####################################################################
// Function World_Space_Simplex
//#####################################################################
template<class TV> typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX RIGID_BODY<TV>::
World_Space_Simplex(const int id) const
{
    const VECTOR<int,TV::dimension>& elements=simplicial_object->mesh.elements(id);
    VECTOR<TV,TV::dimension> pts;
    for(int i=0;i<TV::dimension;i++) pts(i)=World_Space_Point(simplicial_object->particles.X(elements(i)));
    return typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX(pts);
}
//#####################################################################
// Function World_Space_Simplex_Bounding_Box
//#####################################################################
template<class TV> RANGE<TV> RIGID_BODY<TV>::
World_Space_Simplex_Bounding_Box(const int id,const FRAME<TV>& frame) const
{
    const VECTOR<int,TV::dimension>& elements=simplicial_object->mesh.elements(id);
    VECTOR<TV,TV::dimension> pts;
    for(int i=0;i<TV::dimension;i++) pts(i)=frame*simplicial_object->particles.X(elements(i));
    return RANGE<TV>::Bounding_Box(pts);
}
template<class TV> void RIGID_BODY<TV>::
Update_Bounding_Box()
{
    if(bounding_box_up_to_date && is_static) return;
    bounding_box_up_to_date=true;
    const RANGE<TV>& box=Object_Space_Bounding_Box();
    oriented_box=T_ORIENTED_BOX(box,Frame());
    axis_aligned_bounding_box=oriented_box.Axis_Aligned_Bounding_Box();
}
template<class TV> void RIGID_BODY<TV>::
Update_Bounding_Box_From_Implicit_Geometry()
{
    assert(implicit_object);
    bounding_box_up_to_date=true;
    oriented_box=T_ORIENTED_BOX(implicit_object->object_space_implicit_object->box,Frame());
    axis_aligned_bounding_box=oriented_box.Axis_Aligned_Bounding_Box();
}
//#####################################################################
// Function Compute_Velocity_Between_States
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Compute_Velocity_Between_States(const RIGID_BODY_STATE<TV>& state1,const RIGID_BODY_STATE<TV>& state2,RIGID_BODY_STATE<TV>& result_state)
{
    RIGID_BODY_STATE<TV>::Compute_Velocity_Between_States(state1,state2,result_state);
    Update_Angular_Momentum(result_state); // Assumes result_state has a valid orientation
}
//#####################################################################
// Function Interpolate_Between_States
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Interpolate_Between_States(const RIGID_BODY_STATE<TV>& state1,const RIGID_BODY_STATE<TV>& state2,const T time,RIGID_BODY_STATE<TV>& interpolated_state)
{
    T alpha=(time-state1.time)/(state2.time-state1.time);alpha=clamp(alpha,(T)0,(T)1);
    PHYSBAM_ASSERT(((time>=state1.time && time<=state2.time) || (time<=state1.time && time>=state2.time)) && state1.time<state2.time);
    interpolated_state.frame.t=state1.frame.t+alpha*(state2.frame.t-state1.frame.t);
    interpolated_state.frame.r=ROTATION<TV>::Spherical_Linear_Interpolation(state1.frame.r,state2.frame.r,alpha);
    interpolated_state.twist.linear=(1-alpha)*state1.twist.linear+alpha*state2.twist.linear;
    interpolated_state.time=time;
    interpolated_state.angular_momentum=(1-alpha)*state1.angular_momentum+alpha*state2.angular_momentum;
    Update_Angular_Velocity(interpolated_state);
}
//#####################################################################
// Function Print_Pairs
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Print_Pairs(const RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const ARRAY<VECTOR<int,2> >& pairs)
{
    LOG::cout<<"{";
    for(int i=0;i<pairs.m;i++){
        LOG::cout<<"(\""<<rigid_body_collection.Rigid_Body(pairs(i)(0)).name<<"\", \""<<rigid_body_collection.Rigid_Body(pairs(i)(1)).name<<"\")";if(i<pairs.m) LOG::cout<<", ";}
    LOG::cout<<"}";
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Initialize_CFL()
{
    // Assumes bounding box is up to date
    const RANGE<TV>& box=Object_Space_Bounding_Box();
    bounding_box_radius=TV::Componentwise_Max(abs(box.Minimum_Corner()),abs(box.Maximum_Corner())).Magnitude();
    CFL_initialized=true;
}
//#####################################################################
// Function CFL
//#####################################################################
// set max_distance_per_time_step and/or max_rotation_per_time_step to zero to ignore those checks
template<class TV> typename TV::SCALAR RIGID_BODY<TV>::
CFL(const T max_distance_per_time_step,const T max_rotation_per_time_step,const bool verbose)
{
    if(is_static) return FLT_MAX;
    if(!CFL_initialized) Initialize_CFL();

    T angular_speed=Twist().angular.Magnitude(),linear_speed=Twist().linear.Magnitude()+bounding_box_radius*angular_speed;
    T linear_dt=max_distance_per_time_step?Robust_Divide(max_distance_per_time_step,linear_speed):FLT_MAX;
    T angular_dt=max_rotation_per_time_step?Robust_Divide(max_rotation_per_time_step,angular_speed):FLT_MAX;

    T dt=min(linear_dt,angular_dt);
    return dt;
}
//#####################################################################
// Function Apply_Impulse_To_Body
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Impulse_To_Body(const TV& location,const TV& impulse,const T_SPIN& angular_impulse,const bool half_impulse_for_accumulator)
{
    if(Has_Infinite_Inertia()) return;
    T_SPIN total_angular_impulse=TV::Cross_Product(location-Frame().t,impulse)+angular_impulse;
    Twist().linear+=impulse/Mass();
    Angular_Momentum()+=total_angular_impulse;
    Update_Angular_Velocity();
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& impulse,const T_SPIN& angular_impulse,const bool half_impulse_for_accumulator)
{
    body0.Apply_Impulse_To_Body(location,impulse,angular_impulse,half_impulse_for_accumulator);
    body1.Apply_Impulse_To_Body(location,-impulse,-angular_impulse,half_impulse_for_accumulator);
}
//#####################################################################
// Function Apply_Clamped_Impulse
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Compute_Clamped_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,TWIST<TV>& impulse,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2)
{
    RIGID_BODY<TV>* bodies[2]={&body0,&body1};
    const ROTATION<TV>* saved_rotation[2]={&saved_rotation_1,&saved_rotation_2};
    T_SPIN jr[2],I_inverse_jr[2];
    T impulse_ratio_num[2]={0},impulse_ratio_denom[2]={0},impulse_linear_magnitude_squared=impulse.linear.Magnitude_Squared();
    TV velocity_adjustment;
    for(int i=0;i<2;i++) if(bodies[i]->rigid_body_collection.rigid_body_particles.kinematic(bodies[i]->particle_index)) velocity_adjustment=bodies[i]->Pointwise_Object_Velocity(location);
    for(int i=0;i<2;i++) if(!bodies[i]->Has_Infinite_Inertia()){
        jr[i]=TV::Cross_Product(location-bodies[i]->Frame().t,impulse.linear)+impulse.angular;
        I_inverse_jr[i]=bodies[i]->Rigid_Mass().World_Space_Inertia_Tensor_Inverse_Times(*saved_rotation[i],jr[i]);
        impulse_ratio_num[i]=Dot_Product(impulse.linear,bodies[i]->Twist().linear-velocity_adjustment)+Dot_Product(I_inverse_jr[i],bodies[i]->Angular_Momentum());
        impulse_ratio_denom[i]=impulse_linear_magnitude_squared/bodies[i]->Mass()+Dot_Product(jr[i],I_inverse_jr[i]);}
    T total_impulse_ratio_num=-2*(impulse_ratio_num[0]-impulse_ratio_num[1]);
    if(total_impulse_ratio_num>0){
        T total_impulse_ratio_denom=impulse_ratio_denom[0]+impulse_ratio_denom[1];
        T scale=total_impulse_ratio_num<total_impulse_ratio_denom?total_impulse_ratio_num/total_impulse_ratio_denom:1;
        impulse*=scale;}
    else impulse=TWIST<TV>();
}
//#####################################################################
// Function Compute_Collision_Impulse
//#####################################################################
// clamp friction magnitude: should be (true) in the elastic collision case, (false) in the inelastic collision case
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Compute_Collision_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2,const TV& location,const TV& normal,
    const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction,const bool clamp_friction_magnitude,const bool rolling_friction,const bool clamp_energy)
{
    if(body0.Has_Infinite_Inertia() && body1.Has_Infinite_Inertia()) return TWIST<TV>();
    TWIST<TV> impulse;
    bool sticking_impulse;
    impulse.linear=PhysBAM::Compute_Collision_Impulse(normal,Impulse_Factor(body0,body1,location),relative_velocity,coefficient_of_restitution,coefficient_of_friction,&sticking_impulse);
    if(rolling_friction && sticking_impulse) impulse+=Apply_Rolling_Friction(body0,body1,location,normal,impulse.linear.Dot(normal));
    if(clamp_energy) Compute_Clamped_Impulse(body0,body1,location,impulse,saved_rotation_1,saved_rotation_2);
    return impulse;
}
//#####################################################################
// Function Apply_Collision_Impulse
//#####################################################################
// clamp friction magnitude: should be (true) in the elastic collision case, (false) in the inelastic collision case
template<class TV> void RIGID_BODY<TV>::
Apply_Collision_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2,const TV& location,const TV& normal,
    const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction,const bool clamp_friction_magnitude,const bool rolling_friction,const bool clamp_energy,
    const bool half_impulse_for_accumulator)
{
    TWIST<TV> impulse=Compute_Collision_Impulse(body0,body1,saved_rotation_1,saved_rotation_2,location,normal,relative_velocity,coefficient_of_restitution,
        coefficient_of_friction,clamp_friction_magnitude,rolling_friction,clamp_energy);
    Apply_Impulse(body0,body1,location,impulse.linear,impulse.angular,half_impulse_for_accumulator);
}
//#####################################################################
// Function Apply_Sticking_And_Angular_Sticking_Impulse
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Sticking_And_Angular_Sticking_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TWIST<TV>& delta_relative_twist,const MATRIX_MXN<T>& angular_constraint_matrix,
    const MATRIX_MXN<T>& prismatic_constraint_matrix)
{
    //body0.Update_Angular_Velocity();body1.Update_Angular_Velocity();
    TWIST<TV> impulse;

    // faster version for fully constrained prismatic and angular impulse
    if(angular_constraint_matrix.Columns()==T_SPIN::dimension && prismatic_constraint_matrix.Columns()==TV::dimension)
        impulse=Find_Impulse_And_Angular_Impulse(body0,body1,location,delta_relative_twist);
    else if(angular_constraint_matrix.Columns()==0 && prismatic_constraint_matrix.Columns()==TV::dimension)
        impulse.linear=Impulse_Factor(body0,body1,location).Inverse()*delta_relative_twist.linear;
    else impulse=Find_Impulse_And_Angular_Impulse(body0,body1,location,delta_relative_twist,angular_constraint_matrix,prismatic_constraint_matrix);
    Apply_Impulse(body0,body1,location,impulse.linear,impulse.angular);
}
//#####################################################################
// Function Apply_Rolling_Friction
//#####################################################################
// location is a point in world space about which body is rolling
template<class T,class TV> static TWIST<TV> Apply_Rolling_Friction_Helper(RIGID_BODY<VECTOR<T,1> >& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& normal,const T normal_impulse)
{
    return TWIST<TV>(); // TODO: implement
}
template<class T,class TV> static TWIST<TV> Apply_Rolling_Friction_Helper(RIGID_BODY<VECTOR<T,2> >& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& normal,const T normal_impulse)
{
    return TWIST<TV>(); // TODO: implement
}
template<class T,class TV> static TWIST<TV> Apply_Rolling_Friction_Helper(RIGID_BODY<VECTOR<T,3> >& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& normal,const T normal_impulse)
{
    PHYSBAM_ASSERT(!body0.Has_Infinite_Inertia() || !body1.Has_Infinite_Inertia());PHYSBAM_ASSERT(normal_impulse>=0);
    T coefficient_of_rolling_friction=RIGID_BODY<TV>::Coefficient_Of_Rolling_Friction(body0,body1);if(!coefficient_of_rolling_friction) return TWIST<TV>();
    body0.Update_Angular_Velocity();body1.Update_Angular_Velocity();
    TV relative_angular_velocity=RIGID_BODY<TV>::Relative_Angular_Velocity(body0,body1);
    T normal_component=TV::Dot_Product(relative_angular_velocity,normal),normal_magnitude=abs(normal_component);
    TV tangential_component=relative_angular_velocity-normal_component*normal;T tangential_magnitude=tangential_component.Magnitude();
    TV tangential_direction;if(tangential_magnitude!=0) tangential_direction=tangential_component/tangential_magnitude;
    normal_magnitude-=coefficient_of_rolling_friction*normal_impulse;tangential_magnitude-=coefficient_of_rolling_friction*normal_impulse;
    TV new_relative_angular_velocity=sign(normal_component)*max((T)0,normal_magnitude)*normal+max((T)0,tangential_magnitude)*tangential_direction;
    return RIGID_BODY<TV>::Find_Impulse_And_Angular_Impulse(body0,body1,location,TWIST<TV>(TV(),new_relative_angular_velocity-relative_angular_velocity));
}
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Apply_Rolling_Friction(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& normal,const T normal_impulse)
{
    return Apply_Rolling_Friction_Helper(body0,body1,location,normal,normal_impulse);
}
//#####################################################################
// Function Find_Impulse_And_Angular_Impulse
//#####################################################################
template<class T,class TV> TWIST<TV>
Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<VECTOR<T,1> >& body0,const RIGID_BODY<VECTOR<T,1> >& body1,const TV& location,const TWIST<TV>& delta_relative_twist_at_location,
    const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class T,class TV> typename enable_if<(TV::m>1),TWIST<TV> >::type
Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1,const TV& location,const TWIST<TV>& delta_relative_twist_at_location,
    const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix)
{
    // compute blocks of constrained matrix
    SYMMETRIC_MATRIX<T,TV::SPIN::m> I_inverse_1,I_inverse_2;T m1_inv_plus_m2_inv=0;
    if(!body0.Has_Infinite_Inertia()){I_inverse_1=body0.World_Space_Inertia_Tensor_Inverse();m1_inv_plus_m2_inv=1/body0.Mass();}
    if(!body1.Has_Infinite_Inertia()){I_inverse_2=body1.World_Space_Inertia_Tensor_Inverse();m1_inv_plus_m2_inv+=1/body1.Mass();}

    // fill in NXN constrained matrix C
    TV r1=location-body0.Frame().t,r2=location-body1.Frame().t;
    int angular_constrained_axes=angular_constraint_matrix.Columns(),prismatic_constrained_axes=prismatic_constraint_matrix.Columns();
    MATRIX_MXN<T> r_cross_P_1=prismatic_constraint_matrix.Cross_Product_Matrix_Times(r1),r_cross_P_2=prismatic_constraint_matrix.Cross_Product_Matrix_Times(r2);
    MATRIX_MXN<T> P_T_r_cross_T_I_inverse_1=r_cross_P_1.Transpose_Times(I_inverse_1),P_T_r_cross_T_I_inverse_2=r_cross_P_2.Transpose_Times(I_inverse_2);
    MATRIX_MXN<T> C(prismatic_constrained_axes+angular_constrained_axes);
    C.Set_Submatrix(0,0,P_T_r_cross_T_I_inverse_1*r_cross_P_1+P_T_r_cross_T_I_inverse_2*r_cross_P_2+m1_inv_plus_m2_inv*prismatic_constraint_matrix.Transpose_Times(prismatic_constraint_matrix));
    if(angular_constrained_axes){
        C.Set_Submatrix(prismatic_constrained_axes,prismatic_constrained_axes,angular_constraint_matrix.Transpose_Times((I_inverse_1+I_inverse_2)*angular_constraint_matrix));
        MATRIX_MXN<T> c01=(P_T_r_cross_T_I_inverse_1+P_T_r_cross_T_I_inverse_2)*angular_constraint_matrix;
        C.Set_Submatrix(0,prismatic_constrained_axes,c01);C.Set_Submatrix(prismatic_constrained_axes,0,c01.Transposed());}

    ARRAY<T> b(prismatic_constraint_matrix.Transpose_Times(delta_relative_twist_at_location.linear));
    if(angular_constrained_axes) b.Append_Elements(angular_constraint_matrix.Transpose_Times(delta_relative_twist_at_location.angular));

    TWIST<TV> impulse;
    ARRAY<T> x=C.Cholesky_Solve(b),linear(x.Array_View(0,prismatic_constrained_axes));
    impulse.linear=prismatic_constraint_matrix*linear;
    if(angular_constrained_axes){
        ARRAY<T> angular(x.Array_View(prismatic_constrained_axes,angular_constrained_axes));
        impulse.angular=angular_constraint_matrix*angular;}
    return impulse;
}
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Find_Impulse_And_Angular_Impulse(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1,const TV& location,const TWIST<TV>& delta_relative_twist_at_location,
    const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix)
{
    PHYSBAM_ASSERT(!body0.Has_Infinite_Inertia() || !body1.Has_Infinite_Inertia());
    return Find_Impulse_And_Angular_Impulse_Helper(body0,body1,location,delta_relative_twist_at_location,angular_constraint_matrix,prismatic_constraint_matrix);
}
//#####################################################################
// Function Find_Impulse_And_Angular_Impulse
//#####################################################################
template<class T,class TV> static TWIST<TV> Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<VECTOR<T,1> >& body0,const RIGID_BODY<TV>& body1,const TV& location,
    const TWIST<TV>& delta_relative_twist_at_location)
{
    return delta_relative_twist_at_location;
}
template<class T,class TV> static TWIST<TV> Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<VECTOR<T,2> >& body0,const RIGID_BODY<TV>& body1,const TV& location,
    const TWIST<TV>& delta_relative_twist_at_location)
{
    SYMMETRIC_MATRIX<T,1> I_inverse_1=body0.World_Space_Inertia_Tensor_Inverse(),I_inverse_2=body1.World_Space_Inertia_Tensor_Inverse();
    SYMMETRIC_MATRIX<T,2> c00=I_inverse_1.Conjugate_With_Cross_Product_Matrix(location-body0.Frame().t)+
        I_inverse_2.Conjugate_With_Cross_Product_Matrix(location-body1.Frame().t)+1/body0.Mass()+1/body1.Mass();
    MATRIX<T,1,2> c01=I_inverse_1.Times_Cross_Product_Matrix(location-body0.Frame().t)+I_inverse_2.Times_Cross_Product_Matrix(location-body1.Frame().t);
    SYMMETRIC_MATRIX<T,1> c11=I_inverse_1+I_inverse_2;
    SYMMETRIC_MATRIX<T,3> A(c00.x00,c00.x10,c01(0,0),c00.x11,c01(0,1),c11.x00);

    VECTOR<T,3> b(delta_relative_twist_at_location.linear.x,delta_relative_twist_at_location.linear.y,delta_relative_twist_at_location.angular.x);
    VECTOR<T,3> all_impulses=A.Inverse()*b;
    return TWIST<TV>(VECTOR<T,2>(all_impulses.x,all_impulses.y),VECTOR<T,1>(all_impulses.z));
}
template<class T,class TV> static TWIST<TV> Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<VECTOR<T,3> >& body0,const RIGID_BODY<TV>& body1,const TV& location,
    const TWIST<TV>& delta_relative_twist_at_location)
{
    SYMMETRIC_MATRIX<T,3> I_inverse_1=body0.World_Space_Inertia_Tensor_Inverse(),I_inverse_2=body1.World_Space_Inertia_Tensor_Inverse();
    MATRIX<T,3> r_cross_1=MATRIX<T,3>::Cross_Product_Matrix(location-body0.Frame().t),r_cross_I_inverse_1=r_cross_1*I_inverse_1,
                r_cross_2=MATRIX<T,3>::Cross_Product_Matrix(location-body1.Frame().t),r_cross_I_inverse_2=r_cross_2*I_inverse_2;
    SYMMETRIC_MATRIX<T,3> c00=SYMMETRIC_MATRIX<T,3>::Times_Transpose_With_Symmetric_Result(r_cross_I_inverse_1,r_cross_1)+
                              SYMMETRIC_MATRIX<T,3>::Times_Transpose_With_Symmetric_Result(r_cross_I_inverse_2,r_cross_2)+1/body0.Mass()+1/body1.Mass();
    MATRIX<T,3> c01=-(r_cross_I_inverse_1+r_cross_I_inverse_2);
    SYMMETRIC_MATRIX<T,3> c11=I_inverse_1+I_inverse_2,c22_inverse=c11.Inverse();
    MATRIX<T,3> c12_c22_inverse=c01*c22_inverse;
    SYMMETRIC_MATRIX<T,3> A=c00-SYMMETRIC_MATRIX<T,3>::Times_Transpose_With_Symmetric_Result(c12_c22_inverse,c01);
    TV b=delta_relative_twist_at_location.linear-c12_c22_inverse*delta_relative_twist_at_location.angular;
    TV linear_impulse=A.Inverse()*b;
    TV angular_impulse=c22_inverse*(delta_relative_twist_at_location.angular-c01.Transpose_Times(linear_impulse));
    return TWIST<TV>(linear_impulse,angular_impulse);
}
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Find_Impulse_And_Angular_Impulse(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1,const TV& location,const TWIST<TV>& delta_relative_twist_at_location)
{
    PHYSBAM_ASSERT(!body0.Has_Infinite_Inertia() || !body1.Has_Infinite_Inertia());

    if(body0.Has_Infinite_Inertia() || body1.Has_Infinite_Inertia()){
        const RIGID_BODY<TV>& body=(!body0.Has_Infinite_Inertia())?body0:body1;TV r(location-body.Frame().t);
        TV impulse_linear=body.Mass()*(delta_relative_twist_at_location.linear+TV::Cross_Product(r,delta_relative_twist_at_location.angular));
        T_SPIN impulse_angular=body.World_Space_Inertia_Tensor()*delta_relative_twist_at_location.angular-TV::Cross_Product(r,impulse_linear);
        return TWIST<TV>(impulse_linear,impulse_angular);}

    // c00*impulse.linear+c01*impulse.angular=delta_relative_twist_at_location.linear
    // c10*impulse.linear+c11*impulse.angular=delta_relative_twist_at_location.angular
    // Note: c10=c01^T
    return Find_Impulse_And_Angular_Impulse_Helper(body0,body1,location,delta_relative_twist_at_location);
}
//#####################################################################
// Function Apply_Push
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Push_To_Body(const TV& location,const TV& impulse,const T_SPIN& angular_impulse)
{
    if(Has_Infinite_Inertia()) return;
    TV velocity=impulse/Mass();
    T_SPIN angular_velocity=World_Space_Inertia_Tensor_Inverse()*TV::Cross_Product(location-Frame().t,impulse)+angular_impulse;
    Frame().t+=velocity;
    Frame().r=ROTATION<TV>::From_Rotation_Vector(angular_velocity)*Frame().r;Frame().r.Normalize();
    Update_Angular_Velocity();
}
//#####################################################################
// Function Apply_Push
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Push(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& normal,const T distance)
{
    if(body0.Has_Infinite_Inertia() && body1.Has_Infinite_Inertia()) return;
    TV impulse=Impulse_Factor(body0,body1,location).Inverse()*(distance*normal);
    body0.Apply_Push_To_Body(location,impulse);
    body1.Apply_Push_To_Body(location,-impulse);
}
//#####################################################################
// Function Volume
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY<TV>::
Volume() const
{
    PHYSBAM_ASSERT(simplicial_object);return simplicial_object->Volumetric_Volume();
}
//#####################################################################
// Function Volumetric_Density
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY<TV>::
Volumetric_Density() const
{
    PHYSBAM_ASSERT(simplicial_object);return Mass()/Volume();
}
//#####################################################################
// Function Diagonalize_Inertia_Tensor
//#####################################################################
template<class T> void
Diagonalize_Inertia_Tensor_Helper(const SYMMETRIC_MATRIX<T,0>& inertia_input,DIAGONAL_MATRIX<T,0>& inertia,ROTATION<VECTOR<T,1> >& rotation)
{
}
//#####################################################################
// Function Diagonalize_Inertia_Tensor
//#####################################################################
template<class T> void
Diagonalize_Inertia_Tensor_Helper(const SYMMETRIC_MATRIX<T,1>& inertia_input,DIAGONAL_MATRIX<T,1>& inertia,ROTATION<VECTOR<T,2> >& rotation)
{
    inertia=DIAGONAL_MATRIX<T,1>(inertia_input.x00);
}
//#####################################################################
// Function Diagonalize_Inertia_Tensor
//#####################################################################
template<class T> void
Diagonalize_Inertia_Tensor_Helper(const SYMMETRIC_MATRIX<T,3>& inertia_input,DIAGONAL_MATRIX<T,3>& inertia,ROTATION<VECTOR<T,3> >& rotation)
{
    MATRIX<T,3> rotation_matrix;
    inertia_input.Solve_Eigenproblem(inertia,rotation_matrix);
    rotation=ROTATION<VECTOR<T,3> >(rotation_matrix);
}
//#####################################################################
// Function Diagonalize_Inertia_Tensor
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Diagonalize_Inertia_Tensor(const SYMMETRIC_MATRIX<T,TV::SPIN::m>& inertia_tensor_at_center_of_mass)
{
    Diagonalize_Inertia_Tensor_Helper(inertia_tensor_at_center_of_mass,Inertia_Tensor(),Frame().r);
    Update_Angular_Velocity();
}
//#####################################################################
// Function Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface
//#####################################################################
template<class TV> template<class T2> void RIGID_BODY<TV>::
Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface(TETRAHEDRALIZED_VOLUME<T2>& tetrahedralized_volume,TRIANGULATED_SURFACE<T>& triangulated_surface,const T cell_size,
    const int subdivision_loops,const bool (*create_levelset_test)(TETRAHEDRALIZED_VOLUME<T>&),const bool use_implicit_surface_maker,const T shrink_levelset_amount)
{
    Add_Structure(tetrahedralized_volume);
    TRIANGULATED_SURFACE<T>* triangulated_surface_condensed=TRIANGULATED_SURFACE<T>::Create();
    triangulated_surface_condensed->mesh.Initialize_Mesh(triangulated_surface.mesh);
    triangulated_surface_condensed->particles.Initialize(triangulated_surface.particles); 
    triangulated_surface_condensed->Discard_Valence_Zero_Particles_And_Renumber();
    for(int i=0;i<subdivision_loops;i++) triangulated_surface_condensed->Linearly_Subdivide();
    Add_Structure(*triangulated_surface_condensed);
    if(!create_levelset_test || (*create_levelset_test)(tetrahedralized_volume)){
        LEVELSET_IMPLICIT_OBJECT<TV>* levelset=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        levelset->levelset.grid=GRID<TV>::Create_Grid_Given_Cell_Size(*triangulated_surface_condensed->bounding_box,cell_size,false,2); 
        if(use_implicit_surface_maker){
            LEVELSET_MAKER<T> levelset_maker;
            levelset_maker.Compute_Signed_Distance_Function(); 
            levelset_maker.Compute_Level_Set(*triangulated_surface_condensed,levelset->levelset.grid,levelset->levelset.phi);
            levelset->levelset.phi+=shrink_levelset_amount;}
        else SIGNED_DISTANCE::Calculate(*triangulated_surface_condensed,levelset->levelset.grid,levelset->levelset.phi);
        IMPLICIT_OBJECT<TV>* implicit_surface=0;
            implicit_surface=levelset;
        if(implicit_surface) Add_Structure(*implicit_surface);}
}
//#####################################################################
// Function Effective_Inertia_Inverse
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Effective_Inertia_Inverse(MATRIX<T,mm>& extended_mass_inverse,const TV& location) const
{
    if(Has_Infinite_Inertia()){extended_mass_inverse=MATRIX<T,mm>();return;}

    TV r=location-Frame().t;
    MATRIX<T,T_SPIN::m> Ii=World_Space_Inertia_Tensor_Inverse();
    MATRIX<T,T_SPIN::m,TV::m> Ii_rs=Ii.Times_Cross_Product_Matrix(r);
    MATRIX<T,TV::m> mi_rst_Ii_rs=Ii_rs.Cross_Product_Matrix_Transpose_Times(r)+1/Mass();
    extended_mass_inverse.Set_Submatrix(0,0,mi_rst_Ii_rs);
    extended_mass_inverse.Set_Submatrix(0,TV::m,Ii_rs.Transposed());
    extended_mass_inverse.Set_Submatrix(TV::m,0,Ii_rs);
    extended_mass_inverse.Set_Submatrix(TV::m,TV::m,Ii);
}
//#####################################################################
// Function Effective_Inertia_At_Point
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Effective_Inertia(MATRIX<T,mm>& extended_mass_inverse,const TV& location) const
{
    PHYSBAM_ASSERT(!Has_Infinite_Inertia());

    TV r=location-Frame().t;
    T m=Mass();
    MATRIX<T,TV::m> A=MATRIX<T,TV::m>()+m;
    MATRIX<T,T_SPIN::m,TV::m> B=-m*MATRIX<T,T_SPIN::m,TV::m>::Cross_Product_Matrix(r);
    MATRIX<T,T_SPIN::m> C=World_Space_Inertia_Tensor_Inverse()-B.Times_Cross_Product_Matrix_Transpose(r);
    extended_mass_inverse.Set_Submatrix(0,0,A);
    extended_mass_inverse.Set_Submatrix(0,TV::m,B.Transposed());
    extended_mass_inverse.Set_Submatrix(TV::m,0,B);
    extended_mass_inverse.Set_Submatrix(TV::m,TV::m,C);
}
//#####################################################################
// Function Effective_Inertia_Inverse_Times
//#####################################################################
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Effective_Inertia_Inverse_Times(const TWIST<TV>& wrench,const TV& location) const
{
    if(Has_Infinite_Inertia()) return TWIST<TV>();
    return Scatter(Inertia_Inverse_Times(Gather(wrench,location)),location);
}
//#####################################################################
// Function Effective_Inertia_Times
//#####################################################################
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Effective_Inertia_Times(const TWIST<TV>& twist,const TV& location) const
{
    PHYSBAM_ASSERT(!Has_Infinite_Inertia());

    TV r=location-Frame().t;
    T_SPIN torque=twist.angular-TV::Cross_Product(r,twist.linear);
    T_SPIN omega=World_Space_Inertia_Tensor_Times(torque);
    return TWIST<TV>(twist.linear*Mass()-TV::Cross_Product(omega,r),omega);
}
template<class TV> void RIGID_BODY<TV>::
Gather_Matrix(MATRIX<T,mm>& gather,const TV& location) const
{
    gather.Set_Identity_Matrix();
    gather.Set_Submatrix(TV::m,0,MATRIX<T,T_SPIN::m,TV::m>::Cross_Product_Matrix(location-Frame().t));
}
template<class TV> void RIGID_BODY<TV>::
Scatter_Matrix(MATRIX<T,mm>& scatter,const TV& location) const
{
    scatter.Set_Identity_Matrix();
    scatter.Set_Submatrix(0,TV::m,MATRIX<T,T_SPIN::m,TV::m>::Cross_Product_Matrix(location-Frame().t).Transposed());
}
//#####################################################################
template class RIGID_BODY<VECTOR<float,1> >;
template class RIGID_BODY<VECTOR<float,2> >;
template class RIGID_BODY<VECTOR<float,3> >;
template class RIGID_BODY<VECTOR<double,1> >;
template class RIGID_BODY<VECTOR<double,2> >;
template class RIGID_BODY<VECTOR<double,3> >;
}
