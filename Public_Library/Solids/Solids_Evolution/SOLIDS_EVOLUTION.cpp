//#####################################################################
// Copyright 2006-2007, Ron Fedkiw, Eran Guendalman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EVOLUTION
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <Rigids/Muscles/MUSCLE_LIST.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <limits>
using namespace PhysBAM;
template<class TV> SOLIDS_EVOLUTION_CALLBACKS<TV> SOLIDS_EVOLUTION<TV>::solids_evolution_callbacks_default;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_EVOLUTION<TV>::
SOLIDS_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities)
    :solid_body_collection(solid_body_collection_input),solids_parameters(solids_parameters_input),rigid_body_collisions(0),rigid_deformable_collisions(0),time(0),
    solids_evolution_callbacks(&solids_evolution_callbacks_default),fully_implicit(false),
    GV_B(solid_body_collection.New_Generalized_Velocity()),
    kinematic_evolution(solid_body_collection.rigid_body_collection,example_forces_and_velocities),
    example_forces_and_velocities(example_forces_and_velocities)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_EVOLUTION<TV>::
~SOLIDS_EVOLUTION()
{
    delete rigid_body_collisions;
    delete rigid_deformable_collisions;
    krylov_vectors.Delete_Pointers_And_Clean_Memory();
    delete &GV_B;
}
//#####################################################################
// Function Euler_Step_Position
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Euler_Step_Position(const T dt,const T time)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    deformable_body_collection.particles.Euler_Step_Position(deformable_body_collection.dynamic_particles,dt);
    for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++)
        Euler_Step_Position(dt,time,solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i)); // TODO: avoid unnecessary Update_Angular_Velocity
    Set_External_Positions(deformable_body_collection.particles.X,time+dt);
    kinematic_evolution.Set_External_Positions(solid_body_collection.rigid_body_collection.rigid_body_particles.frame,time+dt);
    solid_body_collection.rigid_body_collection.Update_Angular_Velocity(); // Note: Possibly remove as we restore velocities right after this function.
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
     if(mpi_solids){
         mpi_solids->Exchange_Force_Boundary_Data_Global(deformable_body_collection.particles.X);
         mpi_solids->Exchange_Binding_Boundary_Data_Global(deformable_body_collection.particles.X);
/*      TODO: exchange rigid body particles data
        mpi_solids->Exchange_Boundary_Data_Global(solid_body_collection.rigid_body_collection.rigid_body_particles.frame);
        TODO: update angular velocity for received boundary data*/
}
}
//#####################################################################
// Function Save_Position
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Save_Position(ARRAY<TV>& X,ARRAY<FRAME<TV> >& rigid_frame)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    const ARRAY<int>& simulated_particles=solid_body_collection.deformable_body_collection.simulated_particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    const ARRAY<int>& simulated_rigid_body_particles=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles;
    X.Resize(particles.Size(),no_init);
    X.Subset(simulated_particles)=particles.X.Subset(simulated_particles);
    rigid_frame.Resize(rigid_body_collection.rigid_body_particles.Size(),no_init);
    rigid_frame.Subset(simulated_rigid_body_particles)=rigid_body_collection.rigid_body_particles.frame.Subset(simulated_rigid_body_particles);
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i)){
        if(!rigid_body_collection.Rigid_Body(i).Is_Simulated())rigid_frame(i)=rigid_body_collection.rigid_body_particles.frame(i);}
}
//#####################################################################
// Function Restore_Position
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Restore_Position(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const FRAME<TV> > rigid_frame)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    const ARRAY<int>& simulated_particles=solid_body_collection.deformable_body_collection.simulated_particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    const ARRAY<int>& simulated_rigid_body_particles=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles;
    PHYSBAM_ASSERT(X.Size()==particles.Size());PHYSBAM_ASSERT(rigid_frame.Size()==rigid_body_collection.rigid_body_particles.Size());
    particles.X.Subset(simulated_particles)=X.Subset(simulated_particles);
    rigid_body_collection.rigid_body_particles.frame.Subset(simulated_rigid_body_particles)=rigid_frame.Subset(simulated_rigid_body_particles);
    for(int i=0;i<simulated_rigid_body_particles.m;i++) rigid_body_collection.Rigid_Body(simulated_rigid_body_particles(i)).Update_Angular_Velocity();
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){rigid_body_collection.rigid_body_particles.frame(i)=rigid_frame(i);body.Update_Angular_Velocity();}}
}
//#####################################################################
// Function Restore_Position
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Restore_Position_After_Hypothetical_Position_Evolution(ARRAY<TV>& X_save,ARRAY<FRAME<TV> >& rigid_frame_save)
{
    Restore_Position(X_save,rigid_frame_save);
}
//#####################################################################
// Function Initialize_Deformable_Objects
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Initialize_Deformable_Objects(const T frame_rate,const bool restart)
{
    if(!restart){
        solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
        solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
        solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
        solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();}

    solid_body_collection.Update_Simulated_Particles();

    solid_body_collection.deformable_body_collection.collisions.Initialize_Object_Collisions(solids_parameters.deformable_object_collision_parameters.collide_with_interior,
        solids_parameters.deformable_object_collision_parameters.collision_tolerance,
        solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects,
        solids_parameters.deformable_object_collision_parameters.disable_multiple_levelset_collisions,
        solids_parameters.deformable_object_collision_parameters.maximum_levelset_collision_projection_velocity);
    solid_body_collection.collision_body_list.Update_Spatial_Partition(solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_heuristic,
        solids_parameters.deformable_object_collision_parameters.spatial_partition_number_of_cells,
        solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_scale_factor);

    solid_body_collection.Set_CFL_Number(solids_parameters.cfl);
    example_forces_and_velocities.Update_Time_Varying_Material_Properties(time);
    solid_body_collection.Update_Position_Based_State(time,true,true);
}
//#####################################################################
// Function Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions
//#####################################################################
template<class TV> bool SOLIDS_EVOLUTION<TV>::
Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(const T dt,const T time,int& repulsions_found,int& collisions_found,const bool exit_early)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry=solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry;
    ARRAY<bool>& modified=geometry.modified_full;ARRAY<TV>& X_self_collision_free=geometry.X_self_collision_free;

    repulsions_found=0;collisions_found=0; // important to initialize these in case repulsions/collisions are turned off
    if(!solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions && solids_parameters.triangle_collision_parameters.turn_off_all_collisions){
        LOG::cout<<"all repulsions and collisions are turned off - nothing to do"<<std::endl;return false;}

    // update soft bound particle positions and velocities
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true);
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        LOG::SCOPE scope("Collisions All Gather","Collisions All Gather");
        solid_body_collection.deformable_body_collection.mpi_solids->All_Gather_Particles(particles.X,particles.V);}

    modified=CONSTANT_ARRAY<bool>(particles.Size(),false);

    // compute average velocity, but save final velocities in case of no repulsion or collisions
    ARRAY<TV> V_save(particles.V);

    ARRAY<TV> V_averaged(1/dt*(particles.X-X_self_collision_free));
    particles.V=V_averaged;

    // check for triangle pairs that are already intersecting, and label them to be ignored
    if(geometry.allow_intersections) geometry.Compute_Intersecting_Segment_Face_Pairs();

    // self repulsion - adjust time n+1/2 velocity and time n+1 position
    if(solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions){
        if(solid_body_collection.deformable_body_collection.mpi_solids) PHYSBAM_NOT_IMPLEMENTED("per collision step repulsions under mpi");
        repulsions_found=solid_body_collection.deformable_body_collection.triangle_repulsions.Adjust_Velocity_For_Self_Repulsion(dt,false);}

    // self collisions - adjust time n+1/2 velocity and time n+1 position
    if(!solids_parameters.triangle_collision_parameters.turn_off_all_collisions)
        collisions_found=solid_body_collection.deformable_body_collection.triangle_collisions.Adjust_Velocity_For_Self_Collisions(dt,time,exit_early);

    // update velocities
    if(collisions_found<0) // failed to resolve collisions (should only happen when exit_early=true) - otherwise any collisions that were found have been resolved
        PHYSBAM_NOT_IMPLEMENTED("exit_early=true");
    else if(collisions_found>0){ // restore unmodified velocities
        for(int p=0;p<particles.Size();p++) particles.V(p)=modified(p)?V_save(p)+particles.V(p)-V_averaged(p):V_save(p);}
    else if(repulsions_found){ // repulsions only, restore velocities for unmodified and apply velocity delta otherwise
        for(int p=0;p<particles.Size();p++) particles.V(p)=modified(p)?V_save(p)+particles.V(p)-V_averaged(p):V_save(p);}
    else{particles.V=V_save;return false;} // restore all the unmodified velocities
    example_forces_and_velocities.Update_Time_Varying_Material_Properties(time);
    solid_body_collection.Update_Position_Based_State(time,false,true);

    return true;
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Postprocess_Frame(const int frame)
{
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    if(articulated_rigid_body.use_muscle_actuators) for(int i=0;i<articulated_rigid_body.muscle_list->muscles.m;i++){
        articulated_rigid_body.muscle_list->muscles(i)->Set_Segment_Activations(articulated_rigid_body.muscle_activations(i));
        articulated_rigid_body.muscle_list->muscles(i)->Update_Segments();}
    if(solids_parameters.triangle_collision_parameters.perform_self_collision && solids_parameters.triangle_collision_parameters.check_mesh_for_self_intersection
        && (!solid_body_collection.deformable_body_collection.mpi_solids || solid_body_collection.deformable_body_collection.mpi_solids->rank==0)){
        if(solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Check_For_Intersection(false,solids_parameters.triangle_collision_parameters.collisions_small_number)){
            PHYSBAM_DEBUG_WRITE_SUBSTEP("Intersections Found",0);
            throw std::runtime_error("Intersections Found");}}
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    example_forces_and_velocities.Set_External_Positions(X,time);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    if(!solids_parameters.use_rigid_deformable_contact && solid_body_collection.deformable_body_collection.collisions.collisions_on) solid_body_collection.deformable_body_collection.collisions.Set_Collision_Velocities(V);
    example_forces_and_velocities.Set_External_Velocities(V,velocity_time,current_position_time);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    if(!solids_parameters.use_rigid_deformable_contact && solid_body_collection.deformable_body_collection.collisions.collisions_on) solid_body_collection.deformable_body_collection.collisions.Zero_Out_Collision_Velocities(V);
    example_forces_and_velocities.Zero_Out_Enslaved_Velocity_Nodes(V,velocity_time,current_position_time);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
    twist.Subset(solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies).Fill(TWIST<TV>());
    example_forces_and_velocities.Zero_Out_Enslaved_Velocity_Nodes(twist,velocity_time,current_position_time);
}
template<class T>
inline int Correct_Orientation_For_Kinetic_Energy_Using_Direction(RIGID_BODY<VECTOR<T,3> >& rigid_body,const T KE0,const VECTOR<T,3>& direction,const bool use_extrema)
{
    typedef VECTOR<T,3> TV;ROTATION<TV>& R=rigid_body.Frame().r;const TV &w=rigid_body.Twist().angular,&L=rigid_body.Angular_Momentum();
    T KE=(T).5*TV::Dot_Product(w,L),k=KE-KE0,error=abs(k);
    if(error<=4*std::numeric_limits<T>::epsilon()*KE0) return 0;
    if(KE0<(T)sqr(std::numeric_limits<T>::epsilon())) return -1;
    TV s=direction.Normalized();T a=TV::Dot_Product(s,TV::Cross_Product(L,w));
    if(a<0){a=-a;s=-s;}
    if(a<16*std::numeric_limits<T>::epsilon()*KE) return 1;
    TV s_cross_L=TV::Cross_Product(s,L);
    T b=(T).5*TV::Dot_Product(rigid_body.World_Space_Inertia_Tensor_Inverse_Times(s_cross_L),s_cross_L);
    T c=KE-b,r=sqrt(sqr(a)+sqr(c)),n=c-2*k,q2=sqr(a)+4*k*(c-k);
    if(q2>=0){
        T q=a*sqrt(q2),p=sqr(a)+q2+2*q;
        T x=(T).5*sqrt(p+sqr(c+n))/r,y=2*fabs(k)/sqrt(p+sqr(c-n));
        if(k-c*sqr(y)<0) y=-y;
        R=ROTATION<TV>::From_Rotation_Vector(atan2(y,x)*s)*R;R.Normalize();
        rigid_body.Update_Angular_Velocity();
        return -2;}
    if(use_extrema){
        if(error<=16*std::numeric_limits<T>::epsilon()*KE0) return 2;
        T sk=sign(k),e=sk*(sk*c>0?atan2(r+sk*c,a):atan2(a,r-sk*c));
        R=ROTATION<TV>::From_Rotation_Vector(e*s)*R;R.Normalize();
        rigid_body.Update_Angular_Velocity();
        return 3;}
    return 4;
}
//#####################################################################
// Function Correct_Orientation_For_Kinetic_Energy
//#####################################################################
template<class T>
inline bool Correct_Orientation_For_Kinetic_Energy(RIGID_BODY<VECTOR<T,3> >& rigid_body,const T KE0,int iteration=1)
{
    const int max_orientation_correction_iterations=40;
    if(iteration>max_orientation_correction_iterations) PHYSBAM_FATAL_ERROR("Exceeded maximum number of iterations during solids evolution energy clamping");
    typedef VECTOR<T,3> TV;const TV &w=rigid_body.Twist().angular,&L=rigid_body.Angular_Momentum();
    int status=Correct_Orientation_For_Kinetic_Energy_Using_Direction(rigid_body,KE0,TV::Cross_Product(L,w),false);
    if(status<=0) return true;
    int component;
    if((T).5*TV::Dot_Product(w,L)>KE0) component=rigid_body.Inertia_Tensor().To_Vector().Arg_Max();
    else component=rigid_body.Inertia_Tensor().To_Vector().Arg_Min();
    TV axis=rigid_body.Frame().r.Rotate(TV::Axis_Vector(component));
    status=Correct_Orientation_For_Kinetic_Energy_Using_Direction(rigid_body,KE0,TV::Cross_Product(axis,L),true);
    if(status==3) return Correct_Orientation_For_Kinetic_Energy(rigid_body,KE0,iteration+1);
    return status<=0;
}
//#####################################################################
// Function Update_Rotation_Helper
//#####################################################################
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_BODY<VECTOR<T,3> >& rigid_body,bool correct_evolution_energy)
{
    typedef VECTOR<T,3> TV;ROTATION<TV>& R=rigid_body.Frame().r;const TV &w=rigid_body.Twist().angular,&L=rigid_body.Angular_Momentum();
    T KE0=(T).5*TV::Dot_Product(w,L);
    TV rotate_amount=w-(T).5*dt*rigid_body.World_Space_Inertia_Tensor_Inverse_Times(TV::Cross_Product(w,L));
    R=ROTATION<TV>::From_Rotation_Vector(dt*rotate_amount)*R;R.Normalize();
    rigid_body.Update_Angular_Velocity(); // Note that the value of w changes here.
    if(correct_evolution_energy) Correct_Orientation_For_Kinetic_Energy(rigid_body,KE0);
}
//#####################################################################
// Function Update_Rotation_Helper
//#####################################################################
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_BODY<VECTOR<T,2> >& rigid_body,bool correct_evolution_energy)
{
    rigid_body.Frame().r=ROTATION<VECTOR<T,2> >::From_Rotation_Vector(dt*rigid_body.Twist().angular)*rigid_body.Frame().r;rigid_body.Frame().r.Normalize();
}
//#####################################################################
// Function Update_Rotation_Helper
//#####################################################################
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_BODY<VECTOR<T,1> >& rigid_body,bool correct_evolution_energy)
{}
//#####################################################################
// Euler_Step_Position
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Euler_Step_Position(const T dt,const T time,const int p)
{
    RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(p);
    if(rigid_body.is_static) return;
    else if(solid_body_collection.rigid_body_collection.rigid_body_particles.kinematic(p)){
        kinematic_evolution.Set_External_Positions(rigid_body.Frame(),time+dt,p);}
    else{
        rigid_body.Frame().t+=dt*rigid_body.Twist().linear;
        Update_Rotation_Helper(dt,rigid_body,solids_parameters.rigid_body_evolution_parameters.correct_evolution_energy);}
}
//#####################################################################
// Function Clamp_Velocities
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Clamp_Velocities()
{
    const ARRAY<int>& dynamic_rigid_body_particles=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    T max_linear_velocity_squared=sqr(solids_parameters.rigid_body_evolution_parameters.max_rigid_body_linear_velocity),max_angular_velocity_squared=sqr(solids_parameters.rigid_body_evolution_parameters.max_rigid_body_angular_velocity);
    for(int i=0;i<dynamic_rigid_body_particles.m;i++){int p=dynamic_rigid_body_particles(i);
        T magnitude_squared=rigid_body_particles.twist(p).linear.Magnitude_Squared();
        if(magnitude_squared>max_linear_velocity_squared) rigid_body_particles.twist(p).linear*=solids_parameters.rigid_body_evolution_parameters.max_rigid_body_linear_velocity/sqrt(magnitude_squared);
        magnitude_squared=rigid_body_particles.twist(p).angular.Magnitude_Squared();
        if(magnitude_squared>max_angular_velocity_squared){
            rigid_body_particles.twist(p).angular*=solids_parameters.rigid_body_evolution_parameters.max_rigid_body_angular_velocity/sqrt(magnitude_squared);
            solid_body_collection.rigid_body_collection.Rigid_Body(p).Update_Angular_Momentum();}}
}
//#####################################################################
// Function Initialize_World_Space_Masses
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION<TV>::
Initialize_World_Space_Masses()
{
    world_space_rigid_mass.Resize(solid_body_collection.rigid_body_collection.rigid_body_particles.Size(),no_init);
    world_space_rigid_mass_inverse.Resize(solid_body_collection.rigid_body_collection.rigid_body_particles.Size(),no_init);
    for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
        world_space_rigid_mass(p)=solid_body_collection.rigid_body_collection.State(p).World_Space_Rigid_Mass(RIGID_BODY_MASS<TV>(solid_body_collection.rigid_body_collection.rigid_body_particles.mass(p),solid_body_collection.rigid_body_collection.rigid_body_particles.inertia_tensor(p)));
        world_space_rigid_mass_inverse(p)=solid_body_collection.rigid_body_collection.State(p).World_Space_Rigid_Mass_Inverse(RIGID_BODY_MASS<TV>(solid_body_collection.rigid_body_collection.rigid_body_particles.mass(p),solid_body_collection.rigid_body_collection.rigid_body_particles.inertia_tensor(p)));}
}
//#####################################################################
namespace PhysBAM{
template class SOLIDS_EVOLUTION<VECTOR<float,1> >;
template class SOLIDS_EVOLUTION<VECTOR<float,2> >;
template class SOLIDS_EVOLUTION<VECTOR<float,3> >;
template class SOLIDS_EVOLUTION<VECTOR<double,1> >;
template class SOLIDS_EVOLUTION<VECTOR<double,2> >;
template class SOLIDS_EVOLUTION<VECTOR<double,3> >;
}
