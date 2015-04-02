//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_OBJECT_COLLISIONS
//#####################################################################
#include <Tools/Arrays/CONSTANT_ARRAY.h>
#include <Tools/Arrays/PROJECTED_ARRAY.h>
#include <Tools/Utilities/Find_Type.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <Rigids/Collisions/COLLISION_HELPER.h>
#include <Rigids/Collisions/COLLISION_PARTICLE_STATE.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISIONS<TV>::
DEFORMABLE_OBJECT_COLLISIONS(DEFORMABLE_PARTICLES<TV>& particles,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection,ARRAY<STRUCTURE<TV>*>& deformable_object_structures,
    COLLISION_BODY_COLLECTION<TV>& collision_body_list)
    :particles(particles),deformable_body_collection(deformable_body_collection),deformable_object_structures(deformable_object_structures),collision_body_list(collision_body_list),
    collision_tolerance((T)1e-6),use_spatial_partition(true),disable_multiple_levelset_collisions(true),use_protectors(false),maximum_levelset_collision_projection_velocity(FLT_MAX),
    protection_thickness((T)1e-3),ignore_priorities(false),collisions_on(false),collision_tolerances(0),thickness_table(0),friction_table(0),use_structure_collide_collision_body(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISIONS<TV>::
~DEFORMABLE_OBJECT_COLLISIONS()
{}
//#####################################################################
// Function Initialize_Object_Collisions
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Initialize_Object_Collisions(const bool collide_with_interior,const T collision_tolerance_input,
    const bool use_spatial_partition_for_levelset_collisions,const bool disable_multiple_levelset_collisions_input,const T maximum_levelset_collision_projection_velocity_input)
{
    collision_tolerance=collision_tolerance_input;
    use_spatial_partition=use_spatial_partition_for_levelset_collisions;
    disable_multiple_levelset_collisions=disable_multiple_levelset_collisions_input;
    maximum_levelset_collision_projection_velocity=maximum_levelset_collision_projection_velocity_input;
    check_collision=CONSTANT_ARRAY<bool>(particles.Size(),false);
    particle_states.Resize(particles.Size());
    particle_to_collision_body_id.Resize(particles.Size());
    Reset_Object_Collisions(); // in case collisions already exist
    for(int c=0;c<collision_structures.m;c++){
        if(TRIANGULATED_AREA<T>* triangulated_area=dynamic_cast<TRIANGULATED_AREA<T>*>(collision_structures(c)))
            Add_Collision_Mesh(triangulated_area->mesh,collide_with_interior);
        else if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(collision_structures(c)))
            Add_Collision_Mesh(triangulated_surface->mesh,true);
        else if(SEGMENTED_CURVE<TV>* segmented_curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(collision_structures(c))){
            SEGMENT_MESH& mesh=segmented_curve->mesh;
            for(int s=0;s<mesh.elements.m;s++){
                check_collision.Subset(mesh.elements(s)).Fill(true);}}
        else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(collision_structures(c)))
            Add_Collision_Mesh(tetrahedralized_volume->mesh,collide_with_interior);
        else if(HEXAHEDRALIZED_VOLUME<T>* hexahedralized_volume=dynamic_cast<HEXAHEDRALIZED_VOLUME<T>*>(collision_structures(c)))
            Add_Collision_Mesh(hexahedralized_volume->mesh,collide_with_interior);
        else if(EMBEDDING<TV>* embedding=dynamic_cast<EMBEDDING<TV>*>(collision_structures(c)))
            if(collide_with_interior && !dynamic_cast<EMBEDDED_MATERIAL_SURFACE<VECTOR<T,3>,2>*>(embedding)) PHYSBAM_FATAL_ERROR();
            else Add_Collision_Mesh(embedding->material_surface_mesh,true);
        else if(FREE_PARTICLES<TV>* free_particles=dynamic_cast<FREE_PARTICLES<TV>*>(collision_structures(c))){
            check_collision.Subset(free_particles->nodes).Fill(true);}
        else if(BEZIER_SPLINE<TV,3>* spline=dynamic_cast<BEZIER_SPLINE<TV,3>*>(collision_structures(c))){}
        else if(B_SPLINE<TV,3>* spline=dynamic_cast<B_SPLINE<TV,3>*>(collision_structures(c))){}
        else PHYSBAM_NOT_IMPLEMENTED(LOG::sprintf("Collisions with %s",typeid(*collision_structures(c)).name()));}
    check_collision.Subset(ignored_nodes).Fill(false);
}
//#####################################################################
// Function Use_Structure_Skip_Collision_Body
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Use_Structure_Collide_Collision_Body(const bool value)
{
    use_structure_collide_collision_body=value;
    structure_collide_collision_body.Resize(deformable_object_structures.m);
}
//#####################################################################
// Function Compute_Candidate_Nodes_For_Collision_Body_Collisions
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Compute_Candidate_Nodes_For_Collision_Body_Collisions(const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies)
{
    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
    collision_body_candidate_nodes.Resize(collision_body_list.Size(),false);
    for(COLLISION_GEOMETRY_ID id(0);id<collision_body_list.Size();id++) collision_body_candidate_nodes(id).Remove_All();
    if(!bodies.m) return;
    if(use_spatial_partition){
        ARRAY<COLLISION_GEOMETRY_ID> collision_body_indices;
        collision_body_indices.Preallocate(100);
        for(int i=0;i<deformable_body_collection.simulated_particles.m;i++){
            int p=deformable_body_collection.simulated_particles(i);
            if(check_collision(p)){
                if(thickness_table && thickness_table->Contains(p))
                    collision_body_list.spatial_partition->collision_body_thickness=thickness_table->Get(p);
                collision_body_list.spatial_partition->Get_Potential_Collisions(particles.X(p),collision_body_indices);
                for(int j=0;j<collision_body_indices.m;j++)
                    collision_body_candidate_nodes(collision_body_list.bodies(collision_body_indices(j))->collision_geometry_id).Append(p);}}}
    else{
        ARRAY<int> general_collision_body_candidate_nodes;
        for(int i=0;i<deformable_body_collection.simulated_particles.m;i++){
            int p=deformable_body_collection.simulated_particles(i);
            if(check_collision(p))
                general_collision_body_candidate_nodes.Append(p);}
        for(COLLISION_GEOMETRY_ID i(0);i<bodies.m;i++)
            if(bodies(i))
                collision_body_candidate_nodes(i)=general_collision_body_candidate_nodes;}
    if(use_structure_collide_collision_body) for(COLLISION_GEOMETRY_ID body_id(0);body_id<bodies.m;body_id++) if(bodies(body_id)){
        for(int i=collision_body_candidate_nodes(body_id).m-1;i>=0;i--){
            int p=collision_body_candidate_nodes(body_id)(i),structure=particle_to_structure(p);
            if(structure && !structure_collide_collision_body(structure).Contains(body_id))
                collision_body_candidate_nodes(body_id).Remove_Index_Lazy(i);}}
}
//#####################################################################
// Function Adjust_Nodes_For_Collision_Body_Collisions
//#####################################################################
// TODO: should be accelerated by using the hierarchy to eliminate tests on parts of the surface that aren't close to collision bodies.
template<class TV> int DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Nodes_For_Collision_Body_Collisions(BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,ARRAY_VIEW<const TV> X_old,const T dt,
    const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>* bodies)
{
    binding_list.Clamp_Particles_To_Embedded_Positions();
    binding_list.Clamp_Particles_To_Embedded_Velocities();
    soft_bindings.Clamp_Particles_To_Embedded_Positions(true);
    soft_bindings.Clamp_Particles_To_Embedded_Velocities(true); // TODO: move this elsewhere
    if(!bodies) bodies=&collision_body_list.bodies;
    if(!bodies->m) return 0;
    int interactions=0;
    Compute_Candidate_Nodes_For_Collision_Body_Collisions(*bodies);

    if(use_protectors){
        ARRAY<bool> is_protected(check_collision.m);
        ARRAY<TV> X_save(particles.X),V_old(particles.V);
        for(COLLISION_GEOMETRY_ID body_id(0);body_id<bodies->m;body_id++) if((*bodies)(body_id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(*bodies)(body_id);
            for(int j=collision_body_candidate_nodes(body_id).m-1;j>=0;j--){
                int node=collision_body_candidate_nodes(body_id)(j);
                if(protecting_bodies_of_nodes(node).Contains(body_id)){
                    if(!is_protected(node) && collision_body.Implicit_Geometry_Lazy_Inside_Extended_Levelset(X_save(node),(T)protection_thickness)){
                        is_protected(node)=true;
                        particles.X(node)=X_save(node);
                        particles.V(node)=V_old(node);
                        particle_states(node).enforce=false;}}
                else if(is_protected(node)) collision_body_candidate_nodes(body_id).Remove_Index_Lazy(j);}
            interactions+=Adjust_Nodes_For_Collisions(collision_body,particles,collision_body_candidate_nodes(body_id),check_collision,particle_states,
                particle_to_collision_body_id,friction_table,thickness_table);}
        ARRAY<int> collision_count(check_collision.m);
        for(COLLISION_GEOMETRY_ID body_id(0);body_id<bodies->m;body_id++) if((*bodies)(body_id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(*bodies)(body_id);
            for(int j=0;j<collision_body_candidate_nodes(body_id).m;j++){
                int node=collision_body_candidate_nodes(body_id)(j);
                if(particle_states(node).enforce && (!is_protected(node) || protecting_bodies_of_nodes(node).Contains(body_id))
                    && collision_body.Implicit_Geometry_Lazy_Inside(particles.X(node),(thickness_table?thickness_table->Get_Default(node,0):0)-(T)1e-5)){
                    collision_count(node)++;
                    if(collision_count(node)>1){
                        particle_states(node).enforce=false;
                        particles.X(node)=X_save(node);
                        particles.V(node)=V_old(node);}}}}}
    else if(disable_multiple_levelset_collisions){
        ARRAY<TV> X_save(particles.X),V_old(particles.V);
        for(COLLISION_GEOMETRY_ID body_id(0);body_id<bodies->m;body_id++) if((*bodies)(body_id))
            interactions+=Adjust_Nodes_For_Collisions(*(*bodies)(body_id),particles,collision_body_candidate_nodes(body_id),
                check_collision,particle_states,particle_to_collision_body_id,friction_table,thickness_table);
        ARRAY<int> collision_count(check_collision.m);
        for(COLLISION_GEOMETRY_ID body_id(0);body_id<bodies->m;body_id++) if((*bodies)(body_id))
            for(int j=0;j<collision_body_candidate_nodes(body_id).m;j++){int node=collision_body_candidate_nodes(body_id)(j);
                if(particle_states(node).enforce && (*bodies)(body_id)->Implicit_Geometry_Lazy_Inside(particles.X(node),(thickness_table?thickness_table->Get_Default(node,0):0)-(T)1e-5)){
                    collision_count(node)++;
                    if(collision_count(node)>1){
                        particle_states(node).enforce=false;
                        particles.X(node)=X_save(node);
                        particles.V(node)=V_old(node);}}}}
    else for(COLLISION_GEOMETRY_ID body_id=bodies->m-1;body_id>=COLLISION_GEOMETRY_ID(0);body_id--)
        interactions+=Adjust_Nodes_For_Collisions(*(*bodies)(body_id),particles,collision_body_candidate_nodes(body_id),
            check_collision,particle_states,particle_to_collision_body_id,friction_table,thickness_table);

    binding_list.Clamp_Particles_To_Embedded_Positions();
    binding_list.Clamp_Particles_To_Embedded_Velocities();
    return interactions;
}
//#####################################################################
// Function Adjust_Existing_Nodes_For_Collision_Body_Collisions
//#####################################################################
template<class TV> int DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Existing_Nodes_For_Collision_Body_Collisions(const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>* bodies)
{
    int interactions=0;
    for(int i=0;i<deformable_body_collection.simulated_particles.m;i++){
        int p=deformable_body_collection.simulated_particles(i);
        COLLISION_PARTICLE_STATE<TV>& collision=particle_states(p);
        if(!collision.enforce) continue;
        interactions++;
        Adjust_Point_For_Collision(dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>&>(collision_body_list(particle_to_collision_body_id(p))),p,collision,collision.friction);}
    return interactions;
}
//#####################################################################
// Function Set_Collision_Velocities
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Set_Collision_Velocities(ARRAY_VIEW<TV> V) // for external forces and velocities
{
    enforced_particles.Remove_All();
    for(int i=0;i<deformable_body_collection.simulated_particles.m;i++){
        int p=deformable_body_collection.simulated_particles(i);
        COLLISION_PARTICLE_STATE<TV>& collision=particle_states(p);
        if(collision.enforce){
            T VN=TV::Dot_Product(V(p),collision.normal);
            V(p)+=(collision.VN-VN)*collision.normal;}}
// TODO: put this back for efficiency and test
//            enforced_particles.Append(PAIR<int,COLLISION_PARTICLE_STATE<TV> >(p,collision));
}
//#####################################################################
// Function Zero_Out_Collision_Velocities
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Zero_Out_Collision_Velocities(ARRAY_VIEW<TV> V) // for external forces and velocities
{
    for(int p=0;p<particle_states.m;p++){COLLISION_PARTICLE_STATE<TV>& collision=particle_states(p);
        if(collision.enforce){
            V(p)-=TV::Dot_Product(V(p),collision.normal)*collision.normal;}}
// TODO: put this back for efficiency and test
//    for(int i=0;i<enforced_particles.m;i++){
//        const int p=enforced_particles(i).x;const TV& collision_normal=enforced_particles(i).y.normal;
//        V(p)-=TV::Dot_Product(V(p),collision_normal)*collision_normal;}
}
//#####################################################################
// Function Reset_Object_Collisions
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Reset_Object_Collisions() // The unusual signature is not great
{
    collisions_on=false;
    if(particle_states.m)
        particle_states.Subset(deformable_body_collection.simulated_particles).
            template Project<bool,&COLLISION_PARTICLE_STATE<TV>::enforce>().Fill(false);
}
//#####################################################################
// Function Add_Collision_Mesh
//#####################################################################
template<class TV> template<class T_MESH> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Add_Collision_Mesh(T_MESH& mesh,const bool collide_with_interior)
{
    if(collide_with_interior) check_collision.Subset(mesh.elements.Flattened()).Fill(true);
    else{bool boundary_nodes_defined=mesh.boundary_nodes!=0;
        if(!boundary_nodes_defined) mesh.Initialize_Boundary_Nodes();
        check_collision.Subset(*mesh.boundary_nodes).Fill(true);
        if(!boundary_nodes_defined){delete mesh.boundary_nodes;mesh.boundary_nodes=0;}}
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Update_Simulated_Particles()
{
    particle_to_structure.Resize(particles.Size(),false,false);
    particle_to_structure.Fill(0);
    for(int s=0;s<deformable_object_structures.m;s++) deformable_object_structures(s)->Mark_Nodes_Referenced(particle_to_structure,s);
}
//#####################################################################
// Function Adjust_Nodes_For_Collisions_Helper
//#####################################################################
template<class TV,class T> int
Adjust_Nodes_For_Collisions_Helper(COLLISION_GEOMETRY<TV>& body,SOFT_BINDINGS<TV>& soft_bindings,DEFORMABLE_PARTICLES<TV>& collision_particles,const ARRAY<int>& nodes_to_check,
    const ARRAY<bool>& particle_on_surface,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const HASHTABLE<int,T> *friction_table,const HASHTABLE<int,T> *thickness_table)
{ 
    return 0;
}
//#####################################################################
// Function Adjust_Nodes_For_Collisions_Helper
//#####################################################################
template<class T> int
Adjust_Nodes_For_Collisions_Helper(COLLISION_GEOMETRY<VECTOR<T,3> >& body,SOFT_BINDINGS<VECTOR<T,3> >& soft_bindings,DEFORMABLE_PARTICLES<VECTOR<T,3> >& collision_particles,const ARRAY<int>& nodes_to_check,
    const ARRAY<bool>& particle_on_surface,ARRAY<COLLISION_PARTICLE_STATE<VECTOR<T,3> > >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const HASHTABLE<int,T> *friction_table,const HASHTABLE<int,T> *thickness_table)
{
    if(TETRAHEDRON_COLLISION_BODY<T>* tetrahedron_collision_body=dynamic_cast<TETRAHEDRON_COLLISION_BODY<T>*>(&body))
        return tetrahedron_collision_body->Adjust_Nodes_For_Collisions(collision_particles,soft_bindings,nodes_to_check,particle_on_surface,
            collision_particle_state,particle_to_collision_geometry_id,friction_table,thickness_table);
    return 0;
}
//#####################################################################
// Function Adjust_Nodes_For_Collisions
//#####################################################################
template<class TV> int DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Nodes_For_Collisions(COLLISION_GEOMETRY<TV>& body, DEFORMABLE_PARTICLES<TV>& collision_particles,const ARRAY<int>& nodes_to_check,
    const ARRAY<bool>& particle_on_surface,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const HASHTABLE<int,T> *friction_table,const HASHTABLE<int,T> *thickness_table)
{
    if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_body=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(&body))
        return Adjust_Nodes_For_Collisions(*rigid_collision_body,collision_particles,nodes_to_check,
            collision_particle_state,particle_to_collision_geometry_id,friction_table,thickness_table);
    return Adjust_Nodes_For_Collisions_Helper(body,deformable_body_collection.soft_bindings,collision_particles,nodes_to_check,particle_on_surface,collision_particle_state,
        particle_to_collision_geometry_id,friction_table,thickness_table);
}
//#####################################################################
// Function Adjust_Nodes_For_Collisions
//#####################################################################
template<class TV> int DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Nodes_For_Collisions(RIGID_COLLISION_GEOMETRY<TV>& body,DEFORMABLE_PARTICLES<TV>& collision_particles,
    const ARRAY<int>& nodes_to_check,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const HASHTABLE<int,T> *friction_table,const HASHTABLE<int,T> *thickness_table)
{
    // TODO: Add some sort of push out
    int interactions=0;

    T depth; 
    for(int pp=0;pp<nodes_to_check.m;pp++){
        int p=nodes_to_check(pp);
        T thickness=thickness_table?thickness_table->Get_Default(p,0):0;
        COLLISION_PARTICLE_STATE<TV>& collision=collision_particle_state(p);
        if(body.Implicit_Geometry_Lazy_Inside_And_Value(deformable_body_collection.binding_list.X(p),depth,thickness)){
            collision.enforce=true;
            interactions++;
            particle_to_collision_geometry_id(p)=body.collision_geometry_id;
            Adjust_Point_For_Collision(body,p,collision,friction_table?friction_table->Get_Default(p,-1):-1);}}
    return interactions;
}
//#####################################################################
// Function Adjust_Point_For_Collision
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Point_For_Collision(RIGID_COLLISION_GEOMETRY<TV>& body,int p,COLLISION_PARTICLE_STATE<TV>& collision,T local_coefficient_of_friction)
{
    TV X=deformable_body_collection.binding_list.X(p),normal=body.Implicit_Geometry_Normal(X);
    TV V_rel=deformable_body_collection.binding_list.V(p)-body.Pointwise_Object_Velocity(X);
    if(local_coefficient_of_friction<0) local_coefficient_of_friction=body.rigid_body.coefficient_of_friction;
    TV impulse=PhysBAM::Compute_Collision_Impulse(normal,deformable_body_collection.binding_list.Impulse_Factor(p),V_rel,(T)0,local_coefficient_of_friction,0);
    deformable_body_collection.binding_list.Apply_Impulse(p,impulse);
    // set collision state
    collision.normal=body.Implicit_Geometry_Normal(X);
    collision.VN=TV::Dot_Product(body.Pointwise_Object_Velocity(X)-deformable_body_collection.binding_list.V(p),collision.normal);
    collision.friction=local_coefficient_of_friction;
}
//#####################################################################
// Function Adjust_Nodes_For_Push_Out_Helper
//#####################################################################
template<class TV,class T> void
Adjust_Nodes_For_Push_Out_Helper(COLLISION_GEOMETRY<TV>& body,SOFT_BINDINGS<TV>& soft_bindings,DEFORMABLE_PARTICLES<TV>& collision_particles,const ARRAY<int>& nodes_to_check,
    const ARRAY<bool>& particle_on_surface,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const HASHTABLE<int,T> *thickness_table)
{
}
//#####################################################################
// Function Adjust_Nodes_For_Push_Out_Helper
//#####################################################################
template<class T> void
Adjust_Nodes_For_Push_Out_Helper(COLLISION_GEOMETRY<VECTOR<T,3> >& body,SOFT_BINDINGS<VECTOR<T,3> >& soft_bindings,DEFORMABLE_PARTICLES<VECTOR<T,3> >& collision_particles,const ARRAY<int>& nodes_to_check,
    const ARRAY<bool>& particle_on_surface,ARRAY<COLLISION_PARTICLE_STATE<VECTOR<T,3> > >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const HASHTABLE<int,T> *thickness_table)
{
    if(TETRAHEDRON_COLLISION_BODY<T>* tetrahedron_collision_body=dynamic_cast<TETRAHEDRON_COLLISION_BODY<T>*>(&body))
        tetrahedron_collision_body->Adjust_Nodes_For_Push_Out(collision_particles,soft_bindings,nodes_to_check,particle_on_surface,
            collision_particle_state,particle_to_collision_geometry_id,thickness_table);
}
//#####################################################################
// Function Adjust_Nodes_For_Push_Out
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Nodes_For_Push_Out(COLLISION_GEOMETRY<TV>& body, DEFORMABLE_PARTICLES<TV>& collision_particles,const ARRAY<int>& nodes_to_check,
    const ARRAY<bool>& particle_on_surface,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const HASHTABLE<int,T> *thickness_table)
{
    if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_body=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(&body))
        return Adjust_Nodes_For_Push_Out(*rigid_collision_body,collision_particles,nodes_to_check,
            collision_particle_state,particle_to_collision_geometry_id,thickness_table);
    return Adjust_Nodes_For_Push_Out_Helper(body,deformable_body_collection.soft_bindings,collision_particles,nodes_to_check,particle_on_surface,collision_particle_state,
        particle_to_collision_geometry_id,thickness_table);
}
//#####################################################################
// Function Adjust_Nodes_For_Push_Out
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Nodes_For_Push_Out(RIGID_COLLISION_GEOMETRY<TV>& body,DEFORMABLE_PARTICLES<TV>& collision_particles,
    const ARRAY<int>& nodes_to_check,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const HASHTABLE<int,T> *thickness_table)
{
    T depth; 
    for(int pp=0;pp<nodes_to_check.m;pp++){
        int p=nodes_to_check(pp);
        T thickness=thickness_table?thickness_table->Get_Default(p,0):0;
        if(body.Implicit_Geometry_Lazy_Inside_And_Value(deformable_body_collection.binding_list.X(p),depth,thickness)){
            Adjust_Point_For_Push_Out(body,p,depth);}}
}
//#####################################################################
// Function Adjust_Point_For_Collision
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Point_For_Push_Out(RIGID_COLLISION_GEOMETRY<TV>& body,int p,T depth)
{
    TV X=deformable_body_collection.binding_list.X(p),normal=body.Implicit_Geometry_Normal(X),dX=normal*depth;
    if(depth>=0) return;
    TV impulse=-depth/deformable_body_collection.binding_list.One_Over_Effective_Mass(p,normal)*normal;
    deformable_body_collection.binding_list.Apply_Push(p,impulse);
}
//#####################################################################
// Function Adjust_Point_For_Collision
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Process_Push_Out()
{
    Compute_Candidate_Nodes_For_Collision_Body_Collisions(collision_body_list.bodies);
    for(COLLISION_GEOMETRY_ID body_id=collision_body_list.bodies.m-1;body_id>=COLLISION_GEOMETRY_ID(0);body_id--)
        Adjust_Nodes_For_Push_Out(*collision_body_list.bodies(body_id),particles,collision_body_candidate_nodes(body_id),
            check_collision,particle_states,particle_to_collision_body_id,thickness_table);
}
//#####################################################################
namespace PhysBAM{
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<float,1> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<float,2> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<float,3> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<double,1> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<double,2> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<double,3> >;
}
