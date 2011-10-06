#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include "DATA_EXCHANGE_CONVERSION.h"
#include "DEFORMABLE_EXAMPLE.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> DEFORMABLE_EXAMPLE<T>::
DEFORMABLE_EXAMPLE(const STREAM_TYPE stream_type)
    :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),fully_implicit(false),scene_file("scene"),read_text(false),write_text(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> DEFORMABLE_EXAMPLE<T>::
~DEFORMABLE_EXAMPLE()
{
    simulation_object_data.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Option_Argument("-fully_implicit","use fully implicit forces");
    parse_args->Add_String_Argument("-scene","file specifying simulation scene");
    parse_args->Add_Option_Argument("-read_text","scene written in text mode");
    parse_args->Add_Option_Argument("-write_text","scene written in text mode");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    fully_implicit=parse_args->Is_Value_Set("-fully_implicit");
    scene_file=parse_args->Get_String_Value("-scene");
    read_text=parse_args->Is_Value_Set("-read_text");
    write_text=parse_args->Is_Value_Set("-write_text");
 }
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    fluids_parameters.simulate=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness=(T)1e-3;
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    solids_parameters.verbose_dt=true;
    solids_parameters.triangle_collision_parameters.total_collision_loops=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=20;
    solids_parameters.enforce_repulsions_in_cg=false;
    solids_parameters.enforce_poststabilization_in_cg=true;
    solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=30;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness*=(T)5;
    solid_body_collection.collision_body_list.Set_Collision_Body_Thickness((T).1);
    solids_parameters.triangle_collision_parameters.collisions_collision_thickness=(T)1e-3;
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
//    solids_parameters.collisions_small_number=(T)1e-5;
//    solids_parameters.self_collision_friction_coefficient=(T)0;

    data_exchange::deformable_body_simulation dbs;
    std::ifstream ifs(scene_file.c_str());
    if(read_text){
        boost::archive::text_iarchive ia(ifs);
        ia >> dbs;}
    else{
        boost::archive::binary_iarchive ia(ifs);
        ia >> dbs;}

    for(size_t i=0; i<dbs.simulation_objects.size(); i++){
        simulation_object_data.Append(new SIMULATION_OBJECT_DATA);
        if(data_exchange::deformable_body* body=dynamic_cast<data_exchange::deformable_body*>(dbs.simulation_objects[i]))
            Add_Deformable_Body(*body,i);
        else if(data_exchange::ground_plane* body=dynamic_cast<data_exchange::ground_plane*>(dbs.simulation_objects[i]))
            Add_Rigid_Body(*body,i);
        else if(data_exchange::scripted_geometry* body=dynamic_cast<data_exchange::scripted_geometry*>(dbs.simulation_objects[i]))
            Add_Rigid_Body(*body,i);}

    for(int i=1;i<=soft_bindings.use_impulses_for_collisions.m;i++)
        if(soft_bindings.use_impulses_for_collisions(i))
            soft_bindings.bindings_using_impulses_for_collisions.Append(i);

    binding_list.Clamp_Particles_To_Embedded_Positions();
    binding_list.Clamp_Particles_To_Embedded_Velocities();
    soft_bindings.Clamp_Particles_To_Embedded_Positions(true);
    soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    for(size_t i=0;i<dbs.simulation_forces.size();i++){
        if(data_exchange::gravity_force* force=dynamic_cast<data_exchange::gravity_force*>(dbs.simulation_forces[i]))
            Add_Force(*force,i);
        else if(data_exchange::volumetric_force* force=dynamic_cast<data_exchange::volumetric_force*>(dbs.simulation_forces[i]))
            Add_Force(*force,i);}

    last_frame=dbs.number_frames;
    this->frame_rate=1/dbs.dt;
    output_directory=dbs.output_directory;
    filename_base=dbs.filename_base;

    // add structures and rigid bodies to collisions

    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++){
        if(dynamic_cast<FREE_PARTICLES<TV>*>(deformable_body_collection.deformable_geometry.structures(i))){
            deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures(i));}
        if(solids_parameters.triangle_collision_parameters.perform_self_collision && dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.deformable_geometry.structures(i)))
            solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.deformable_geometry.structures(i));}

    // disable strain rate CFL for all forces
    for(int i=1;i<=solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->limit_time_step_by_strain_rate=false;
    for(int i=1;i<=deformable_body_collection.deformables_forces.m;i++) deformable_body_collection.deformables_forces(i)->limit_time_step_by_strain_rate=false;
    for(int i=1;i<=solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;

    // Set fully_implicit flag
    for(int i=1;i<=solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=fully_implicit;
    for(int i=1;i<=deformable_body_collection.deformables_forces.m;i++) deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=fully_implicit;
    for(int i=1;i<=solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=fully_implicit;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Write_Output_Files(const int frame) const
{
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    binding_list.Clamp_Particles_To_Embedded_Velocities();
    BASE::Write_Output_Files(frame);

    FILE_UTILITIES::Create_Directory(output_directory);
    std::string filename=STRING_UTILITIES::string_sprintf("%s/%s.%d",output_directory.c_str(),filename_base.c_str(),frame);

    data_exchange::deformable_body_simulation_output dbof;
    for(int i=1;i<=simulation_object_data.m;i++){
        if(simulation_object_data(i)->structure_index){
            data_exchange::deformable_body_output* dbo=new data_exchange::deformable_body_output;
            const ARRAY<int>& particle_map=simulation_object_data(i)->particle_map;
            From_Pb(dbo->position,particles.X.Subset(particle_map));
            From_Pb(dbo->velocity,particles.V.Subset(particle_map));
            dbof.simulation_data.push_back(dbo);}
        else dbof.simulation_data.push_back(0);}

    std::ofstream ofs(filename.c_str());
    if(write_text){
        boost::archive::text_oarchive oa(ofs);
        oa << dbof;}
    else{
        boost::archive::binary_oarchive oa(ofs);
        oa << dbof;}
}
//#####################################################################
// Function Add_Deformable_Body
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Add_Deformable_Body(const data_exchange::deformable_body& body,int body_index)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    BINDING_LIST<TV>& binding_list=deformable_body_collection.binding_list;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    ARRAY<int> parent_index;
    Triangulated_Surface_From_Data_Exchange(*surface,body.mesh,body.position,&parent_index);

    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/geometry.%d.tri",output_directory.c_str(),body_index),*surface);
    T density=body.mass?body.mass/surface->Volumetric_Volume():1000;

    TETRAHEDRALIZED_VOLUME<T>* tet_volume=0;
    TRIANGULATED_SURFACE<T>* tri_surface=0;
    ARRAY<int>& particle_map=simulation_object_data(body_index+1)->particle_map;
    tests.Create_Regular_Embedded_Surface(binding_list,soft_bindings,*surface,density,1000,(T)1e-5,particle_map,&tri_surface,&tet_volume,true);

    int structure=deformable_body_collection.deformable_geometry.structures.m;
    int enclosing_structure=structure-2;
    simulation_object_data(body_index+1)->structure_index=structure;
    simulation_object_data(body_index+1)->encosing_structure_index=enclosing_structure;

    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(*tri_surface,soft_bindings);
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Add_Rigid_Body(const data_exchange::ground_plane& body,int body_index)
{
    RIGID_BODY<TV>& ground=tests.Add_Ground();
    ground.X()=TV(To_Pb(body.position));
    ground.Rotation()=ROTATION<TV>::From_Rotated_Vector(TV(0,1,0),TV(To_Pb(body.normal)));
    simulation_object_data(body_index+1)->rigid_index=ground.particle_index;
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Add_Rigid_Body(const data_exchange::scripted_geometry& body,int body_index)
{
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create(particles);
    Triangulated_Surface_From_Data_Exchange(*surface,body.mesh,body.position,0);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/geometry.%d.tri",output_directory.c_str(),body_index+1),*surface);
    RIGID_BODY<TV>& rigid_body=*tests.Create_Rigid_Body_From_Triangulated_Surface(*surface,solid_body_collection.rigid_body_collection,1);
    LEVELSET_IMPLICIT_OBJECT<TV>* implicit_object=tests.Initialize_Implicit_Surface(*surface,20);
    rigid_body.Add_Structure(*implicit_object);
    rigid_body.is_static=true;
    solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    simulation_object_data(body_index+1)->rigid_index=rigid_body.particle_index;
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Add_Force(const data_exchange::gravity_force& force,int force_index)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    PHYSBAM_ASSERT(force.bodies_affected.size());
    ARRAY<int>* influenced_particles=new ARRAY<int>(particles.X.m);
    HASHTABLE<int> hash_particles;
    for(size_t i=0;i<force.bodies_affected.size();i++){
        int body=force.bodies_affected[i];
        int structure=simulation_object_data(body+1)->encosing_structure_index;
        PHYSBAM_ASSERT(structure);
        if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.deformable_geometry.structures(structure)))
            hash_particles.Set_All(volume->mesh.elements.Flattened());
        else PHYSBAM_FATAL_ERROR();}
    hash_particles.Get_Keys(*influenced_particles);
    DEFORMABLE_GRAVITY<TV>* gravity=new DEFORMABLE_GRAVITY<TV>(particles,influenced_particles,force.magnitude,TV(To_Pb(force.direction)));
    gravity->Own_Influenced_Particles();
    deformable_body_collection.Add_Force(gravity);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Add_Force(const data_exchange::volumetric_force& force,int force_index)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PHYSBAM_ASSERT(force.bodies_affected.size());
    for(size_t i=0;i<force.bodies_affected.size();i++){
        int body=force.bodies_affected[i];
        int structure=simulation_object_data(body+1)->encosing_structure_index;
        PHYSBAM_ASSERT(structure);
        if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.deformable_geometry.structures(structure))){
            solid_body_collection.Add_Force(Create_Finite_Volume(*volume,new COROTATED<T,3>(force.stiffness,force.poissons_ratio,force.damping),false));}
        else PHYSBAM_FATAL_ERROR();}
}
template class DEFORMABLE_EXAMPLE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLE_EXAMPLE<double>;
#endif
