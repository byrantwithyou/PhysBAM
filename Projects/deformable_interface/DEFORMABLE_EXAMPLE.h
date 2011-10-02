//#####################################################################
// Copyright 2011, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_EXAMPLE
//#####################################################################
#ifndef __DEFORMABLE_EXAMPLE__
#define __DEFORMABLE_EXAMPLE__

#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
#include "../data_exchange/deformable_body_simulation.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
namespace PhysBAM{

template<class T,int d>
VECTOR<T,d> To_Pb(const data_exchange::fixed_vector<T,d>& v)
{
    VECTOR<T,d> x;
    for(size_t i=0; i<d; i++) x(i+1)=v.data[i];
    return x;
}

template<class TV>
void To_Pb(ARRAY_VIEW<TV> view,const std::vector<data_exchange::vf3>& v)
{
    for(size_t i=0; i<v.size(); i++) view(i+1)=To_Pb(v[i]);
}

template<class TV>
void To_Pb(ARRAY<TV>& array,const std::vector<data_exchange::vf3>& v)
{
    array.Resize(v.size());
    To_Pb(array, v);
}

void Triangle_Mesh_From_Data_Exchange(TRIANGLE_MESH& mesh, const data_exchange::polygon_mesh& poly, ARRAY<int>* parent_index)
{
    for(size_t i=0, k=0; i<poly.polygon_counts.size(); i++){
        int n=poly.polygon_counts[i];
        for(int j=2; j<n; j++){
            mesh.elements.Append(VECTOR<int,3>(poly.polygons[k], poly.polygons[k+j-1], poly.polygons[k+j]));
            if(parent_index) parent_index->Append(i);}
        k+=n;}
}

template<class TV>
void Triangulated_Surface_From_Data_Exchange(TRIANGULATED_SURFACE<TV>& surface, const data_exchange::polygon_mesh& poly, const std::vector<data_exchange::vf3>& pos, ARRAY<int>* parent_index)
{
    Triangle_Mesh_From_Data_Exchange(surface.mesh, poly, parent_index);
    surface.particles.array_collection->Add_Elements(pos.size());
    To_Pb(surface.particles.X, pos);
}

template<class T_input>
class DEFORMABLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    T stiffness_multiplier;
    T damping_multiplier;
    T bending_stiffness_multiplier;
    T bending_damping_multiplier;
    RANDOM_NUMBERS<T> random_numbers;
    bool fully_implicit;
    std::string scene_file;
    std::string filename_base;
    bool read_text;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;using BASE::parse_args;

    DEFORMABLE_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),fully_implicit(false),scene_file("scene"),read_text(false)
    {
    }

    ~DEFORMABLE_EXAMPLE()
    {}

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_End_Callback(const T dt,const T time) {}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return true;}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Option_Argument("-fully_implicit","use fully implicit forces");
    parse_args->Add_String_Argument("-scene","file specifying simulation scene");
    parse_args->Add_String_Argument("-read_text","scene written in text mode");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    fully_implicit=parse_args->Is_Value_Set("-fully_implicit");
    scene_file=parse_args->Get_String_Value("-scene");
    read_text=parse_args->Is_Value_Set("-read_text");
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    random_numbers.Set_Seed(1234);
    fluids_parameters.simulate=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness=(T)1e-3;
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    solids_parameters.verbose_dt=true;
    solids_parameters.triangle_collision_parameters.total_collision_loops=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=20;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
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

    ARRAY<int> body_index_to_structure_index;
    ARRAY<int> body_index_to_encosing_structure_index;
    ARRAY<int> body_index_to_rigid_index;

    for(size_t i=0; i<dbs.simulation_objects.size(); i++){
        if(data_exchange::deformable_body* body=dynamic_cast<data_exchange::deformable_body*>(dbs.simulation_objects[i])){
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            ARRAY<int> parent_index;
            Triangulated_Surface_From_Data_Exchange(*surface,body->mesh,body->position,&parent_index);
            T density=body->mass?body->mass/surface->Volumetric_Volume():1000;

            ARRAY<int> particle_map;
            TETRAHEDRALIZED_VOLUME<T>* tet_volume=0;
            TRIANGULATED_SURFACE<T>* tri_surface=0;
            tests.Create_Regular_Embedded_Surface(binding_list,*surface,density,512,(T)1e-5,&particle_map,&tri_surface,&tet_volume);

            int enclosing_structure=deformable_body_collection.deformable_geometry.Add_Structure(tet_volume);
            int structure=deformable_body_collection.deformable_geometry.Add_Structure(tri_surface);
            body_index_to_rigid_index.Append(0);
            body_index_to_structure_index.Append(structure);
            body_index_to_encosing_structure_index.Append(enclosing_structure);}
        else if(data_exchange::ground_plane* body=dynamic_cast<data_exchange::ground_plane*>(dbs.simulation_objects[i])){
            RIGID_BODY<TV>& ground=tests.Add_Ground();
            ground.X()=To_Pb(body->position);
            ground.Rotation()=ROTATION<TV>::From_Rotated_Vector(TV(0,1,0),To_Pb(body->normal));
            body_index_to_rigid_index.Append(ground.particle_index);
            body_index_to_structure_index.Append(0);
            body_index_to_encosing_structure_index.Append(0);}
        else if(data_exchange::scripted_geometry* body=dynamic_cast<data_exchange::scripted_geometry*>(dbs.simulation_objects[i])){
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            Triangulated_Surface_From_Data_Exchange(*surface,body->mesh,body->position,0);
            RIGID_BODY<TV>& rigid_body=*tests.Create_Rigid_Body_From_Triangulated_Surface(*surface,solid_body_collection.rigid_body_collection,1);
            body_index_to_rigid_index.Append(rigid_body.particle_index);
            body_index_to_structure_index.Append(0);
            body_index_to_encosing_structure_index.Append(0);}}

    PHYSBAM_ASSERT(body_index_to_rigid_index.m==(int)dbs.simulation_objects.size());
    PHYSBAM_ASSERT(body_index_to_structure_index.m==(int)dbs.simulation_objects.size());
    PHYSBAM_ASSERT(body_index_to_encosing_structure_index.m==(int)dbs.simulation_objects.size());

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    for(size_t i=0;i<dbs.simulation_forces.size();i++){
        if(data_exchange::gravity_force* force=dynamic_cast<data_exchange::gravity_force*>(dbs.simulation_forces[i])){
            PHYSBAM_ASSERT(force->bodies_affected.size());
            ARRAY<int>* influenced_particles=new ARRAY<int>(particles.X.m);
            HASHTABLE<int> hash_particles;
            for(size_t i=0;i<force->bodies_affected.size();i++){
                int body=force->bodies_affected[i];
                int structure=body_index_to_encosing_structure_index(body);
                PHYSBAM_ASSERT(structure);
                if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.deformable_geometry.structures(structure)))
                    hash_particles.Set_All(volume->mesh.elements.Flattened());
                else PHYSBAM_FATAL_ERROR();}
            hash_particles.Get_Keys(*influenced_particles);
            DEFORMABLE_GRAVITY<TV>* gravity=new DEFORMABLE_GRAVITY<TV>(particles,influenced_particles,force->magnitude,To_Pb(force->direction));
            gravity->Own_Influenced_Particles();
            deformable_body_collection.Add_Force(gravity);}
        else if(data_exchange::volumetric_force* force=dynamic_cast<data_exchange::volumetric_force*>(dbs.simulation_forces[i])){
            PHYSBAM_ASSERT(force->bodies_affected.size());
            for(size_t i=0;i<force->bodies_affected.size();i++){
                int body=force->bodies_affected[i];
                int structure=body_index_to_encosing_structure_index(body);
                PHYSBAM_ASSERT(structure);
                if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.deformable_geometry.structures(structure))){
                    solid_body_collection.Add_Force(Create_Finite_Volume(*volume,new NEO_HOOKEAN<T,3>(force->stiffness,force->poissons_ratio,force->damping),false));}
                else PHYSBAM_FATAL_ERROR();}}}

    last_frame=dbs.number_frames;
    this->frame_rate=1/dbs.dt;
    output_directory=dbs.output_directory;
    filename_base=dbs.filename_base;

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

////////////////

    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++){
        deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures(i));
        if(solids_parameters.triangle_collision_parameters.perform_self_collision && !dynamic_cast<FREE_PARTICLES<TV>*>(deformable_body_collection.deformable_geometry.structures(i)))
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
};
}
#endif

