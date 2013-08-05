//#####################################################################
// Copyright 2007, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_HARD_BINDING_EXAMPLE
//#####################################################################
#ifndef __SIMPLE_HARD_BINDING_EXAMPLE__
#define __SIMPLE_HARD_BINDING_EXAMPLE__

#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
#include <Deformables/Bindings/PARTICLE_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/IMPLICIT_ZERO_LENGTH_SPRINGS.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Dynamics/Meshing/RED_GREEN_TRIANGLES.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SIMPLE_HARD_BINDING_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    TRIANGULATED_SURFACE<T>* surface;
    RED_GREEN_TRIANGLES<TV>* redgreen;
    int subsamples;
    T sphere_scale;
    RANDOM_NUMBERS<T> random_numbers;
    ARRAY<ARRAY<int> > triangle_free_particles;
    bool dynamic_subsampling;
    int refinement_level;
    ARRAY<BINDING<TV>*> redgreen_bindings;
    ARRAY<T> refinement_distance;

    SIMPLE_HARD_BINDING_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,data_directory,solid_body_collection),surface(0),redgreen(0),subsamples(6),sphere_scale((T).05),dynamic_subsampling(false),
        refinement_level(3),refinement_distance(2)
    {
    }

    virtual ~SIMPLE_HARD_BINDING_EXAMPLE()
    {delete surface;delete redgreen;}

    // Unused callbacks
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}    
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
    parse_args->Add("-subsamples",&subsamples,"value","subsamples");
    parse_args->Add("-sphere_scale",&sphere_scale,"value","sphere scale");
    parse_args->Add("-radius",&refinement_distance(1),"value","radius");
    parse_args->Add("-dynamic",&dynamic_subsampling,"dynamic");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
    tests.data_directory=data_directory;
}    
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Redgreen
//#####################################################################
void Initialize_Redgreen()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    output_directory="Simple_Hard_Binding/output";
    last_frame=120;
    frame_rate=24;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-5;
    solids_parameters.implicit_solve_parameters.cg_iterations=1000;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.write_static_variables_every_frame=true;
    refinement_distance(1)=refinement_distance(2)=(T).04;
    random_numbers.Set_Seed(1234);
    fluids_parameters.simulate=false;

    redgreen=new RED_GREEN_TRIANGLES<TV>(*surface);
    redgreen->object.Update_Number_Nodes();
    redgreen->Initialize();
    for(int i=0;i<refinement_level;i++) redgreen->Refine_Simplex_List(ARRAY<int>(IDENTITY_ARRAY<>(surface->mesh.elements.m)));
    redgreen->Initialize_Segment_Index_From_Midpoint_Index();
    triangle_free_particles.Resize(surface->mesh.elements.m,true);

    ARRAY<int> surface_particles;surface->mesh.elements.Flattened().Get_Unique(surface_particles);
    for(int i=0;i<surface_particles.m;i++){int p=surface_particles(i);
        ARRAY<int> parents_list;ARRAY<T> weights_list;
        redgreen->Unrefined_Parents(p,parents_list,weights_list);
        if(parents_list.m<1||parents_list.m>3) PHYSBAM_FATAL_ERROR();
        if(parents_list.m==1){
            int parent=parents_list(1);
            redgreen_bindings.Append(new PARTICLE_BINDING<TV>(particles,p,parent));}
        else if(parents_list.m==2){
            VECTOR<int,2> parents(parents_list(1),parents_list(2));VECTOR<T,2> weights(weights_list(1),weights_list(2));
            redgreen_bindings.Append(new LINEAR_BINDING<TV,2>(particles,p,parents,weights));}
        else{
            VECTOR<int,3> parents(parents_list(1),parents_list(2),parents_list(3));VECTOR<T,3> weights(weights_list(1),weights_list(2),weights_list(3));
            redgreen_bindings.Append(new LINEAR_BINDING<TV,3>(particles,p,parents,weights));}}
}
//#####################################################################
// Function Subsample_Surface
//#####################################################################
void Subsample_Surface()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    FREE_PARTICLES<TV>& free_particles=deformable_body_collection.template Find_Structure<FREE_PARTICLES<TV>&>();

    // Sprinkle more points on the boundary
    for(int i=0;i<subsamples*surface->mesh.elements.m;i++){
        int triangle=random_numbers.Get_Uniform_Integer(1,surface->mesh.elements.m);
        VECTOR<T,3> weights;do{weights=random_numbers.Get_Uniform_Vector(TV(),TV(1,1,1));}while(weights.x+weights.y>=1);weights.z=(T)1-weights.x-weights.y;
        int hard_bound_particle=particles.Add_Element();
        binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(particles,hard_bound_particle,surface->mesh.elements(triangle),weights));
        int soft_bound_particle=particles.Add_Element();
        soft_bindings.Add_Binding(VECTOR<int,2>(soft_bound_particle,hard_bound_particle),true);
        free_particles.nodes.Append(soft_bound_particle);}
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    // Initialize simulation geometry
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
        RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1.2,0))),true,true,1000);
    tetrahedralized_volume.Update_Number_Nodes();tetrahedralized_volume.Initialize_Triangulated_Surface();
    tetrahedralized_volume.triangulated_surface->mesh.Initialize_Incident_Elements();
    surface=tetrahedralized_volume.triangulated_surface;
    tetrahedralized_volume.triangulated_surface=0;tetrahedralized_volume.mesh.boundary_mesh=0;
//    deformable_body_collection.Add_Structure(&tetrahedralized_volume);

    FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create(particles);
    deformable_body_collection.Add_Structure(&free_particles);
    
    if(!dynamic_subsampling) Subsample_Surface();
    else Initialize_Redgreen();

    tests.Add_Ground();
    tests.Add_Rigid_Body("sphere",sphere_scale,(T).1);

    // add structures and rigid bodies to collisions
    if(subsamples==0) deformable_body_collection.collisions.collision_structures.Append(&tetrahedralized_volume);
    else deformable_body_collection.collisions.collision_structures.Append(&free_particles);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    Get_Initial_Data();

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,tetrahedralized_volume.mesh,0));
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)1e6,(T).475,(T).01,(T).25),
        true,(T).1));
    solid_body_collection.Update_Simulated_Particles();
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(false);
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Velocities(false);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(id==int(2)) twist.linear=TV((T)-2,0,0);
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(id==int(2)) frame.t=TV((T)5-2*time,(T).75,0);
}
//#####################################################################
// Function Parent_Nodes
//#####################################################################
VECTOR<int,3> Parent_Nodes(const int leaf_triangle) const
{
    int level,triangle;redgreen->leaf_levels_and_indices(leaf_triangle).Get(level,triangle);
    for(;level>1;level--) triangle=(*redgreen->parent(level))(triangle);
    return redgreen->meshes(1)->elements(triangle);
}
//#####################################################################
// Function Add_Subsamples
//#####################################################################
void Add_Subsamples(const int triangle,ARRAY<BINDING<TV>*>& new_binding_list,ARRAY<VECTOR<int,2> >& new_soft_bindings)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    for(int i=0;i<subsamples;i++){
        VECTOR<T,3> weights;do{weights=random_numbers.Get_Uniform_Vector(TV(),TV(1,1,1));}while(weights.x+weights.y>=1);weights.z=(T)1-weights.x-weights.y;
        // add hard binding
        int hard_bound_particle=particles.Add_Element_From_Deletion_List();
        // get the parents in the coarse mesh
        VECTOR<int,3> parents=Parent_Nodes(triangle);
        TV location=TRIANGLE_3D<T>(particles.X.Subset(surface->mesh.elements(triangle))).Point_From_Barycentric_Coordinates(weights);
        // get the weights relative to the parents
        weights=TRIANGLE_3D<T>(particles.X.Subset(parents)).Barycentric_Coordinates(location);
        new_binding_list.Append(new LINEAR_BINDING<TV,3>(particles,hard_bound_particle,parents,weights));
        // add soft binding
        int soft_bound_particle=particles.Add_Element_From_Deletion_List();
        new_soft_bindings.Append(VECTOR<int,2>(soft_bound_particle,hard_bound_particle));
        // add to triangle_free_particles
        triangle_free_particles(triangle).Append(soft_bound_particle);}
}
//#####################################################################
// Function Delete_Subsamples
//#####################################################################
void Delete_Subsamples(const int triangle)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    
    // only need to delete the particles
    for(int i=0;i<triangle_free_particles(triangle).m;i++){
        int soft_bound_particle=triangle_free_particles(triangle)(i);
        int hard_bound_particle=soft_bindings.Parent(soft_bound_particle);
        particles.Add_To_Deletion_List(soft_bound_particle);
        particles.Add_To_Deletion_List(hard_bound_particle);}

    triangle_free_particles(triangle).Clean_Memory();
}
//#####################################################################
// Function Persist_Subsamples
//#####################################################################
void Persist_Subsamples(const int triangle,ARRAY<BINDING<TV>*>& new_binding_list,ARRAY<VECTOR<int,2> >& new_soft_bindings)
{
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    
    for(int i=0;i<triangle_free_particles(triangle).m;i++){
        int soft_bound_particle=triangle_free_particles(triangle)(i);
        // save soft binding
        const VECTOR<int,2>& soft_binding=soft_bindings.bindings(soft_bindings.Soft_Binding(soft_bound_particle));
        new_soft_bindings.Append(soft_binding);
        // save hard binding
        int hard_binding_index=binding_list.binding_index_from_particle_index(soft_binding.y);
        new_binding_list.Append(binding_list.bindings(hard_binding_index));
        binding_list.bindings(hard_binding_index)=0;} // so that binding won't be deleted when binding_list is rebuilt
}
//#####################################################################
// Function Preprocess_Solids_Substep
//#####################################################################
void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    if(!dynamic_subsampling) return;

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    COLLISION_BODY_COLLECTION<TV>& collision_body_list=deformable_body_collection.collisions.collision_body_list;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    FREE_PARTICLES<TV>& free_particles=deformable_body_collection.template Find_Structure<FREE_PARTICLES<TV>&>();
    ARRAY<VECTOR<int,3> >& surface_elements=surface->mesh.elements;

    ARRAY<BINDING<TV>*> new_binding_list;
    ARRAY<VECTOR<int,2> > new_soft_bindings;

    // surface particles
    ARRAY<int> surface_particles;surface_elements.Flattened().Get_Unique(surface_particles);

    // clamp binding list and redgreen bindings
    binding_list.Clamp_Particles_To_Embedded_Positions();
    for(int b=0;b<redgreen_bindings.m;b++) redgreen_bindings(b)->Clamp_To_Embedded_Position();

    // the minimum distance of each particle to a collision object
    ARRAY<ARRAY<T>,COLLISION_GEOMETRY_ID> particle_distances(collision_body_list.Size());
    for(COLLISION_GEOMETRY_ID i(0);i<particle_distances.Size();i++){
        particle_distances(i).Resize(particles.Size());
        INDIRECT_ARRAY<ARRAY<T>,ARRAY<int>&> subset=particle_distances(i).Subset(surface_particles);subset.Fill((T)FLT_MAX);}
    for(int i=0;i<surface_particles.m;i++){int p=surface_particles(i);
        for(COLLISION_GEOMETRY_ID body(0);body<collision_body_list.Size();body++)
            particle_distances(body)(p)=PhysBAM::min(particle_distances(body)(p),collision_body_list(body).Implicit_Geometry_Extended_Value(particles.X(p)));}

    // iterate over surface elements
    for(int t=0;t<surface_elements.m;t++){
        T triangle_distance1=particle_distances(COLLISION_GEOMETRY_ID(1)).Subset(surface_elements(t)).Min();
        T triangle_distance2=particle_distances(COLLISION_GEOMETRY_ID(2)).Subset(surface_elements(t)).Min();
        if(triangle_distance1<refinement_distance(1)||triangle_distance2<refinement_distance(2)) // triangle should be subsampled
            if(triangle_free_particles(t).m) Persist_Subsamples(t,new_binding_list,new_soft_bindings);
            else Add_Subsamples(t,new_binding_list,new_soft_bindings);
        else Delete_Subsamples(t);}

    // rebuild free particles
    free_particles.nodes.Clean_Memory();
    for(int t=0;t<triangle_free_particles.m;t++) free_particles.nodes.Append_Elements(triangle_free_particles(t));

    // rebuild hard bindings
    binding_list.Clean_Memory();
    for(int b=0;b<new_binding_list.m;b++) binding_list.Add_Binding(new_binding_list(b));

    // rebuild soft bindings
    soft_bindings.Clean_Memory();
    for(int b=0;b<new_soft_bindings.m;b++) soft_bindings.Add_Binding(new_soft_bindings(b),true);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();
    
    // correct mass
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // reinitialize collisions
    deformable_body_collection.collisions.Initialize_Object_Collisions(solids_parameters.deformable_object_collision_parameters.collide_with_interior,
        solids_parameters.deformable_object_collision_parameters.collision_tolerance,
        solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects,
        solids_parameters.deformable_object_collision_parameters.disable_multiple_levelset_collisions,
        solids_parameters.deformable_object_collision_parameters.maximum_levelset_collision_projection_velocity);

    // reinitialize fragments
    solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
};
}
#endif
