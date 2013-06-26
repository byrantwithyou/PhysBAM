//#####################################################################
// Copyright 2006-2007, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_COLLISIONS_EXAMPLE
//#####################################################################
#ifndef __EMBEDDED_COLLISIONS_EXAMPLE__
#define __EMBEDDED_COLLISIONS_EXAMPLE__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
#include <Deformables/Bindings/PARTICLE_BINDING.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Fracture/EMBEDDING.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Dynamics/Meshing/RED_GREEN_TRIANGLES.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class EMBEDDED_COLLISIONS_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::solid_body_collection;using BASE::Set_External_Positions;using BASE::parse_args; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    RED_GREEN_TRIANGLES<TV>* redgreen;

    int maximum_number_of_boundary_refinements;
    T refinement_ratio;
    T sphere_scale;

    EMBEDDED_COLLISIONS_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,data_directory,solid_body_collection),redgreen(0),maximum_number_of_boundary_refinements(4),refinement_ratio(.5),
        sphere_scale(.5)
    {
    }

    virtual ~EMBEDDED_COLLISIONS_EXAMPLE()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
    parse_args->Add("-sphere_scale",&sphere_scale,"scale","sphere scale");
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
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    bool automatically_add_to_collision_structures=true;
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
        RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2,0))),true,true,1000);
    tetrahedralized_volume.Update_Number_Nodes();
    // initialize and refine boundary
    EMBEDDING<TV>& embedding=*new EMBEDDING<TV>(particles);
    tetrahedralized_volume.Initialize_Triangulated_Surface();
    embedding.material_surface_mesh.Initialize_Mesh(tetrahedralized_volume.triangulated_surface->mesh);
    delete tetrahedralized_volume.triangulated_surface;tetrahedralized_volume.triangulated_surface=0;
    TRIANGULATED_SURFACE<T>& triangulated_surface=embedding.material_surface;
    redgreen=new RED_GREEN_TRIANGLES<TV>(triangulated_surface);
    deformable_body_collection.Add_Structure(&embedding);

    tests.Add_Ground();
//    tests.Add_Rigid_Body("sphere",(T)sphere_scale,(T).5);
    tests.Add_Rigid_Body("sphere",(T)sphere_scale,(T).1);

    // add structures and rigid bodies to collisions
    if(automatically_add_to_collision_structures) deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    output_directory="Embedded_Collisions/output";
    last_frame=180;
    frame_rate=24;
    
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.implicit_solve_parameters.cg_iterations=500;
//        solids_parameters.cfl=.5;
    solids_parameters.write_static_variables_every_frame=true;

    Get_Initial_Data();

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,tetrahedralized_volume.mesh,0));
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).05,(T).25),
        true,(T).1));
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    COLLISION_BODY_COLLECTION<TV>& collision_body_list=deformable_body_collection.collisions.collision_body_list;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    EMBEDDING<TV>& embedding=deformable_body_collection.template Find_Structure<EMBEDDING<TV>&>();

    ARRAY<int>& segment_midpoints=redgreen->segment_midpoints;
    SEGMENT_MESH& segment_mesh=redgreen->segment_mesh;
    ARRAY<int> old_midpoints;
    for(int s=0;s<segment_midpoints.m;s++) if(segment_midpoints(s)){
        redgreen->Add_Free_Segment_Midpoint(segment_mesh.elements(s),segment_midpoints(s));
        old_midpoints.Append(segment_midpoints(s));}
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    tetrahedralized_volume.Initialize_Triangulated_Surface();
    redgreen->object.mesh.Initialize_Mesh(tetrahedralized_volume.triangulated_surface->mesh);
    redgreen->object.Update_Number_Nodes();
    redgreen->Initialize();
    delete tetrahedralized_volume.triangulated_surface;tetrahedralized_volume.triangulated_surface=0;
    TRIANGULATED_SURFACE<T>& triangulated_surface=redgreen->object;
    ARRAY<int> surface_particles;
    while(1){
        triangulated_surface.mesh.elements.Flattened().Get_Unique(surface_particles);
        ARRAY<T> particle_distances(particles.Size());
        particle_distances.Subset(surface_particles).Fill((T)FLT_MAX);
        for(int i=0;i<surface_particles.m;i++){int p=surface_particles(i);
            for(COLLISION_GEOMETRY_ID body(0);body<collision_body_list.bodies.m;body++)
                particle_distances(p)=PhysBAM::min(particle_distances(p),collision_body_list.bodies(body)->Implicit_Geometry_Extended_Value(particles.X(p)));}
        ARRAY<int> triangles_to_refine;
        for(int t=0;t<triangulated_surface.mesh.elements.m;t++){
            T triangle_distance=particle_distances.Subset(triangulated_surface.mesh.elements(t)).Min();
            T triangle_size=TRIANGLE_3D<T>(particles.X.Subset(triangulated_surface.mesh.elements(t))).Maximum_Edge_Length();
            int level=redgreen->leaf_levels_and_indices(t)(1);
            if(level<maximum_number_of_boundary_refinements && triangle_distance<triangle_size*refinement_ratio && redgreen) triangles_to_refine.Append(t);}
        if(!triangles_to_refine.m) break;
        redgreen->Refine_Simplex_List(triangles_to_refine);}
    LOG::cout<<"===== Boundary triangles : "<<triangulated_surface.mesh.elements.m<<" ====="<<std::endl;
    triangulated_surface.mesh.Initialize_Incident_Elements();

    // reinitialize bindings
    redgreen->Initialize_Segment_Index_From_Midpoint_Index();
    for(int i=0;i<old_midpoints.m;i++) if(!(*redgreen->segment_index_from_midpoint_index)(old_midpoints(i))) particles.Add_To_Deletion_List(old_midpoints(i));
    delete redgreen->free_segment_midpoints;redgreen->free_segment_midpoints=0;
    ARRAY<int> parents;ARRAY<T> weights;
    binding_list.Clean_Memory();
    triangulated_surface.mesh.elements.Flattened().Get_Unique(surface_particles);
    for(int i=0;i<surface_particles.m;i++){
        redgreen->Unrefined_Parents(surface_particles(i),parents,weights);
        switch(parents.m){
          case 1:
//            binding_list.Add_Binding(new PARTICLE_BINDING<TV>(particles,particles.Add_Element(),surface_particles(i)));
            break;
          case 2: 
            binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,surface_particles(i),VECTOR<int,2>(parents(1),parents(2)),VECTOR<T,2>(weights(1),weights(2))));
            break;
          case 3:
            binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(particles,surface_particles(i),VECTOR<int,3>(parents(1),parents(2),parents(3)),
                    VECTOR<T,3>(weights(1),weights(2),weights(3))));
            break;
          default: PHYSBAM_FATAL_ERROR();}} // should have strictly between 1-3 parents
    ARRAY<int> old_bound_particles;
    HASHTABLE<int,int> free_soft_bound_particles;
    for(int b=0;b<soft_bindings.bindings.m;b++){
        int bound_node,embedded_node;
        soft_bindings.bindings(b).Get(bound_node,embedded_node);
        free_soft_bound_particles.Insert(embedded_node,bound_node);
        old_bound_particles.Append(bound_node);}
    soft_bindings.Clean_Memory();
    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,soft_bindings,&free_soft_bound_particles);
    for(int i=0;i<old_bound_particles.m;i++) if(!soft_bindings.Particle_Is_Bound(old_bound_particles(i))) particles.Add_To_Deletion_List(old_bound_particles(i));

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
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(id==2) twist.linear=TV((T)-2,0,0);
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(id==2) frame.t=TV((T)5-2*time,(T).75,0);
}
//#####################################################################
};
}
#endif
