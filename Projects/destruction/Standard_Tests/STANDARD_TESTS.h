//#####################################################################
// Copyright 2008, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Prescored cube
//   2. Prescored pillar
//   3. Single cube
//   4. Cube stack
//   5. Single cylinder
//   6. Sphere pillar
//   7. Taller cube stack
//   8. Ball hitting wall
//   9. Fracture pattern viewer
//  10. Bunny
//  11. Gaussian points
//  12. Raining spheres
//  13. Friction test
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_PATTERN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_REGION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_CLUSTER_CONSTITUENT_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <climits>
#include <cstdio>

namespace PhysBAM{
template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;
    int parameter;
    bool prune_stacks_from_contact;
    bool prune_contact_using_velocity;
    bool use_nonanalytic_levelsets;

    TV sphere_initial_location;
    TV sphere2_initial_location;
    TV sphere_velocity;

    int fracture_pattern_index;
    T initial_speed_down_ramp;
    int kinematic_id;
    T ground_angle_rad;
    T mu;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    FRACTURE_PATTERN<T> fracture_pattern;
    ARRAY<PAIR<int,T> > rigid_body_clamp_time;
    T maximum_fall_speed;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::frame_rate;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),rigid_body_collection(solid_body_collection.rigid_body_collection)
    {
    }

    ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T curent_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");
    parse_args->Add_Integer_Argument("-fp",2,"specify fracture pattern");
    parse_args->Add_Option_Argument("-prunestacks","Do something quick for stacks during contact");
    parse_args->Add_Option_Argument("-velocityprune","Let collisions handle pairs with high velocity");
    parse_args->Add_Option_Argument("-noanalytic","disable analytic collisions");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Option_Argument("-noanalyticlevelsets","prevent usage of analytic levelsets");
    parse_args->Add_Option_Argument("-createpattern","create a fracture pattern");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);

    if(parse_args->Get_Option_Value("-createpattern")){
        if(test_number==1) Create_Box_Split_Pattern();
        else if(test_number==2) Create_Crossing_Planes_Pattern();
        else if(test_number==3) Create_Grain_Boundary_Surfaces();
        else if(test_number==4) Create_Pyramid_Pattern(RANGE<TV>(TV(-1,-1,-1),TV(1,1,1)),100,2);
        else if(test_number==5) Create_Wall_Pattern();
        else if(test_number==6) Create_Bunny_Pattern();
        else if(test_number==7) Create_Cylinder_Pattern();
        else if(test_number==8) Create_Raining_Spheres_Pattern();
        else Create_Pattern(test_number);
        exit(0);}

    frame_rate=30;

    parameter=parse_args->Get_Integer_Value("-parameter");
    if(parameter) output_directory+=STRING_UTILITIES::string_sprintf("_param%i",parameter);
    if(parse_args->Is_Value_Set("-fp")) fracture_pattern_index=parse_args->Get_Integer_Value("-fp");
    prune_stacks_from_contact=parse_args->Get_Option_Value("-prunestacks");
    prune_contact_using_velocity=parse_args->Get_Option_Value("-velocityprune");
    use_nonanalytic_levelsets=parse_args->Get_Option_Value("-noanalyticlevelsets");

    solids_parameters.rigid_body_collision_parameters.use_analytic_collisions=!parse_args->Get_Option_Value("-noanalytic");
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    solids_parameters.rigid_body_evolution_parameters.correct_evolution_energy=true;
    solids_parameters.rigid_body_collision_parameters.collision_iterations=10;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;
    solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_particle_partition=false;
    solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_triangle_hierarchy=false;
    solids_parameters.rigid_body_collision_parameters.use_fracture_particle_optimization=false;
    fracture_pattern.use_particle_partitions=solids_parameters.rigid_body_collision_parameters.use_fracture_particle_optimization;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    if(test_number==12) solids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
    else solids_parameters.rigid_body_collision_parameters.use_legacy_push_out=false;
    solids_parameters.use_rigid_deformable_contact=true;
    fluids_parameters.simulate=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    last_frame=240;
 
    switch(test_number){
        case 1: Prescore_Cube();break;
        case 2: Prescore_Pillar();break;
        case 3: Single_Cube();break;
        case 4: Cube_Stack();break;
        case 5: Single_Cylinder();break;
        case 6: Sphere_Pillar();break;
        case 7: Tall_Cube_Stack();break;
        case 8: Ball_Hitting_Wall();break;
        case 9: Display_Pattern();break;
        case 10: Poor_Bunny();break;
        case 11: Grain_Points();break;
        case 12: Raining_Spheres();break;
        case 13: Friction_Test();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct mass
    //binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    if(test_number!=9) tests.Add_Gravity();
}
//#####################################################################
// Function Display_Pattern
//#####################################################################
void Display_Pattern()
{
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;

    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fracture_pattern_index),fracture_pattern);
    for(int i=0;i<fracture_pattern.regions.m;i++){
        RIGID_BODY<TV>* new_body=new RIGID_BODY<TV>(rigid_body_collection,true);
        new_body->Add_Structure(*fracture_pattern.regions(i)->triangulated_surface);
        new_body->Add_Structure(*fracture_pattern.regions(i)->implicit_object);
        new_body->Frame().t-=TV(fracture_pattern.regions(i)->fracture_offset)*fracture_pattern.regions(i)->implicit_object->levelset.grid.DX();
        rigid_body_collection.Add_Rigid_Body_And_Geometry(new_body);}
}
//#####################################################################
// Function Prescore_Cube
//#####################################################################
void Prescore_Cube()
{
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;

    for(int i=0;i<89;i++){
        if(i==36 || i==71 || i==30) continue;
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body(STRING_UTILITIES::string_sprintf("Fractured_Cube/fragment.%02d",i),1,(T).5);(void)rigid_body;}


    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;
}
//#####################################################################
// Function Prescore_Pillar
//#####################################################################
void Prescore_Pillar()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;

    for(int i=0;i<94;i++){
        if(i==90) continue;
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body(STRING_UTILITIES::string_sprintf("Fractured_Pillar/fragment.%02d",i),1,(T).5);(void)rigid_body;
        if(rigid_body.Mass()<(T)1e-5) rigid_body_collection.rigid_body_particle.Remove_Body(rigid_body.particle_index);}

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;
}
//#####################################################################
// Function Single_Cube
//#####################################################################
void Single_Cube()
{
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;
    rigid_body_collection.rigid_geometry_collection.always_create_structure=true;
    solids_parameters.rigid_body_collision_parameters.use_fracture_pattern=true;

    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,(T).5,true,use_nonanalytic_levelsets);(void)rigid_body;
    rigid_body.Frame()=FRAME<TV>();
    rigid_body.fracture_threshold=(T)3.0;

    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fracture_pattern_index),fracture_pattern);

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;
}
//#####################################################################
// Function Cube_Stack
//#####################################################################
void Cube_Stack()
{
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;

    RIGID_BODY<TV>& rigid_body_1=tests.Add_Rigid_Body("box",1,(T).5,true,use_nonanalytic_levelsets);
    rigid_body_1.Frame().t=TV(0,(T)-.5,0);
    RIGID_BODY<TV>& rigid_body_2=tests.Add_Rigid_Body("box",(T).5,(T).5,true,use_nonanalytic_levelsets);
    rigid_body_2.Frame().t=TV(0,(T)1,0);
    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;
}
//#####################################################################
// Function Single_Cylinder
//#####################################################################
void Single_Cylinder()
{
    fracture_pattern_index=7;
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;
    solids_parameters.rigid_body_collision_parameters.use_fracture_pattern=true;

    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("cyllink",1,(T).5,true,use_nonanalytic_levelsets);(void)rigid_body;
    rigid_body.fracture_threshold=(T)1.0;

    RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere",(T).35,(T).5,true,true);
    sphere_initial_location=TV((T)15,(T)0,(T)0);
    sphere_velocity=TV((T)-10,(T)0,(T)0);
    sphere_body.Frame().t=sphere_initial_location;
    sphere_body.Is_Kinematic()=true;

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;

    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fracture_pattern_index),fracture_pattern);
}
//#####################################################################
// Function Sphere_Pillar
//#####################################################################
void Sphere_Pillar()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;

    for(int i=0;i<94;i++){
        if(i==90) continue;
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body(STRING_UTILITIES::string_sprintf("Fractured_Pillar/fragment.%02d",i),1,(T).5);(void)rigid_body;
        rigid_body.Update_Bounding_Box();T_ORIENTED_BOX oriented_box=rigid_body.Oriented_Bounding_Box();
        T min_side_length=FLT_MAX;TV saved_frame=rigid_body.Frame().t;
        for(int dim=0;dim<TV::m;dim++) min_side_length=min(min_side_length,oriented_box.edges.Column(dim).Magnitude());
        rigid_body_collection.rigid_body_particle.Remove_Body(rigid_body.particle_index);
        RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",min_side_length/(T)2.5,(T).5,true,use_nonanalytic_levelsets);sphere.Frame().t=saved_frame;
        if(sphere.Mass()<(T)1e-5) rigid_body_collection.rigid_body_particle.Remove_Body(sphere.particle_index);}
    rigid_body_collection.rigid_geometry_collection.Destroy_Unreferenced_Geometry();

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;
    RIGID_BODY<TV>& wall1=tests.Add_Ground((T).5,-2,1);
    wall1.Frame()=FRAME<TV>(TV((T)1.5,0,0),ROTATION<TV>((T)(pi/2),TV(0,0,1)));
    RIGID_BODY<TV>& wall2=tests.Add_Ground((T).5,-2,1);
    wall2.Frame()=FRAME<TV>(TV((T)-.5,0,0),ROTATION<TV>((T)(1.5*pi),TV(0,0,1)));
    RIGID_BODY<TV>& wall3=tests.Add_Ground((T).5,-2,1);
    wall3.Frame()=FRAME<TV>(TV(0,0,(T)-.5),ROTATION<TV>((T)(pi/2),TV(1,0,0)));
    RIGID_BODY<TV>& wall4=tests.Add_Ground((T).5,-2,1);
    wall4.Frame()=FRAME<TV>(TV(0,0,(T)1.5),ROTATION<TV>((T)(1.5*pi),TV(1,0,0)));
    (void)ground;
}
//#####################################################################
// Function Tall_Cube_Stack
//#####################################################################
void Tall_Cube_Stack()
{
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;

    for(int i=0;i<10;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",(T).5,(T).5,true,use_nonanalytic_levelsets);
        rigid_body.Frame().t=TV(0,(T)-.5+i,0);}
    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;
}
//#####################################################################
// Function Ball_Hitting_Wall
//#####################################################################
void Ball_Hitting_Wall()
{
    fracture_pattern_index=5;
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;
    solids_parameters.rigid_body_collision_parameters.use_fracture_pattern=true;

    TV edge_lengths((T)6,(T).4,(T)4);
    TV_INT dimensions(31,3,21);RANGE<TV> box((T)-.5*edge_lengths,(T).5*edge_lengths);TV dx=edge_lengths/TV(dimensions-1);
    TV_INT levelset_resolution(151,11,101);TV levelset_dx=edge_lengths/TV(levelset_resolution-1);int ghost_cells=2;
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);
    TRIANGULATED_SURFACE<T>& simplicial_object=*TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=simplicial_object.particles;particles.array_collection->Add_Elements(2*(dimensions.x*dimensions.z+dimensions.x*(dimensions.y-2)+(dimensions.z-2)*(dimensions.y-2)));
    TRIANGLE_MESH& triangle_mesh=simplicial_object.mesh;triangle_mesh.number_nodes=particles.array_collection->Size();
    int particle_index=1;
    for(int side=0;side<=1;side++){
        // Construct top and bottom
        for(int x=0;x<dimensions.x;x++) for(int z=0;z<dimensions.z;z++){
            particles.X(particle_index)=box.min_corner+dx*TV((T)x-1,side?(T)dimensions.y-1:0,(T)z-1);particle_index++;}
        for(int x=0;x<dimensions.x;x++) for(int y=2;y<=dimensions.y-1;y++){
            particles.X(particle_index)=box.min_corner+dx*TV((T)x-1,(T)y-1,side?(T)dimensions.z-1:0);particle_index++;}
        for(int z=2;z<=dimensions.z-1;z++) for(int y=2;y<=dimensions.y-1;y++){
            particles.X(particle_index)=box.min_corner+dx*TV(side?(T)dimensions.x-1:0,(T)y-1,(T)z-1);particle_index++;}}
    rigid_body->Add_Structure(simplicial_object);
    LEVELSET_IMPLICIT_OBJECT<TV>& implicit_object=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    implicit_object.levelset.grid.Initialize(levelset_resolution+2*ghost_cells,box.Thickened(levelset_dx.x*ghost_cells),false);
    implicit_object.levelset.phi.Resize(implicit_object.levelset.grid.Domain_Indices());implicit_object.levelset.phi.Fill(FLT_MAX);
    implicit_object.Update_Box();implicit_object.Update_Minimum_Cell_Size();
    for(NODE_ITERATOR iterator(implicit_object.levelset.grid);iterator.Valid();iterator.Next())
        implicit_object.levelset.phi(iterator.index)=box.Signed_Distance(iterator.Location());
    rigid_body->Add_Structure(implicit_object);
    rigid_body->Frame().r=ROTATION<TV>((T)-pi/2,TV(0,0,(T)1))*ROTATION<TV>((T)pi/2,TV(0,(T)1,0));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);

    RIGID_BODY<TV>* sphere_body;
    RIGID_BODY<TV>* sphere_body2;
    switch(parameter){
        case 1:
            rigid_body->fracture_threshold=(T)1.0;
            sphere_body=&tests.Add_Rigid_Body("sphere",(T).5,(T).5,true,true);
            sphere_initial_location=TV((T)20,(T).5,(T)1.5);
            break;
        case 2:
            rigid_body->fracture_threshold=(T)1.0;
            sphere_body=&tests.Add_Rigid_Body("sphere",(T).5,(T).5,true,true);
            sphere_initial_location=TV((T)20,(T)-.5,(T).25);
            break;
        case 3:
            rigid_body->fracture_threshold=(T)1.0;
            if(fracture_pattern_index==5) rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,0,(T)1))*ROTATION<TV>((T)pi/2,TV(0,(T)1,0));
            sphere_body=&tests.Add_Rigid_Body("sphere",(T).5,(T).5,true,true);
            sphere_initial_location=TV((T)15,(T)-1,(T)-2);
            break;
        case 4:
            rigid_body->fracture_threshold=FLT_MAX;
            sphere_body=&tests.Add_Rigid_Body("sphere_refined",(T).25,(T).5,true,true);
            sphere_body->Frame().r=ROTATION<TV>((T)-pi/2,TV(0,0,(T)1));
            sphere_initial_location=TV((T)20,(T)0,(T)0);
            sphere_body->fracture_threshold=(T).1;
            break;
        case 5:
            rigid_body->fracture_threshold=FLT_MAX;
            sphere_body=&tests.Add_Rigid_Body("sphere_refined",(T).25,(T).5,true,true);
            sphere_initial_location=TV((T)20,(T)0,(T)0);
            sphere_body->fracture_threshold=(T).1;
            break;
        case 6:
            rigid_body->fracture_threshold=(T)1.0;
            sphere_body=&tests.Add_Rigid_Body("sphere",(T).5,(T).5,true,true);
            sphere_initial_location=TV((T)10,(T)-1,(T)-2);
            sphere_body2=&tests.Add_Rigid_Body("sphere",(T).35,(T).5,true,true);
            sphere2_initial_location=TV((T)22,(T).5,(T)1.5);
            sphere_body2->Frame().t=sphere2_initial_location;
            sphere_body2->Is_Kinematic()=true;
            solids_parameters.rigid_body_collision_parameters.allow_refracturing=true;
            break;
        case 7:
            rigid_body->fracture_threshold=(T)3.0;
            sphere_body=&tests.Add_Rigid_Body("sphere",(T).5,(T).5,true,true);
            sphere_initial_location=TV((T)10,(T)0,(T)0);
            sphere_body2=&tests.Add_Rigid_Body("sphere",(T).5,(T).5,true,true);
            sphere2_initial_location=TV((T)-20,(T)0,(T)-1.5);
            sphere_body2->Frame().t=sphere2_initial_location;
            sphere_body2->Is_Kinematic()=true;
            solids_parameters.rigid_body_collision_parameters.allow_refracturing=true;
            break;
        default:
            PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized param number %d",parameter));}
    if(fracture_pattern_index==8) sphere_velocity=TV((T)-10,(T)0,(T)0);
    else sphere_velocity=TV((T)-15,(T)0,(T)0);
    sphere_body->Frame().t=sphere_initial_location;
    sphere_body->Is_Kinematic()=true;
    
    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;

    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fracture_pattern_index),fracture_pattern);
}
//#####################################################################
// Function Poor_Bunny
//#####################################################################
void Poor_Bunny()
{
    fracture_pattern_index=6;
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;

    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("bunny",3,(T).5,true,use_nonanalytic_levelsets);(void)rigid_body;
    rigid_body.Frame().t=TV(0,(T)-1,0);
    rigid_body.fracture_threshold=(T)10.0;

    RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere",(T).2,(T).5,true,true);
    sphere_initial_location=TV((T)0,(T)-1,(T)35);
    sphere_velocity=TV((T)0,(T)0,(T)-20);
    sphere_body.Frame().t=sphere_initial_location;
    sphere_body.Is_Kinematic()=true;

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-2,1);
    (void)ground;

    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fracture_pattern_index),fracture_pattern);
}
//#####################################################################
// Function Raining_Spheres
//#####################################################################
void Raining_Spheres()
{
    if(parameter==25) last_frame=240; // 25
    else last_frame=400*(parameter/100); // every 100 takes 400 seconds to fall
    fracture_pattern_index=8;
    RANDOM_NUMBERS<T> random;random.Set_Seed(100);
    solid_body_collection.print_diagnostics=false;
    solid_body_collection.print_residuals=false;
    solids_parameters.rigid_body_collision_parameters.use_fracture_pattern=true;

    TV edge_lengths((T)6,(T).4,(T)4);
    TV_INT dimensions(31,3,21);RANGE<TV> box((T)-.5*edge_lengths,(T).5*edge_lengths);TV dx=edge_lengths/TV(dimensions-1);
    TV_INT levelset_resolution(151,11,101);TV levelset_dx=edge_lengths/TV(levelset_resolution-1);int ghost_cells=2;
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);
    TRIANGULATED_SURFACE<T>& simplicial_object=*TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=simplicial_object.particles;particles.array_collection->Add_Elements(2*(dimensions.x*dimensions.z+dimensions.x*(dimensions.y-2)+(dimensions.z-2)*(dimensions.y-2)));
    TRIANGLE_MESH& triangle_mesh=simplicial_object.mesh;triangle_mesh.number_nodes=particles.array_collection->Size();
    int particle_index=1;
    for(int side=0;side<=1;side++){
        // Construct top and bottom
        for(int x=0;x<dimensions.x;x++) for(int z=0;z<dimensions.z;z++){
            particles.X(particle_index)=box.min_corner+dx*TV((T)x-1,side?(T)dimensions.y-1:0,(T)z-1);particle_index++;}
        for(int x=0;x<dimensions.x;x++) for(int y=2;y<=dimensions.y-1;y++){
            particles.X(particle_index)=box.min_corner+dx*TV((T)x-1,(T)y-1,side?(T)dimensions.z-1:0);particle_index++;}
        for(int z=2;z<=dimensions.z-1;z++) for(int y=2;y<=dimensions.y-1;y++){
            particles.X(particle_index)=box.min_corner+dx*TV(side?(T)dimensions.x-1:0,(T)y-1,(T)z-1);particle_index++;}}
    rigid_body->Add_Structure(simplicial_object);
    LEVELSET_IMPLICIT_OBJECT<TV>& implicit_object=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    implicit_object.levelset.grid.Initialize(levelset_resolution+2*ghost_cells,box.Thickened(levelset_dx.x*ghost_cells),false);
    implicit_object.levelset.phi.Resize(implicit_object.levelset.grid.Domain_Indices());implicit_object.levelset.phi.Fill(FLT_MAX);
    implicit_object.Update_Box();implicit_object.Update_Minimum_Cell_Size();
    for(NODE_ITERATOR iterator(implicit_object.levelset.grid);iterator.Valid();iterator.Next())
        implicit_object.levelset.phi(iterator.index)=box.Signed_Distance(iterator.Location());
    rigid_body->Add_Structure(implicit_object);
    rigid_body->Frame().t=TV(0,(T)-.2,0);
    rigid_body->fracture_threshold=(T)FLT_MAX;
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Mass(5);
    rigid_body->is_static=true;
    rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);

    RIGID_BODY<TV>& rigid_body2=tests.Add_Rigid_Body("cyllink",(T).4,(T).5,true,use_nonanalytic_levelsets);(void)rigid_body2;
    rigid_body2.Frame().t=TV((T)-2.6,(T)-1.2,(T)-1.5);
    rigid_body2.fracture_threshold=(T)FLT_MAX;
    rigid_body2.is_static=true;
    rigid_body2.Set_Mass(2.5);
    RIGID_BODY<TV>& rigid_body3=tests.Add_Rigid_Body("cyllink",(T).4,(T).5,true,use_nonanalytic_levelsets);(void)rigid_body3;
    rigid_body3.Frame().t=TV((T)-2.6,(T)-1.2,(T)1.5);
    rigid_body3.fracture_threshold=(T)FLT_MAX;
    rigid_body3.is_static=true;
    rigid_body3.Set_Mass(2.5);
    RIGID_BODY<TV>& rigid_body4=tests.Add_Rigid_Body("cyllink",(T).4,(T).5,true,use_nonanalytic_levelsets);(void)rigid_body4;
    rigid_body4.Frame().t=TV((T)2.6,(T)-1.2,(T)-1.5);
    rigid_body4.fracture_threshold=(T)FLT_MAX;
    rigid_body4.is_static=true;
    rigid_body4.Set_Mass(2.5);
    RIGID_BODY<TV>& rigid_body5=tests.Add_Rigid_Body("cyllink",(T).4,(T).5,true,use_nonanalytic_levelsets);(void)rigid_body5;
    rigid_body5.Frame().t=TV((T)2.6,(T)-1.2,(T)1.5);
    rigid_body5.fracture_threshold=(T)FLT_MAX;
    rigid_body5.is_static=true;
    rigid_body5.Set_Mass(2.5);

    RIGID_BODY<TV>& large_sphere_body=tests.Add_Rigid_Body("sphere",(T)1,(T).2,true,true);
    int world_height;
    if(parameter==25){
        large_sphere_body.Frame().t=TV((T)0,(T)60,(T)0);
        world_height=40;}
    else{
        world_height=100*(parameter/100); // 100 in height for each 100 spheres
        large_sphere_body.Frame().t=TV((T)0,(T)world_height+20,(T)0);}
    large_sphere_body.fracture_threshold=(T)FLT_MAX;
    large_sphere_body.Set_Coefficient_Of_Restitution(0);
    large_sphere_body.Set_Mass((T)10);
    ARRAY<ORIENTED_BOX<TV> > bounding_boxes;
    RANGE<TV> world=RANGE<TV>(TV(-(T)2.6,4,-(T)2.6),TV((T)2.6,(T)world_height,(T)2.6));
    T g=(T)9.8,h=world.min_corner.y,base_t=sqrt(2*h/g);
    int num_spheres=0;

    if(parameter==25){
        num_spheres=parameter;
        rigid_body_clamp_time.Append(PAIR<int,T>(large_sphere_body.particle_index,(T)4));}
    else{
        num_spheres=parameter;
        maximum_fall_speed=sqrt(3*g*h);
        rigid_body_clamp_time.Append(PAIR<int,T>(large_sphere_body.particle_index,(world_height-10)/maximum_fall_speed));}
    for(int i=0;i<num_spheres;i++){
        RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere_refined",(T).35,(T).2,true,true);
        sphere_body.Update_Bounding_Box();
        sphere_body.Frame()=Find_Placement(random,sphere_body.axis_aligned_bounding_box,bounding_boxes,world,false);
        sphere_body.Set_Mass((T)0.15);
        sphere_body.fracture_threshold=(T).25;
        if(maximum_fall_speed){
            T d=max((T)0,bounding_boxes.Last().Center().y-world.min_corner.y),clamp_time=d<h?0:base_t+(d-h)/maximum_fall_speed;
            if(clamp_time>0) rigid_body_clamp_time.Append(PAIR<int,T>(sphere_body.particle_index,clamp_time));}}

    tests.Add_Ground((T).5,-2,1);

    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fracture_pattern_index),fracture_pattern);
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
    if(solids_parameters.rigid_body_collision_parameters.use_fracture_pattern && collisions.fracture_pattern && !collisions.fracture_pattern->regions.m){
        for(int i=0;i<fracture_pattern.regions.m;i++) collisions.fracture_pattern->regions.Append(fracture_pattern.regions(i));}
    if(test_number==2){
        collisions.Set_Contact_Level_Iterations(2);
        collisions.Set_Contact_Pair_Iterations(20);
        solids_parameters.rigid_body_collision_parameters.contact_iterations=10;}
    collisions.prune_stacks_from_contact=prune_stacks_from_contact;
    collisions.prune_contact_using_velocity=prune_contact_using_velocity;
}
//#####################################################################
// Function Create_Pattern
//#####################################################################
void Create_Pattern(const int test_number)
{
    RANGE<TV> domain;int count=100;TV_INT pattern_center;
    rigid_body_collection.rigid_geometry_collection.always_create_structure=true;

    // Load in bodies and compute domain
    switch(test_number){
        case 2:
            {RANDOM_NUMBERS<T> rn;ARRAY<TV> seed_points;TV a(0,0,0);TV b(100,100,100);
            while(seed_points.m<10){
                TV rand_vec=rn.Get_Uniform_Vector(a,b);
                if(rand_vec.Normalize()<(T).01) continue;
                seed_points.Append(rand_vec);}
            GRID<TV> levelset_grid(TV_INT(100,100,100),RANGE<TV>(TV((T)-2.5,(T)-2.5,(T)-2.5),TV((T)2.5,(T)2.5,(T)2.5)),false);
            ARRAY<int,VECTOR<int,3> > regions(levelset_grid.Domain_Indices());regions.Fill(-1);
            ARRAY<TV_INT> neighbor_list;
            for(int s=0;s<seed_points.m;s++){
                TV_INT cell_index=levelset_grid.Cell(seed_points(s),0);
                regions(cell_index)=s;
                for(int axis=0;axis<3;axis++) for(int side=-1;side<=1;side+=2){
                    TV_INT neighbor_index=cell_index+side*TV_INT::Axis_Vector(axis);
                    if(levelset_grid.Domain_Indices().Lazy_Inside(neighbor_index) && regions(neighbor_index)==-1)
                        neighbor_list.Append(neighbor_index);}}
            while(neighbor_list.m>0){
                int rand_int=rn.Get_Uniform_Integer(1,neighbor_list.m);
                TV_INT cell_index=neighbor_list(rand_int);
                if(regions(cell_index)==-1){
                    ARRAY<int> filled_neighbors;
                    for(int axis=0;axis<3;axis++) for(int side=-1;side<=1;side+=2){
                        TV_INT neighbor_index=cell_index+side*TV_INT::Axis_Vector(axis);
                        if(levelset_grid.Domain_Indices().Lazy_Inside(neighbor_index)){
                            if(regions(neighbor_index)!=-1) filled_neighbors.Append(regions(neighbor_index));
                            else neighbor_list.Append(neighbor_index);}}
                    int region_index=rn.Get_Uniform_Integer(1,filled_neighbors.m);
                    regions(cell_index)=filled_neighbors(region_index);}
                neighbor_list.Remove_Index_Lazy(rand_int);}
            for(int r=0;r<seed_points.m;r++){
                LOG::cout << "constructing region " << r << std::endl;
                RANGE<TV_INT> local_counts=RANGE<TV_INT>::Zero_Box().Thickened(-INT_MAX);
                RANGE<TV> local_domain;
                for(NODE_ITERATOR iterator(levelset_grid);iterator.Valid();iterator.Next()){
                    if(regions(iterator.index)==r){
                        local_counts.Enlarge_To_Include_Point(iterator.index);
                        local_domain.Enlarge_To_Include_Point(iterator.Location());
                        for(int axis=0;axis<3;axis++) for(int side=-1;side<=1;side+=2){
                            TV_INT neighbor_index=iterator.index+side*TV_INT::Axis_Vector(axis);
                            local_counts.Enlarge_To_Include_Point(neighbor_index);
                            local_domain.Enlarge_To_Include_Point(levelset_grid.Node(neighbor_index));}}}
                GRID<TV> local_grid(local_counts.Edge_Lengths()+1,local_domain,false);
                ARRAY<T,VECTOR<int,3> > local_phi(local_grid.Domain_Indices());local_phi.Fill(FLT_MAX);
                LOG::cout << "local_grid " << local_grid << std::endl;
                ARRAY<TV_INT> local_initialized_cells;
                for(NODE_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next()){
                    TV_INT global_index=iterator.index+local_counts.min_corner-1;
                    if(levelset_grid.Domain_Indices().Lazy_Inside(global_index) && regions(global_index)==r){
                        local_phi(iterator.index)=-FLT_MAX;
                        for(int axis=0;axis<3;axis++) for(int side=-1;side<=1;side+=2){
                            TV_INT global_neighbor_index=global_index+side*TV_INT::Axis_Vector(axis);
                            TV_INT local_neighbor_index=iterator.index+side*TV_INT::Axis_Vector(axis);
                            if(!levelset_grid.Domain_Indices().Lazy_Inside(global_neighbor_index) || regions(global_neighbor_index)!=r){
                                if(local_phi(iterator.index)==FLT_MAX||local_phi(iterator.index)==-FLT_MAX){
                                    local_initialized_cells.Append(iterator.index);local_phi(iterator.index)=(T)-.5*levelset_grid.dX.x;}
                                if(local_phi(local_neighbor_index)==FLT_MAX||local_phi(local_neighbor_index)==-FLT_MAX){
                                    local_initialized_cells.Append(local_neighbor_index);local_phi(local_neighbor_index)=(T).5*levelset_grid.dX.x;}}}}}
                LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(local_grid,local_phi);
                lio->levelset.Fast_Marching_Method(0,0,&local_initialized_cells);
                TRIANGULATED_SURFACE<T>* trisurf=DUALCONTOUR_3D<T>::Create_Triangulated_Surface_From_Levelset(lio->levelset);
                int region_index=fracture_pattern.regions.Append(new FRACTURE_REGION<T>(trisurf,lio,false));
                fracture_pattern.regions(region_index)->fracture_offset=TV_INT(); // TODO need to find center and then compute offset correctly
                FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("levelset.%d.debug.phi",r),lio->levelset);
                FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("levelset.%d.debug.tri",r),*trisurf);}
            FILE_UTILITIES::Create_Directory(output_directory);
            FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/fracture_pattern.%d",output_directory.c_str(),test_number),fracture_pattern);
            return;}
        case 3:
            for(int i=0;i<89;i++){
                RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body(STRING_UTILITIES::string_sprintf("Fractured_Cube/fragment.%02d",i),1,(T).5);(void)rigid_body;
                rigid_body.Update_Bounding_Box();domain.Enlarge_To_Include_Box(rigid_body.axis_aligned_bounding_box);}
            pattern_center=TV_INT(25,25,25);
            break;
        default:
            PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
    T dx=domain.Edge_Lengths().Max()/count;
    TV expansion=(ceil(domain.Edge_Lengths()/dx)*dx-domain.Edge_Lengths())/2+2*dx;
    domain.Change_Size(expansion);
    TV_INT counts=TV_INT(rint(domain.Edge_Lengths()/dx));
    // TODO offset so that pattern_center is at a cell center (body 54, particle 19 location)

    GRID<TV> grid(counts,domain,false);

    LEVELSET_MAKER_UNIFORM<T> levelset_maker;
    levelset_maker.Verbose_Mode();
    for(int i=int(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(i);
        LEVELSET_IMPLICIT_OBJECT<TV>* fragment_implicit_object=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        TV_INT min_corner_index=grid.Clamped_Index(rigid_body.axis_aligned_bounding_box.min_corner)-TV_INT::All_Ones_Vector();grid.Clamp(min_corner_index);
        TV_INT max_corner_index=grid.Clamped_Index(rigid_body.axis_aligned_bounding_box.max_corner)+TV_INT::All_Ones_Vector();grid.Clamp(max_corner_index);
        RANGE<TV> clamped_domain=RANGE<TV>(grid.Node(min_corner_index),grid.Node(max_corner_index));
        TV_INT local_counts=max_corner_index-min_corner_index+TV_INT::All_Ones_Vector();
        for(int p=0;p<rigid_body.simplicial_object->particles.array_collection->Size();p++)
            rigid_body.simplicial_object->particles.X(p)=rigid_body.World_Space_Point(rigid_body.simplicial_object->particles.X(p))-clamped_domain.Center();
        rigid_body.simplicial_object->Initialize_Hierarchy();
        rigid_body.simplicial_object->Update_Bounding_Box();
        RANGE<TV> local_domain=RANGE<TV>(-TV(local_counts)*dx*(T).5,TV(local_counts)*dx*(T).5);

        fragment_implicit_object->levelset.grid.Initialize(local_counts,local_domain,false);
        fragment_implicit_object->levelset.phi.Resize(fragment_implicit_object->levelset.grid.Domain_Indices());
        fragment_implicit_object->Update_Box();fragment_implicit_object->Update_Minimum_Cell_Size();
        levelset_maker.Compute_Level_Set(*rigid_body.simplicial_object,fragment_implicit_object->levelset.grid,fragment_implicit_object->levelset.phi);
        int region_index=fracture_pattern.regions.Append(new FRACTURE_REGION<T>(rigid_body.simplicial_object,fragment_implicit_object,false));
        fracture_pattern.regions(region_index)->fracture_offset=(i==1)?TV_INT(25,0,0):TV_INT(-25,0,0);// pattern_center-min_corner_index+TV_INT::All_Ones_Vector();
    }

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/fracture_pattern.%d",output_directory.c_str(),test_number),fracture_pattern);
}
//#####################################################################
// Function Shrink_Levelset
//#####################################################################
void Shrink_Levelset(GRID<TV>& grid,ARRAY<T,VECTOR<int,3> >& phi,int boundary,TV_INT& center)
{
    TV_INT min_adjust;
    RANGE<TV_INT> inside=RANGE<TV_INT>::Zero_Box().Thickened(-INT_MAX),domain=grid.Domain_Indices();
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(phi(iterator.index)<0) inside.Enlarge_To_Include_Point(iterator.index);
    inside=RANGE<TV_INT>::Intersect(inside.Thickened(boundary),domain);
    min_adjust=inside.min_corner-domain.min_corner;
    center-=min_adjust;
    GRID<TV> grid_new(inside.Edge_Lengths()+1,RANGE<TV>(grid.Node(inside.min_corner),grid.Node(inside.max_corner)),false);
    ARRAY<T,VECTOR<int,3> > phi_new(grid_new.Domain_Indices());
    for(NODE_ITERATOR iterator(grid_new);iterator.Valid();iterator.Next()) phi_new(iterator.index)=phi(iterator.index+min_adjust);
    grid=grid_new;
    phi=phi_new;
}
//#####################################################################
// Function Create_Wall_Pattern
//#####################################################################
void Create_Wall_Pattern()
{
    FRACTURE_PATTERN<T> fp;
    TV_INT resolution(301,21,201);TV edge_lengths((T)6,(T).4,(T)4);TV dx=edge_lengths/TV(resolution-1);int ghost_cells=2;
    TV original_half_edge_length(edge_lengths/2);
    TV_INT actual_resolution=resolution+ghost_cells*2;TV actual_edge_lengths=TV(actual_resolution-1)*dx;
    TV half_edge_length(actual_edge_lengths/2);

    for(int i=2;i<=26;i++){
        GRID<TV>& local_grid=*new GRID<TV>(actual_resolution,RANGE<TV>(-half_edge_length,half_edge_length),false);
        ARRAY<T,VECTOR<int,3> >& local_phi=*new ARRAY<T,VECTOR<int,3> >(local_grid.Domain_Indices());local_phi.Fill(FLT_MAX);
        LEVELSET_IMPLICIT_OBJECT<TV>* refined_lio=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/wall/fragment.%d.phi",data_directory.c_str(),i),refined_lio->levelset);
        LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(local_grid,local_phi);
        for(NODE_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next())
            local_phi(iterator.index)=refined_lio->levelset.Extended_Phi(iterator.Location());
        TV_INT center_index=local_grid.Domain_Indices().Center();
        Shrink_Levelset(local_grid,local_phi,2,center_index);
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/wall/fragment.%d.tri",data_directory.c_str(),i),*surface);
        surface->Update_Number_Nodes();
        FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
        fr->fracture_offset=center_index;
        fp.regions.Append(fr);}
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),test_number),fp);
    return;
}
//#####################################################################
// Function Create_Bunny_Pattern
//#####################################################################
void Create_Bunny_Pattern()
{
    FRACTURE_PATTERN<T> fp;
    TV_INT resolution(71,71,71);TV edge_lengths((T)1.2,(T)1.2,(T).9);TV dx=edge_lengths/TV(resolution-1);int ghost_cells=2;
    TV original_half_edge_length(edge_lengths/2);
    TV_INT actual_resolution=resolution+ghost_cells*2;TV actual_edge_lengths=TV(actual_resolution-1)*dx;
    TV half_edge_length(actual_edge_lengths/2);

    for(int i=2;i<=49;i++){if(i==42 || i==47) continue;
        LOG::cout << "REGION " << i << std::endl;
        GRID<TV>& local_grid=*new GRID<TV>(actual_resolution,RANGE<TV>(-half_edge_length,half_edge_length),false);
        ARRAY<T,VECTOR<int,3> >& local_phi=*new ARRAY<T,VECTOR<int,3> >(local_grid.Domain_Indices());local_phi.Fill(FLT_MAX);
        LEVELSET_IMPLICIT_OBJECT<TV>* refined_lio=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/bunny/fragment.%d.phi",data_directory.c_str(),i),refined_lio->levelset);
        LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(local_grid,local_phi);
        for(NODE_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next())
            local_phi(iterator.index)=refined_lio->levelset.Extended_Phi(iterator.Location());
        TV_INT center_index=local_grid.Domain_Indices().Center();
        Shrink_Levelset(local_grid,local_phi,2,center_index);
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/bunny/fragment.%d.tri",data_directory.c_str(),i),*surface);
        surface->Update_Number_Nodes();
        FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
        fr->fracture_offset=center_index;
        fp.regions.Append(fr);}
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),test_number),fp);
    return;
}
//#####################################################################
// Function Create_Cylinder_Pattern
//#####################################################################
void Create_Cylinder_Pattern()
{
    FRACTURE_PATTERN<T> fp;
    TV_INT resolution(101,151,101);TV edge_lengths((T)2,(T)4,(T)2);TV dx=edge_lengths/TV(resolution-1);int ghost_cells=2;
    TV original_half_edge_length(edge_lengths/2);
    TV_INT actual_resolution=resolution+ghost_cells*2;TV actual_edge_lengths=TV(actual_resolution-1)*dx;
    TV half_edge_length(actual_edge_lengths/2);

    for(int i=2;i<=55;i++){if(i==33 || i==41 || i==46) continue;
        LOG::cout << "REGION " << i << std::endl;
        GRID<TV>& local_grid=*new GRID<TV>(actual_resolution,RANGE<TV>(-half_edge_length,half_edge_length),false);
        ARRAY<T,VECTOR<int,3> >& local_phi=*new ARRAY<T,VECTOR<int,3> >(local_grid.Domain_Indices());local_phi.Fill(FLT_MAX);
        LEVELSET_IMPLICIT_OBJECT<TV>* refined_lio=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/cylinder/fragment.%d.phi",data_directory.c_str(),i),refined_lio->levelset);
        LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(local_grid,local_phi);
        for(NODE_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next())
            local_phi(iterator.index)=refined_lio->levelset.Extended_Phi(iterator.Location());
        TV_INT center_index=local_grid.Domain_Indices().Center();
        Shrink_Levelset(local_grid,local_phi,2,center_index);
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/cylinder/fragment.%d.tri",data_directory.c_str(),i),*surface);
        surface->Update_Number_Nodes();
        FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
        fr->fracture_offset=center_index;
        fp.regions.Append(fr);}
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),test_number),fp);
    return;
}
//#####################################################################
// Function Create_Raining_Spheres_Pattern
//#####################################################################
void Create_Raining_Spheres_Pattern()
{
    FRACTURE_PATTERN<T> fp;
    TV_INT resolution(251,251,251);TV edge_lengths((T)2,(T)2,(T)2);TV dx=edge_lengths/TV(resolution-1);int ghost_cells=2;
    TV original_half_edge_length(edge_lengths/2);
    TV_INT actual_resolution=resolution+ghost_cells*2;TV actual_edge_lengths=TV(actual_resolution-1)*dx;
    TV half_edge_length(actual_edge_lengths/2);

    for(int i=2;i<=26;i++){if(i==25) continue;
        LOG::cout << "REGION " << i << std::endl;
        GRID<TV>& local_grid=*new GRID<TV>(actual_resolution,RANGE<TV>(-half_edge_length,half_edge_length),false);
        ARRAY<T,VECTOR<int,3> >& local_phi=*new ARRAY<T,VECTOR<int,3> >(local_grid.Domain_Indices());local_phi.Fill(FLT_MAX);
        LEVELSET_IMPLICIT_OBJECT<TV>* refined_lio=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/raining_spheres/fragment.%d.phi",data_directory.c_str(),i),refined_lio->levelset);
        LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(local_grid,local_phi);
        for(NODE_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next())
            local_phi(iterator.index)=refined_lio->levelset.Extended_Phi(iterator.Location());
        TV_INT center_index=local_grid.Domain_Indices().Center();
        Shrink_Levelset(local_grid,local_phi,2,center_index);
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/raining_spheres/fragment.%d.tri",data_directory.c_str(),i),*surface);
        surface->Update_Number_Nodes();
        FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
        fr->fracture_offset=center_index;
        fp.regions.Append(fr);}
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),test_number),fp);
    return;
}
//#####################################################################
// Function Create_Box_Split_Pattern
//#####################################################################
void Create_Box_Split_Pattern()
{
    FRACTURE_PATTERN<T> fp;
    for(int i=0;i<2;i++){
        int resolution=100;T edge_length=4;T dx=edge_length/resolution;int ghost_cells=2;
        TV half_edge_length;half_edge_length.Fill(edge_length/2);
        TV_INT counts;counts.Fill(resolution);half_edge_length.Fill(edge_length/2);
        GRID<TV> original_grid(counts,RANGE<TV>(-half_edge_length,half_edge_length),false);
        int actual_resolution=resolution+ghost_cells*2;T actual_edge_length=actual_resolution*dx;
        counts.Fill(actual_resolution);half_edge_length.Fill(actual_edge_length/2);
        GRID<TV> local_grid(counts,RANGE<TV>(-half_edge_length,half_edge_length),false);
        ARRAY<T,VECTOR<int,3> > local_phi(local_grid.Domain_Indices());local_phi.Fill(FLT_MAX);
        LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(local_grid,local_phi);
        for(NODE_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next())
            local_phi(iterator.index)=RANGE<TV>(original_grid.Domain()).Signed_Distance(iterator.Location());
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        for(NODE_ITERATOR iterator(original_grid,0,GRID<TV>::BOUNDARY_REGION);iterator.Valid();iterator.Next())
            surface->particles.X(surface->particles.array_collection->Add_Element())=iterator.Location();
        surface->Update_Number_Nodes();
        FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
        fr->fracture_offset=TV_INT(i==1?resolution-1+ghost_cells:ghost_cells,resolution/2-1+ghost_cells,resolution/2-1+ghost_cells);
        fp.regions.Append(fr);}
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),test_number),fp);
    return;
}
//#####################################################################
// Function Create_Pyramid_Pattern
//#####################################################################
void Create_Pyramid_Pattern(RANGE<TV> boundary,const int min_resolution,const int subdivision)
{
    T dx=(boundary.Edge_Lengths()/(T)min_resolution).Min();
    boundary.Scale_About_Center((rint(boundary.Edge_Lengths()/dx)*dx)/boundary.Edge_Lengths());
    ARRAY<TV> pts;
    for(int i=0;i<8;i++){
        TV pt=boundary.min_corner;
        for(int j=0;j<3;j++) if(i&(1<<j)) pt(j+1)=boundary.max_corner(j+1);
        pts.Append(pt);}
    TV center=boundary.Center();
    VECTOR<VECTOR<int,4>,6> corners(VECTOR<int,4>(0,1,3,2),VECTOR<int,4>(4,5,7,6),VECTOR<int,4>(0,1,5,4),VECTOR<int,4>(2,3,7,6),VECTOR<int,4>(0,2,6,4),VECTOR<int,4>(1,3,7,5));
    corners+=VECTOR<int,4>::All_Ones_Vector();
    TV jitter((T).001,(T).002,(T).003);
    FRACTURE_PATTERN<T> fp;
    for(int m=0;m<6;m++){
        GRID<TV>& grid=*new GRID<TV>(TV_INT(rint(boundary.Edge_Lengths()/dx))+1,boundary,false);
        ARRAY<T,VECTOR<int,3> >& phi=*new ARRAY<T,VECTOR<int,3> >(grid.Domain_Indices());
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        TV pt_inside=pts.Subset(corners(m)).Average();
        phi.Fill(-FLT_MAX);
        for(int i=0;i<4;i++){
            TRIANGLE_3D<T> tri(pts(corners(m)(i)),pts(corners(m)(i%4+1)),center);
            for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
                phi(iterator.index)=-min(-phi(iterator.index),tri.Distance_To_Triangle(iterator.Location()+jitter));}
        for(int i=0;i<4;i++){
            PLANE<T> plane(pts(corners(m)(i)),pts(corners(m)(i%4+1)),center);
            if(plane.Signed_Distance(pt_inside)>0) plane.normal=-plane.normal;
            for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
                if(plane.Signed_Distance(iterator.Location()+jitter)>0) phi(iterator.index)=abs(phi(iterator.index));}

        surface->particles.X(surface->particles.array_collection->Add_Element())=center;
        for(int i=0;i<4;i++){
            TV A=pts(corners(m)(i)),B=pts(corners(m)(i%4+1)),C=center;
            for(int j=0;j<subdivision;j++){
                for(int k=0;k<j;k++){
                    TV D=C+(B-C)*((T)k/j);
                    TV E=A+(D-A)*((T)j/subdivision);
                    surface->particles.X(surface->particles.array_collection->Add_Element())=E;}}}
        surface->Update_Number_Nodes();

        TV_INT center_index=grid.Domain_Indices().Center();
        Shrink_Levelset(grid,phi,2,center_index);
        LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(grid,phi);

        FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
        fr->fracture_offset=center_index;
        fp.regions.Append(fr);}
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),test_number),fp);
    for(int m=0;m<6;m++){
        delete &fp.regions(m)->implicit_object->levelset.grid;
        delete &fp.regions(m)->implicit_object->levelset.phi;}
}
//#####################################################################
// Function Create_Crossing_Planes_Pattern
//#####################################################################
void Create_Crossing_Planes_Pattern()
{
    FRACTURE_PATTERN<T> fp;
    int planes_per_side=5;
    TV_INT resolution(301,21,201);TV edge_lengths((T)6,(T).4,(T)4);TV dx=edge_lengths/TV(resolution-1);int ghost_cells=2;
    TV original_half_edge_length(edge_lengths/2);
    TV_INT actual_resolution=resolution+ghost_cells*2;TV actual_edge_lengths=TV(actual_resolution-1)*dx;
    TV half_edge_length(actual_edge_lengths/2);

    T x_diff=edge_lengths.x/(planes_per_side-1);
    int nodes_per_long_side=25;int nodes_per_short_side=(int)(nodes_per_long_side*(x_diff/edge_lengths.x));
    int nodes_per_static_side=5;
    TV jitter((T).001,(T).002,(T).003);
    for(int side=-1;side<=1;side+=2){
        for(int x_index=1;x_index<planes_per_side;x_index++){
            RANGE<TV> box=RANGE<TV>(-half_edge_length,half_edge_length);
            TV center=box.Center();
            GRID<TV>& grid=*new GRID<TV>(actual_resolution,box,false);
            ARRAY<T,VECTOR<int,3> >& phi=*new ARRAY<T,VECTOR<int,3> >(grid.Domain_Indices());phi.Fill(FLT_MAX);
            VECTOR<T,2> x_coords=VECTOR<T,2>((x_index-1)*x_diff-original_half_edge_length.x,x_index*x_diff-original_half_edge_length.x);
            PLANE<T> plane1(TV(x_coords(1),-original_half_edge_length.y,side*original_half_edge_length.z),TV(x_coords(1),original_half_edge_length.y,side*original_half_edge_length.z),center);
            PLANE<T> plane2(TV(x_coords(2),original_half_edge_length.y,side*original_half_edge_length.z),TV(x_coords(2),-original_half_edge_length.y,side*original_half_edge_length.z),center);
            TV pt_inside=TV(x_coords.Average(),0,side*original_half_edge_length.z);
            if(plane1.Signed_Distance(pt_inside)>0) plane1.normal=-plane1.normal;
            if(plane2.Signed_Distance(pt_inside)>0) plane2.normal=-plane2.normal;
            for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
                phi(iterator.index)=max(plane1.Signed_Distance(iterator.Location()+jitter),plane2.Signed_Distance(iterator.Location()+jitter));
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            for(int plane=0;plane<2;plane++){
                T min_x_coord=min(center.x,x_coords(plane));
                T max_x_coord=max(center.x,x_coords(plane));
                bool flipped=min_x_coord==center.x;
                T x_step=(max_x_coord-min_x_coord)/nodes_per_long_side;
                T z_step=(flipped?-side:side)*(original_half_edge_length.z-center.z)/nodes_per_long_side;
                T x=min_x_coord;
                T z=flipped?center.z:side*original_half_edge_length.z;
                for(int x_node=0;x_node<nodes_per_long_side;x_node++){
                    for(T y=-original_half_edge_length.y;y<=original_half_edge_length.y;y+=edge_lengths.y/nodes_per_static_side)
                        surface->particles.X(surface->particles.array_collection->Add_Element())=TV(x,y,z);
                    x+=x_step;z-=z_step;}}
            for(T x=x_coords(1);x<=x_coords(2);x+=(x_coords(2)-x_coords(1))/nodes_per_short_side)
                for(T y=-original_half_edge_length.y;y<=original_half_edge_length.y;y+=edge_lengths.y/nodes_per_static_side)
                    surface->particles.X(surface->particles.array_collection->Add_Element())=TV(x,y,side*original_half_edge_length.z);
            surface->Update_Number_Nodes();
            TV_INT center_index=grid.Domain_Indices().Center();
            Shrink_Levelset(grid,phi,5,center_index);
            LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(grid,phi);
            FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
            fr->fracture_offset=center_index;
            fp.regions.Append(fr);}}
    T z_x_ratio=((T)resolution.z)/resolution.x;
    planes_per_side=(int)(z_x_ratio*planes_per_side);
    T z_diff=edge_lengths.z/(planes_per_side-1);
    nodes_per_long_side=(int)(z_x_ratio*50);nodes_per_short_side=(int)(nodes_per_long_side*(z_diff/edge_lengths.z));
    for(int side=-1;side<=1;side+=2){
        for(int z_index=1;z_index<planes_per_side;z_index++){
            RANGE<TV> box=RANGE<TV>(-half_edge_length,half_edge_length);
            TV center=box.Center();
            GRID<TV>& grid=*new GRID<TV>(actual_resolution,box,false);
            ARRAY<T,VECTOR<int,3> >& phi=*new ARRAY<T,VECTOR<int,3> >(grid.Domain_Indices());phi.Fill(FLT_MAX);
            VECTOR<T,2> z_coords=VECTOR<T,2>((z_index-1)*z_diff-original_half_edge_length.z,z_index*z_diff-original_half_edge_length.z);
            PLANE<T> plane1(TV(side*original_half_edge_length.x,-original_half_edge_length.y,z_coords(1)),TV(side*original_half_edge_length.x,original_half_edge_length.y,z_coords(1)),center);
            PLANE<T> plane2(TV(side*original_half_edge_length.x,original_half_edge_length.y,z_coords(2)),TV(side*original_half_edge_length.x,-original_half_edge_length.y,z_coords(2)),center);
            TV pt_inside=TV(side*original_half_edge_length.x,0,z_coords.Average());
            if(plane1.Signed_Distance(pt_inside)>0) plane1.normal=-plane1.normal;
            if(plane2.Signed_Distance(pt_inside)>0) plane2.normal=-plane2.normal;
            for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
                phi(iterator.index)=max(plane1.Signed_Distance(iterator.Location()+jitter),plane2.Signed_Distance(iterator.Location()+jitter));
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            for(int plane=0;plane<2;plane++){
                T min_z_coord=min(center.z,z_coords(plane));
                T max_z_coord=max(center.z,z_coords(plane));
                bool flipped=min_z_coord==center.z;
                T z_step=(max_z_coord-min_z_coord)/nodes_per_long_side;
                T x_step=(flipped?-side:side)*(original_half_edge_length.x-center.x)/nodes_per_long_side;
                T z=min_z_coord;
                T x=flipped?center.x:side*original_half_edge_length.x;
                for(int z_node=0;z_node<nodes_per_long_side;z_node++){
                    for(T y=-original_half_edge_length.y;y<=original_half_edge_length.y;y+=edge_lengths.y/nodes_per_static_side)
                        surface->particles.X(surface->particles.array_collection->Add_Element())=TV(x,y,z);
                    z+=z_step;x-=x_step;}}
            for(T z=z_coords(1);z<=z_coords(2);z+=(z_coords(2)-z_coords(1))/nodes_per_short_side)
                for(T y=-original_half_edge_length.y;y<=original_half_edge_length.y;y+=edge_lengths.y/nodes_per_static_side)
                    surface->particles.X(surface->particles.array_collection->Add_Element())=TV(side*original_half_edge_length.x,y,z);
            surface->Update_Number_Nodes();
            TV_INT center_index=grid.Domain_Indices().Center();
            Shrink_Levelset(grid,phi,5,center_index);
            LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(grid,phi);
            FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
            fr->fracture_offset=center_index;
            fp.regions.Append(fr);}}
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),test_number),fp);
}
//#####################################################################
// Function Create_Grain_Boundary_Surfaces
//#####################################################################
void Create_Grain_Boundary_Surfaces()
{
    RANGE<TV> wall_box(TV((T)-3,(T)-.2,(T)-2),TV((T)3,(T).2,(T)2));
    TRIANGULATED_SURFACE<T>* wall_surface=TESSELLATION::Generate_Triangles(wall_box);
    for(int i=0;i<3;i++)
        wall_surface->Linearly_Subdivide();
    FILE_UTILITIES::Write_To_File(stream_type,"wall.tri",*wall_surface);

    RANGE<TV> cylinder_box(TV((T)-1,(T)-2,(T)-1),TV((T)1,(T)2,(T)1));
    TRIANGULATED_SURFACE<T>* cylinder_surface=TESSELLATION::Generate_Triangles(cylinder_box);
    for(int i=0;i<3;i++)
        cylinder_surface->Linearly_Subdivide();
    FILE_UTILITIES::Write_To_File(stream_type,"cylinder.tri",*cylinder_surface);

    RANGE<TV> bunny_box(TV((T)-.6,(T)-.6,(T)-.45),TV((T).6,(T).6,(T).45));
    TRIANGULATED_SURFACE<T>* bunny_surface=TESSELLATION::Generate_Triangles(bunny_box);
    for(int i=0;i<3;i++)
        bunny_surface->Linearly_Subdivide();
    FILE_UTILITIES::Write_To_File(stream_type,"bunny.tri",*bunny_surface);

    RANGE<TV> raining_box(TV((T)-2,(T)-2,(T)-2),TV((T)2,(T)2,(T)2));
    TRIANGULATED_SURFACE<T>* raining_surface=TESSELLATION::Generate_Triangles(raining_box);
    for(int i=0;i<3;i++)
        raining_surface->Linearly_Subdivide();
    FILE_UTILITIES::Write_To_File(stream_type,"raining.tri",*raining_surface);
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(test_number==8 && (parameter==4 || parameter==5) && frame==59){
        rigid_body_collection.Rigid_Body(2).Is_Kinematic()=false;
        solid_body_collection.Update_Simulated_Particles();}
    else if(test_number==10 && frame==10)
        solids_parameters.rigid_body_collision_parameters.use_fracture_pattern=true;
    else if(test_number==12){
        if(rigid_body_collection.rigid_body_particle.frame(6).t.y<=5 && rigid_body_collection.Rigid_Body(5).is_static){
            for(int id=int(1);id<=int(5);id++)
              rigid_body_collection.Rigid_Body(id).is_static=false;
            rigid_body_collection.Rigid_Body(int(1)).fracture_threshold=(T)3.0;
            solid_body_collection.Update_Simulated_Particles();
        }}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==5 && id==int(2))
        twist.linear=sphere_velocity;
    else if(test_number==8){
        if(id==int(2))
            twist.linear=sphere_velocity;
        else if(parameter==6 && id==int(3))
            twist.linear=sphere_velocity;
        else if(parameter==7 && id==int(3))
            twist.linear=-sphere_velocity;}
    else if(test_number==10 && id==int(2))
        twist.linear=sphere_velocity;
    else if(test_number==13 && id==kinematic_id){
        T magnitude=max((T)0,initial_speed_down_ramp+time*(T)9.8*(sin(ground_angle_rad)-mu*cos(ground_angle_rad)));
        TV direction(-cos(ground_angle_rad),-sin(ground_angle_rad),0);
        twist.linear=magnitude*direction;}
    else return false;
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==5 && id==int(2))
        frame.t=sphere_initial_location+time*sphere_velocity;
    else if(test_number==8){
        if(id==int(2))
            frame.t=sphere_initial_location+time*sphere_velocity;
        else if(parameter==6 && id==int(3))
            frame.t=sphere2_initial_location+time*sphere_velocity;
        else if(parameter==7 && id==int(3))
            frame.t=sphere2_initial_location-time*sphere_velocity;}
    else if(test_number==10 && id==int(2))
        frame.t=sphere_initial_location+time*sphere_velocity;
    else if(test_number==13 && id==kinematic_id){
        T magnitude;
        T ug=(T)9.8*(sin(ground_angle_rad)-mu*cos(ground_angle_rad));
        if(initial_speed_down_ramp+time*ug<0) magnitude=-(T).5*sqr(initial_speed_down_ramp)/ug;
        else magnitude=initial_speed_down_ramp*time+(T).5*sqr(time)*ug;
        TV direction(-cos(ground_angle_rad),-sin(ground_angle_rad),0);
        frame.t=magnitude*direction-TV(0,0,(T).25);}
}
//#####################################################################
// Function Grain_Points
//#####################################################################
void Grain_Points()
{
    last_frame=1;
    RANDOM_NUMBERS<T> random;
    RANGE<TV> box;
    if(parameter==1) box=RANGE<TV>(TV(-3,-(T).2,-2),TV(3,(T).2,2));
    else if(parameter==2) box=RANGE<TV>(TV(-.5,-.5,-.5),TV(.5,(T).5,(T).5));
    else if(parameter==3) box=RANGE<TV>(TV(-(T).2,-(T).2,-(T).15),TV((T).2,(T).2,(T).15));
    else if(parameter==4) box=RANGE<TV>(TV(-(T)2,-(T)2,-(T)2),TV((T)2,(T)2,(T)2));
    else PHYSBAM_FATAL_ERROR();

    // Parameter abuse.
    int num_pts=fracture_pattern_index;
    T sigma=0; // TODO: Use different parameter.

    FILE* F=fopen(STRING_UTILITIES::string_sprintf("%s/points.txt",output_directory.c_str()).c_str(),"w");
    PHYSBAM_ASSERT(F);

    fprintf(F,"%d\n",num_pts);
    for(int i=0;i<num_pts;i++){
        TV x;
        x.x=random.Get_Gaussian();
        x.y=random.Get_Gaussian();
        x.z=random.Get_Gaussian();
        x=x*sigma+box.Center();
        if(box.Lazy_Inside(x)) fprintf(F,"%22.16g %22.16g %22.16g\n",x.x,x.y,x.z);
        else i--;}
    fclose(F);

    exit(0);
}
//#####################################################################
// Function Friction_Test
//#####################################################################
void Friction_Test()
{
    RIGID_BODY<TV>* rigid_body=0;

    T ground_angle_deg=10;
    ground_angle_rad=ground_angle_deg*((T)pi/180);
    mu=(T).2;
    T box_size=(T)0.05;
    initial_speed_down_ramp=(T)1;

    LOG::cout<<"Using angle "<<ground_angle_deg<<" deg ("<<ground_angle_rad<<") (critical coeff="<<tan(ground_angle_rad)<<std::endl;

    FRAME<TV> frame(TV(-box_size*sin(ground_angle_rad),box_size*cos(ground_angle_rad),0),ROTATION<TV>(ground_angle_rad,TV(0,0,1)));
    TWIST<TV> twist(frame.r.Rotate(TV(-initial_speed_down_ramp,0,0)),TV());

    rigid_body=&tests.Add_Rigid_Body("box",box_size,mu);
    rigid_body->Frame()=frame;
    rigid_body->Twist()=twist;
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Name("box");

    frame.t.z+=(T).2;
    rigid_body=&tests.Add_Rigid_Body("sphere",(T).05,0);
    rigid_body->Frame().t=frame.t;
    rigid_body->Set_Name("analytic solution");
    rigid_body->Is_Kinematic()=true;
    kinematic_id=rigid_body->particle_index;

    RIGID_BODY<TV>& ground=tests.Add_Ground(mu,0,0);
    ground.Set_Coefficient_Of_Rolling_Friction(1);
    ground.Frame().r=ROTATION<TV>(ground_angle_rad,TV(0,0,1));
}
//#####################################################################
// Function Find_Placement
//#####################################################################
FRAME<TV> Find_Placement(RANDOM_NUMBERS<T>& random,const RANGE<TV>& bounding_box,ARRAY<ORIENTED_BOX<TV> >& bounding_boxes,const RANGE<TV>& world,bool want_rotate)
{
    for(int i=0;i<10000;i++){
        FRAME<TV> frame;
        if(want_rotate) frame.r=random.template Get_Rotation<TV>();
        ORIENTED_BOX<TV> oriented_box(bounding_box,frame.r);
        RANGE<TV> new_box(oriented_box.Bounding_Box());
        frame.t=random.Get_Uniform_Vector(world.min_corner-new_box.min_corner,world.max_corner-new_box.max_corner);
        oriented_box.corner+=frame.t;
        bool okay=true;
        for(int j=0;j<bounding_boxes.m;j++) if(oriented_box.Intersection(bounding_boxes(j))){okay=false;break;}
        if(okay){
            bounding_boxes.Append(oriented_box);
            return frame;}}
    PHYSBAM_FATAL_ERROR("Could not find suitable placement");
}
//#####################################################################
// Function Advance_One_Time_Step_End_Callback
//#####################################################################
void Advance_One_Time_Step_End_Callback(const T dt,const T time)
{
    if(test_number==12 && maximum_fall_speed){
        for(int i=0;i<rigid_body_clamp_time.m;i++) if(rigid_body_clamp_time(i).y>=time){int id=rigid_body_clamp_time(i).x;
            if(rigid_body_collection.Is_Active(id)){
                RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);
                T speed=rigid_body.Twist().linear.Magnitude();
                if(speed>maximum_fall_speed) rigid_body.Twist().linear*=maximum_fall_speed/speed;}}}
}
//#####################################################################
};
}
#endif
