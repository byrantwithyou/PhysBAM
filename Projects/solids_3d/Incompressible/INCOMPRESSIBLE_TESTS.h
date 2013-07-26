//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_TESTS
//#####################################################################
//   1. Ping pong
//   2. Squish
//   3. Hammer
//   4. Splat!
//   5. Sphere drop
//   6. Flatten
//   7. Buddha gears
//   8. Buddha fills container
//   9. Sphere fills container
//  10. Armadillo puddle
//  11. Expanding mattress
//  12. Sphere and stairs
//  13. Two Spheres
//  14. Armadillo gears
//  15. Sphere gears
//  16. Armadillo and stairs
//  17. Pile of tori
//  18. Pile of tori in a cage
//  19. Torus falling in awesome collisions
//  20. Hanging cloth
//#####################################################################
#ifndef __INCOMPRESSIBLE_TESTS__
#define __INCOMPRESSIBLE_TESTS__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Math_Tools/integer_log.h>
#include <Tools/Matrices/QUATERNION.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/RIGID_BODY_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Constitutive_Models/SPLINE_MODEL.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/BINDING_SPRINGS.h>
#include <Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <climits>
namespace PhysBAM{

template<class T_input>
class INCOMPRESSIBLE_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::Time_At_Frame;using BASE::stream_type;using BASE::solid_body_collection;using BASE::test_number;using BASE::parse_args;using BASE::mpi_world;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    T test_poissons_ratio;
    T hittime;
    T recovery_time;
    T minimum_volume_recovery_time_scale;
    bool high_resolution;
    T stiffen;
    int max_cg_iterations;
    T solids_cg_tolerance;
    ROTATION<TV> initial_orientation;
    TV initial_velocity,initial_angular_velocity;
    ROTATION<TV> roller_orientation;
    T roller_speed,roller_friction;
    T cylinder_time;
    TV cylinder_start,cylinder_velocity;
    T ground_friction;
    bool merge_at_boundary;
    bool use_soft_bindings;
    bool use_tet_collisions;
    bool use_volumetric_self_collisions;
    bool use_neohookean;
    bool use_neumann;
    bool self_collide_with_interior_nodes;
    ARRAY<int> fixed_particles;
    ARRAY<TV> fixed_particle_positions;
    INTERPOLATION_CURVE<T,TV> curve,curve2;

    T youngs_modulus;
    T hardening_deformation;
    T hardening_strength;
    T internal_poissons_ratio;

    int tori_m,tori_n,tori_mn;
    T tori_base_height;
    T tori_max_angular_velocity;
    ARRAY<RIGID_BODY_STATE<TV> > tori_initial_states;

    INCOMPRESSIBLE_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,data_directory,solid_body_collection),test_poissons_ratio((T).5),hittime(1),
        minimum_volume_recovery_time_scale(0),high_resolution(false),stiffen(1),max_cg_iterations(20),solids_cg_tolerance((T)1e-3),
        ground_friction(0),merge_at_boundary(false),tori_m(2),tori_n(2),tori_mn(2),tori_base_height((T)1.5),tori_max_angular_velocity(2)
    {
    }

    virtual ~INCOMPRESSIBLE_TESTS()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
    parse_args->Add("-poisson",&test_poissons_ratio,"","poisson's ratio for test 24");
    parse_args->Add("-stiffen",&stiffen,"","stiffness multiplier for various tests");
    parse_args->Add_Not("-noself",&solids_parameters.triangle_collision_parameters.perform_self_collision,"disable self-collisions");
    parse_args->Add("-hittime",&hittime,"time","time required for splat");
    parse_args->Add("-hires",&high_resolution,"use high resolution objects");
    parse_args->Add("-cgincomp",&max_cg_iterations,"iterations","maximum number of CG iterations for incompressible pressure solves");
    parse_args->Add("-project",&merge_at_boundary,"Combine boundary one-rings with neighbors");
    parse_args->Add("-cgsolids",&solids_cg_tolerance,"tol","CG tolerance for backward Euler");
    parse_args->Add("-min_volume_recovery_time_scale",&minimum_volume_recovery_time_scale,"value","Minimum dt over which entire volume may be recovered");
    parse_args->Add_Not("-noneumann",&use_neumann,"Do not apply Neumann boundary conditions");
    parse_args->Add("-ground_friction",&ground_friction,"value","Ground friction");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
    tests.data_directory=data_directory;
    LOG::cout<<"Running Incompressible Test Number "<<test_number<<std::endl;
    output_directory=STRING_UTILITIES::string_sprintf("Incompressible/Test_%d",test_number);
    if(frame_rate!=24) output_directory+=STRING_UTILITIES::string_sprintf("_fr%g",frame_rate);
    if(minimum_volume_recovery_time_scale) output_directory+=STRING_UTILITIES::string_sprintf("_mvrts%g",minimum_volume_recovery_time_scale);
    
    if(test_poissons_ratio!=(T).5) output_directory+=STRING_UTILITIES::string_sprintf("_p%g",test_poissons_ratio);
    if(hittime!=1) output_directory+=STRING_UTILITIES::string_sprintf("_ht%g",hittime);
    if(high_resolution) output_directory+="_hires";
    LOG::Stat("stiffen",stiffen);
    if(stiffen!=1) output_directory+=STRING_UTILITIES::string_sprintf("_stiff%g",stiffen);
    if(max_cg_iterations!=20) output_directory+=STRING_UTILITIES::string_sprintf("_cgi%d",max_cg_iterations);
    if(abs(solids_cg_tolerance-(T)1e-3)>(T)1e-7) output_directory+=STRING_UTILITIES::string_sprintf("_cgs%g",solids_cg_tolerance);
    if(abs(solids_cg_tolerance-(T)1e-3)>(T)1e-7 && (test_number==7 || test_number==8 || test_number==10)) solids_cg_tolerance=(T)1e-2;
    if(merge_at_boundary) output_directory+="_bound";
    if(!use_neumann) output_directory+="_noneumann";
    if(ground_friction) output_directory+=STRING_UTILITIES::string_sprintf("_fric%g",ground_friction);
    
    if(hittime!=1) tori_n=(int)hittime;

    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    LOG::cout<<"perform self collisions: "<<solids_parameters.triangle_collision_parameters.perform_self_collision<<std::endl;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Collide_With_Soft_Bound_Copy
//#####################################################################
void Collide_With_Soft_Bound_Surface_Copy(TETRAHEDRALIZED_VOLUME<T>& volume,bool rigid_body_soft)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    volume.Update_Number_Nodes();
    volume.Initialize_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create(deformable_body_collection.particles);
    surface->mesh.Initialize_Mesh(volume.triangulated_surface->mesh);
    tests.Substitute_Soft_Bindings_For_Nodes(*surface,solid_body_collection.deformable_body_collection.soft_bindings);
    deformable_body_collection.Add_Structure(surface);
    if(rigid_body_soft) deformable_body_collection.collisions.collision_structures.Append(surface);
    else deformable_body_collection.collisions.collision_structures.Append(&volume);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(surface);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    max_cg_iterations=500;
    solid_body_collection.deformable_body_collection.triangle_repulsions.hierarchy_repulsion_thickness_multiplier=2;
    
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)0;
    solids_parameters.triangle_collision_parameters.repulsions_youngs_modulus=solids_parameters.triangle_collision_parameters.collisions_final_repulsion_youngs_modulus=(T)3e15;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).03;
    solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
    solids_parameters.triangle_collision_parameters.allow_intersections=false;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
    self_collide_with_interior_nodes=false;

    solids_parameters.cfl=(T)2.0;
    solids_parameters.implicit_solve_parameters.cg_tolerance=solids_cg_tolerance;
    merge_at_boundary=false;
    use_soft_bindings=false;
    use_tet_collisions=false;
    use_volumetric_self_collisions=false;
    use_neohookean=false;
    
    youngs_modulus=(T)2e4*stiffen;
    hardening_deformation=(T).5;
    hardening_strength=4;
    internal_poissons_ratio=test_poissons_ratio<(T).5?test_poissons_ratio:0;
    
    // From BUDDHA_EXAMPLE.h
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1;//3e-3
    solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).25;
    solids_parameters.triangle_collision_parameters.collisions_collision_thickness=(T)1e-5;
    solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=8;

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    std::string sphere_filename=data_directory+"/Tetrahedralized_Volumes/";
    sphere_filename+=high_resolution?"sphere_103k.tet":"sphere.tet";
    std::string torus_filename=data_directory+"/Tetrahedralized_Volumes/";
    torus_filename+=high_resolution?"torus_115K.tet":"torus_12K.tet";
    solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
    bool rigid_body_soft=false;
    std::string buddha_filename=data_directory+"/Tetrahedralized_Volumes/";
    T buddha_resize;
    if(high_resolution){
        buddha_filename+="buddha_fine.tet";
        buddha_resize=1;}
    else{
        buddha_filename+="buddha_241k.tet";
        buddha_resize=2;}
    std::string armadillo_filename=data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet";
    if(high_resolution && test_number==14){LOG::cout<<"No high resolution for armadillo simulations."<<std::endl;PHYSBAM_FATAL_ERROR();}
    T armadillo_scale=(T).0085;
    // default material parameters
    T binding_stiffness=(T)1e5;
    bool use_gravity=false;
    RIGID_BODY<TV>* body=0;

    switch(test_number){
        case 1:{
            last_frame=(int)(8*frame_rate);
            body=&tests.Add_Rigid_Body("box",(T)2,(T)0);body->Frame().t.z=4;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            body=&tests.Add_Rigid_Body("box",(T)2,(T)0);body->Frame().t.z=-4;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,false,1000);
            curve.Add_Control_Point((T)0,TV(0,0,(T)4));
            curve.Add_Control_Point((T).2,TV(0,0,(T)0));
            curve.Add_Control_Point((T)1,TV(0,0,(T)2));
            curve.Add_Control_Point((T)1.8,TV(0,0,(T)2));
            curve.Add_Control_Point((T)2,TV(0,0,(T)-4));
            curve.Add_Control_Point((T)3,TV(0,0,(T)2));
            curve.Add_Control_Point((T)3.9,TV(0,0,(T)2));
            curve.Add_Control_Point((T)4,TV(0,0,(T)-7.8));
            curve.Add_Control_Point((T)5,TV(0,0,(T)-2));
            curve.Add_Control_Point((T)6,TV(0,0,(T)16));
            curve.Add_Control_Point((T)7,TV(0,0,(T)-3));
            curve.Add_Control_Point((T)8,TV(0,0,(T)0));
            curve2.Add_Control_Point((T)0,TV(0,0,(T)-12));
            curve2.Add_Control_Point((T)1,TV(0,0,(T)-12));
            curve2.Add_Control_Point((T)2,TV(0,0,(T)-12));
            curve2.Add_Control_Point((T)3,TV(0,0,(T)-12));
            curve2.Add_Control_Point((T)4,TV(0,(T)0,(T)-35));
            curve2.Add_Control_Point((T)4.20,TV(0,(T)0,(T)-35));
            curve2.Add_Control_Point((T)4.6,TV(0,(T)0,(T)-35));
            curve2.Add_Control_Point((T)7,TV(0,0,(T)-12));
            curve2.Add_Control_Point((T)8,TV(0,0,(T)-12));
            break;}
        case 2:{
            last_frame=(int)(8*frame_rate);
            body=&tests.Add_Rigid_Body("box",(T)10,(T)0);body->Frame().t.z=11;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            body=&tests.Add_Rigid_Body("box",(T)10,(T)0);body->Frame().t.z=-11;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,false,1000);
            curve.Add_Control_Point((T)0,TV(0,0,(T)11));
            curve.Add_Control_Point((T)4,TV(0,0,(T)10));
            curve.Add_Control_Point((T)6,TV(0,0,(T)10));
            curve.Add_Control_Point((T)10,TV(0,0,(T)11));
            curve2.Add_Control_Point((T)0,TV(0,0,(T)-11));
            curve2.Add_Control_Point((T)4,TV(0,0,(T)-10));
            curve2.Add_Control_Point((T)6,TV(0,0,(T)-10));
            curve2.Add_Control_Point((T)10,TV(0,0,(T)-11));
            break;}
        case 3:{
            last_frame=(int)(8*frame_rate);
            body=&tests.Add_Rigid_Body("box",(T)10,(T)0);body->Frame().t.z=11;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            body=&tests.Add_Rigid_Body("box",(T)10,(T)0);body->Frame().t.z=-11;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,false,1000);
            int N=32;
            for(int i=0;i<N;i++){
                curve.Add_Control_Point((T)(2*i)/(N/(T)4),TV(0,0,(T)3));
                curve.Add_Control_Point((T)(2*i+1)/(N/(T)4),TV(0,0,(T)1));
                curve2.Add_Control_Point((T)(2*i)/(N/(T)4),TV(0,0,(T)-3));
                curve2.Add_Control_Point((T)(2*i+1)/(N/(T)4),TV(0,0,(T)-1));}
            curve.Add_Control_Point((T)8,TV(0,0,(T)3));
            curve2.Add_Control_Point((T)8,TV(0,0,(T)-3));
            break;}
        case 4:{
            last_frame=(int)(2*frame_rate);
            body=&tests.Add_Rigid_Body("box",(T)1,(T)0);body->Frame().t.z=11;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            body=&tests.Add_Rigid_Body("box",(T)1,(T)0);body->Frame().t.z=-11;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,false,1000);
            curve.Add_Control_Point((T)0,TV(0,0,(T)3));
            curve.Add_Control_Point((T)hittime,TV(0,0,(T)1));
            curve.Add_Control_Point((T)2,TV(0,0,(T)1));
            curve2.Add_Control_Point((T)0,TV(0,0,(T)-3));
            curve2.Add_Control_Point((T)hittime,TV(0,0,(T)-1));
            curve2.Add_Control_Point((T)2,TV(0,0,(T)-1));
            break;}
        case 5:{
            last_frame=(int)(5*frame_rate);
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,false,1000);
            tests.Add_Ground(0);

            use_gravity=true;
            use_neohookean=true;
            hardening_deformation=(T).01;
            hardening_strength=(T).25;
            break;}
        case 6:{
            last_frame=(int)(8*frame_rate);
            body=&tests.Add_Rigid_Body("box",(T)3,(T)0);body->Frame().t.z=5;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            body=&tests.Add_Rigid_Body("box",(T)3,(T)0);body->Frame().t.z=-5;
            rigid_body_collection.rigid_body_particles.kinematic(body->particle_index)=true;
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,false,1000);
            curve.Add_Control_Point((T)0,TV(0,0,(T)5));
            curve.Add_Control_Point((T)3,TV(0,0,(T)3));
            curve.Add_Control_Point((T)5,TV(0,0,(T)3));
            curve.Add_Control_Point((T)8,TV(0,0,(T)5));
            curve2.Add_Control_Point((T)0,TV(0,0,(T)-5));
            curve2.Add_Control_Point((T)3,TV(0,0,(T)-3));
            curve2.Add_Control_Point((T)5,TV(0,0,(T)-3));
            curve2.Add_Control_Point((T)8,TV(0,0,(T)-5));
            break;}
        case 7:{
            last_frame=(int)(10*frame_rate);
            solids_parameters.cfl=(T)10.0;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            initial_orientation=ROTATION<TV>((T)pi/2,TV(0,0,1));
            initial_velocity=TV();
            initial_angular_velocity=TV();
            roller_speed=(T)2;
            roller_friction=(T).1;
            roller_orientation=ROTATION<TV>(0,TV(1,0,0));
            cylinder_time=(T).2;
            cylinder_start=TV(0,2,0);
            cylinder_velocity=TV(0,-5,0);
            RIGID_BODY<TV>&gear1=tests.Add_Rigid_Body("gear",(T).375,roller_friction);
            gear1.Frame().t=TV(-(T).4,(T)1.5,(T)-.75);
            gear1.Frame().r=roller_orientation;
            gear1.Twist().angular=-roller_speed*TV(0,0,1);

            RIGID_BODY<TV>&gear2=tests.Add_Rigid_Body("gear",(T).375,roller_friction);
            gear2.Frame().t=TV((T).4,(T)1.5,(T)-.75);
            gear2.Frame().r=roller_orientation;
            gear2.Twist().angular=roller_speed*TV(0,0,1);

            RIGID_BODY<TV>&cylinder=tests.Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).375/2,0);
            cylinder.Frame().t=cylinder_start;
            cylinder.Twist().linear=cylinder_velocity;
            cylinder.Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));

            TETRAHEDRALIZED_VOLUME<T>& buddha=tests.Create_Tetrahedralized_Volume(buddha_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,false,1000);
            buddha.Update_Bounding_Box();
            TV center(buddha.bounding_box->Center());
            for(int i=0;i<particles.Size();i++){
                particles.V(i)=initial_velocity+TV::Cross_Product(initial_angular_velocity,particles.X(i)-center);
                particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center)*buddha_resize;
                particles.X(i).y+=(T)2.234;}
//            use_soft_bindings=true;
//            rigid_body_soft=true;
            tests.Add_Ground((T)1);

            use_gravity=true;
            youngs_modulus=(T)50000*stiffen;
            break;}
        case 8:{
            last_frame=(int)(5*frame_rate);
            solids_parameters.cfl=(T)10.0;
            RIGID_BODY<TV>& box=tests.Add_Rigid_Body("cutout_box",(T)0.708855526,0);
            box.Frame().r=ROTATION<TV>(-(T)pi/2,TV(1,0,0));

            TETRAHEDRALIZED_VOLUME<T>& buddha=tests.Create_Tetrahedralized_Volume(buddha_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,false,1000);
            buddha.Update_Bounding_Box();
            TV center(buddha.bounding_box->Center());
            for(int i=0;i<particles.Size();i++){
                particles.X(i).y+=(T)2.234;}

            use_gravity=true;
            youngs_modulus=(T)1000*stiffen;
            break;}
        case 9:{
            last_frame=(int)(2*frame_rate);
            RIGID_BODY<TV>& box=tests.Add_Rigid_Body("cutout_box",(T)1.611991954,0);
            box.Frame().r=ROTATION<TV>(-(T)pi/2,TV(1,0,0));
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1.2,0))),true,false,1000);

            use_gravity=true;
            youngs_modulus=(T)10000*stiffen;
            break;}
        case 10:{
            last_frame=(int)(5*frame_rate);
            solids_parameters.cfl=(T)2;
            initial_orientation=ROTATION<TV>((T)pi/2,TV(0,0,1));
            tests.Create_Tetrahedralized_Volume(armadillo_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.234,0)),TWIST<TV>(TV(),TV(0,0,1))),true,false,1000,armadillo_scale);
            tests.Add_Ground(ground_friction);

            use_gravity=true;
            youngs_modulus=(T)4000*stiffen;
            hardening_deformation=1;
            hardening_strength=1;  // 1 or 2
            binding_stiffness=(T)1e4;
            break;}
        case 11:{
            last_frame=(int)(10*frame_rate);
            solids_parameters.cfl=(T)10.0;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            GRID<TV> mattress_grid(TV_INT(5,10,5),RANGE<TV>(TV((T)-.25,(T).8,(T)-.30),TV((T).25,(T)5.10,(T).30)));
            tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(tetrahedralized_volume,1000,true);
            tests.Create_Mattress(mattress_grid);
            solids_parameters.deformable_object_collision_parameters.collide_with_interior=false;
            tests.Add_Ground((T).3);

            use_gravity=true;
            break;}
        case 12:{
            last_frame=(int)(8*frame_rate);
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)1.5,1,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-.5,0,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-2.5,-1,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-4.5,-2,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-6.5,-3,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-8.5,-4,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-10.5,-5,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-12.5,-6,0);body->is_static=true;
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,3,0))),true,false,1000);

            use_gravity=true;
            use_neohookean=true;
            break;}
        case 13:{
            last_frame=(int)(5*frame_rate);
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,false,1000);
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)5.1,0))),true,false,1000);
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)7.2,0))),true,false,1000);
            TETRAHEDRALIZED_VOLUME<T>& merged_surface=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            for(int s=0;s<deformable_body_collection.structures.m;s++){TETRAHEDRALIZED_VOLUME<T>*ts=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(s));
                if(ts) merged_surface.mesh.elements.Append_Elements(ts->mesh.elements);}
            deformable_body_collection.structures.Clean_Memory();deformable_body_collection.Add_Structure(&merged_surface);
            tests.Add_Ground(0);

            use_gravity=true;
            use_neohookean=true;
            hardening_deformation=(T).01;
            hardening_strength=(T).25;
            break;}
        case 14:{
            last_frame=(int)(10*frame_rate);
            solids_parameters.cfl=(T)10.0;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            initial_orientation=ROTATION<TV>((T)pi/2,TV(0,0,1))*ROTATION<TV>(-(T)pi/2,TV(0,1,0));
            initial_velocity=TV();
            initial_angular_velocity=TV();
            roller_speed=(T)2;
            roller_friction=(T).1;
            roller_orientation=ROTATION<TV>(0,TV(1,0,0));
            cylinder_time=(T).2;
            cylinder_start=TV(0,2,0);
            cylinder_velocity=TV(0,-5,0);
            RIGID_BODY<TV>&gear1=tests.Add_Rigid_Body("gear",(T).375,roller_friction);
            gear1.Frame().t=TV(-(T).4,(T)1.5,(T)-.75);
            gear1.Frame().r=roller_orientation;
            gear1.Twist().angular=-roller_speed*TV(0,0,1);

            RIGID_BODY<TV>&gear2=tests.Add_Rigid_Body("gear",(T).375,roller_friction);
            gear2.Frame().t=TV((T).4,(T)1.5,(T)-.75);
            gear2.Frame().r=roller_orientation;
            gear2.Twist().angular=roller_speed*TV(0,0,1);

            RIGID_BODY<TV>&cylinder=tests.Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).375/2,0);
            cylinder.Frame().t=cylinder_start;
            cylinder.Twist().linear=cylinder_velocity;
            cylinder.Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));
            
            tests.Create_Tetrahedralized_Volume(armadillo_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.3,0),initial_orientation)),true,false,1000,armadillo_scale);
            tests.Add_Ground((T)1);

            use_gravity=true;
            youngs_modulus=(T)1e6*stiffen;
            break;}
        case 15:{
            last_frame=(int)(10*frame_rate);
            initial_orientation=ROTATION<TV>((T)pi/2,TV(0,0,1))*ROTATION<TV>(-(T)pi/2,TV(0,1,0));
            initial_velocity=TV();
            initial_angular_velocity=TV();
            roller_speed=(T)2;
            roller_friction=(T).1;
            roller_orientation=ROTATION<TV>(0,TV(1,0,0));
            cylinder_time=(T).2;
            cylinder_start=TV(0,2,0);
            cylinder_velocity=TV(0,-5,0);
            RIGID_BODY<TV>&gear1=tests.Add_Rigid_Body("gear",(T).375,roller_friction);
            gear1.Frame().t=TV(-(T).4,(T)1.5,(T)-.75);
            gear1.Frame().r=roller_orientation;
            gear1.Twist().angular=-roller_speed*TV(0,0,1);

            RIGID_BODY<TV>&gear2=tests.Add_Rigid_Body("gear",(T).375,roller_friction);
            gear2.Frame().t=TV((T).4,(T)1.5,(T)-.75);
            gear2.Frame().r=roller_orientation;
            gear2.Twist().angular=roller_speed*TV(0,0,1);

            RIGID_BODY<TV>&cylinder=tests.Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).375/2,0);
            cylinder.Frame().t=cylinder_start;
            cylinder.Twist().linear=cylinder_velocity;
            cylinder.Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));
            
            tests.Create_Tetrahedralized_Volume(sphere_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.3,0),initial_orientation)),true,false,1000,(T).3);
            tests.Add_Ground((T)1);

            use_gravity=true;
            youngs_modulus=(T)2e4*stiffen;
            break;}
        case 16:{
            last_frame=(int)(8*frame_rate);
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)1.5,1,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-.5,0,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-2.5,-1,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-4.5,-2,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-6.5,-3,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-8.5,-4,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-10.5,-5,0);body->is_static=true;
            body=&tests.Add_Rigid_Body("plank",(T)1,(T)0);body->Frame().t=TV((T)-12.5,-6,0);body->is_static=true;
            tests.Create_Tetrahedralized_Volume(armadillo_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T).33,(T)3,0),initial_orientation)),true,false,1000,armadillo_scale*4);
            tests.Add_Ground(0,(T)-7);

            use_gravity=true;
            break;}
        case 17:
        case 18:{
            last_frame=(int)(10*frame_rate);
            solids_parameters.cfl=(T)10;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            if(test_number==18) max_cg_iterations=200;
            T size=(T)3.01;
            GRID<TV> grid(TV_INT(tori_m,tori_n,tori_mn),RANGE<TV>(TV(0,tori_base_height,0),TV((tori_m-1)*size,tori_base_height+(tori_n-1)*size,(tori_mn-1)*size)));
            RANDOM_NUMBERS<T> random;random.Set_Seed(12321);
            tori_initial_states.Remove_All(); 
            for(int j=0;j<tori_n;j++)for(int i=0;i<tori_m;i++)for(int ij=0;ij<tori_mn;ij++){ // Note nonstandard dimension ordering to get nice layers
                ROTATION<TV> orientation(random.Get_Uniform_Number((T)0,(T)1),random.template Get_Direction<TV>());
                TV new_center=grid.Node(VECTOR<int,3>(i,j,ij));
                TV angular_velocity=random.Get_Uniform_Number((T)0,tori_max_angular_velocity)*random.template Get_Direction<TV>();
                std::cout<<"Adding torus at center "<<new_center<<", orientation "<<orientation<<", angular velocity "<<angular_velocity<<"\n";
                std::cout<<i<<" "<<j<<" "<<ij<<" "<<grid.dX.y<<"\n";
                tori_initial_states.Append(RIGID_BODY_STATE<TV>(FRAME<TV>(new_center,orientation),TWIST<TV>(TV(),angular_velocity)));
                tests.Create_Tetrahedralized_Volume(torus_filename,tori_initial_states.Last(),true,false,1000,1);}
            TETRAHEDRALIZED_VOLUME<T>& merged_surface=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            for(int s=0;s<deformable_body_collection.structures.m;s++){TETRAHEDRALIZED_VOLUME<T>*ts=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(s));
                if(ts) merged_surface.mesh.elements.Append_Elements(ts->mesh.elements);}
            deformable_body_collection.structures.Clean_Memory();deformable_body_collection.Add_Structure(&merged_surface);
            // rigid bodies
            if(test_number==18){
                for(int i=0;i<16;i++){T a=i*(T)pi/8;
                    body=&tests.Add_Rigid_Body("skinnyhexlink",(T)2,(T)0);body->Frame().t=TV((T)cos(a)*6+(T)1.5,(T)2,(T)sin(a)*6+(T)1.5);body->is_static=true;}
                for(int i=0;i<32;i++){T a=i*(T)pi/16;
                    body=&tests.Add_Rigid_Body("skinnyhexlink",(T)2,(T)0);body->Frame().t=TV((T)cos(a)*12+(T)1.5,(T)2,(T)sin(a)*12+(T)1.5);body->is_static=true;}}
            tests.Add_Ground((T).1);

            use_gravity=true;
            break;}
        case 19:{
            last_frame=(int)(5*frame_rate);
            tests.Create_Tetrahedralized_Volume(torus_filename,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,false,1000);
            tests.Add_Ground(0);
            use_gravity=true;
            use_neohookean=true;
            hardening_deformation=(T).01;
            hardening_strength=(T).25;
            break;}
        case 20:{
            frame_rate=24;last_frame=(int)(10*frame_rate);
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            int number_side_panels=50;T aspect_ratio=(T)1;
            TRIANGULATED_SURFACE<T>& surface=tests.Create_Cloth_Panel(number_side_panels,(T)1,aspect_ratio,RIGID_BODY_STATE<TV>());
            T density=TV::dimension==1?1:TV::dimension==2?100:1000;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,density);
            int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
            i=1;j=1;//particles.mass(i+m*(j-1))=FLT_MAX;
            fixed_particles.Append(i+m*(j-1));fixed_particle_positions.Append(particles.X(i+m*(j-1)));
            i=1;j=n;//particles.mass(i+m*(j-1))=FLT_MAX;
            fixed_particles.Append(i+m*(j-1));fixed_particle_positions.Append(particles.X(i+m*(j-1)));
            use_gravity=true;
            use_neohookean=false;
            youngs_modulus=(T)2e4*stiffen;
            hardening_deformation=(T).05;
            hardening_strength=(T).25;
            break;}
         default:
            LOG::cerr<<"(3) Unrecognized test number "<<test_number<<std::endl;exit(1);}
    
    bool automatically_add_to_collision_structures=true;
    solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
    solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
    if(use_soft_bindings && solids_parameters.triangle_collision_parameters.perform_self_collision){
        solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
        solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
        automatically_add_to_collision_structures=false;
        for(int i=0;i<deformable_body_collection.structures.m;i++) if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(i)))
            Collide_With_Soft_Bound_Surface_Copy(*volume,rigid_body_soft);}

    if(use_volumetric_self_collisions && solids_parameters.triangle_collision_parameters.perform_self_collision){
        PHYSBAM_FATAL_ERROR();
        solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
        solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
        solids_parameters.triangle_collision_parameters.turn_off_all_collisions=true;
        solid_body_collection.deformable_body_collection.triangle_repulsions.compute_edge_edge_friction=false;
        solid_body_collection.deformable_body_collection.triangle_repulsions.compute_edge_edge_inelastic_collision_repulsion=false;
        solid_body_collection.deformable_body_collection.triangle_repulsions.compute_edge_edge_repulsion=false;

        solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1e-2;
        solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).25;
        solids_parameters.triangle_collision_parameters.collisions_collision_thickness=(T)1e-5;
        solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=8;}

    // add structures and rigid bodies to collisions
    if(automatically_add_to_collision_structures){
        deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
        solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    LOG::Stat("self_collide_with_interior_nodes",self_collide_with_interior_nodes);

    if(solids_parameters.triangle_collision_parameters.perform_self_collision){
        if(use_tet_collisions && solid_body_collection.deformable_body_collection.mpi_solids){LOG::cout<<"use_tet_collisions and MPI are incompatible"<<std::endl;PHYSBAM_FATAL_ERROR();}
        if(use_soft_bindings+use_tet_collisions+use_volumetric_self_collisions>1) PHYSBAM_FATAL_ERROR();
        if(use_soft_bindings){
            soft_bindings.use_impulses_for_collisions.Fill(false);
            soft_bindings.Initialize_Binding_Mesh();
            LOG::Stat("binding stiffness",binding_stiffness);
            solid_body_collection.Add_Force(Create_Edge_Binding_Springs(particles,*soft_bindings.binding_mesh,binding_stiffness,(T)1));}
        else if(use_tet_collisions){
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            for(int s=0;s<deformable_body_collection.structures.m;s++) if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(s)))
                tests.Initialize_Tetrahedron_Collisions(s,output_directory,*volume,solids_parameters.triangle_collision_parameters);}
        else if(self_collide_with_interior_nodes){
            FREE_PARTICLES<TV>& interior=*FREE_PARTICLES<TV>::Create(particles);
            for(int s=0;s<deformable_body_collection.structures.m;s++) if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(s))){
                if(!volume->mesh.node_on_boundary) volume->mesh.Initialize_Node_On_Boundary();
                ARRAY<int> nodes;volume->mesh.elements.Flattened().Get_Unique(nodes);
                for(int i=0;i<nodes.m;i++)if(!(*volume->mesh.node_on_boundary)(nodes(i))) interior.nodes.Append(nodes(i));}
            LOG::Stat("interior particles",interior.nodes.m);
            solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(&interior);}}

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        if(power_of_two(solid_body_collection.deformable_body_collection.mpi_solids->number_of_ranks))
            solid_body_collection.deformable_body_collection.mpi_solids->KD_Tree_Partition(solid_body_collection.deformable_body_collection,solid_body_collection.rigid_body_collection,ARRAY<TV>(particles.X));
        else{LOG::cout<<"needs 2 or 4 procs"<<std::endl;PHYSBAM_FATAL_ERROR();}}
    
    LOG::cout<<"number of particles: "<<deformable_body_collection.particles.Size()<<std::endl;
    if(use_gravity) solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
    {LOG::SCOPE scope("creating finite volume","creating finite volume");
    LOG::Stat("use_neohookean",use_neohookean);
    LOG::Stat("youngs_modulus",youngs_modulus);
    LOG::Stat("incompressible",test_poissons_ratio==(T).5);
    LOG::Stat("internal_poissons_ratio",internal_poissons_ratio);
    LOG::Stat("hardening_deformation",hardening_deformation);
    LOG::Stat("hardening_strength",hardening_strength);
    LOG::Stat("merge_at_boundary",merge_at_boundary);
    LOG::Stat("use_neumann",use_neumann);
    LOG::Stat("minimum_volume_recovery_time_scale",minimum_volume_recovery_time_scale);
    LOG::Stat("max_cg_iterations",max_cg_iterations);
    for(int i=0;i<deformable_body_collection.structures.m;i++){
        TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(i));
        TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.structures(i));
        if(volume){LOG::SCOPE scope("mesh statistics","mesh statistics");volume->Print_Statistics(LOG::cout);}
        Add_Incompressible_Force(volume);
        Add_Incompressible_Force(surface);}}
    solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
// Function Add_Incompressible_Force
//#####################################################################
template<class T_OBJECT>
void Add_Incompressible_Force(T_OBJECT* object)
{
    if(!object) return;
    enum WORKAROUND {d=T_OBJECT::BASE::MESH_OBJECT::MESH::dimension};
    CONSTITUTIVE_MODEL<T,d>* constitutive_model;
    if(use_neohookean) constitutive_model=new NEO_HOOKEAN<T,d>(youngs_modulus,internal_poissons_ratio,hardening_deformation,hardening_strength);
    else constitutive_model=new SPLINE_MODEL<T,d>(youngs_modulus,internal_poissons_ratio,hardening_deformation,hardening_strength);
    INCOMPRESSIBLE_FINITE_VOLUME<TV,d>& fvm=*Create_Incompressible_Finite_Volume(*object);
    solid_body_collection.Add_Force(&fvm);
    solid_body_collection.Add_Force(Create_Finite_Volume(*object,constitutive_model));
    fvm.mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    fvm.merge_at_boundary=merge_at_boundary;
    fvm.minimum_volume_recovery_time_scale=minimum_volume_recovery_time_scale;
    fvm.max_cg_iterations=max_cg_iterations;
    fvm.use_neumann=use_neumann;
    if(test_poissons_ratio<(T).5) fvm.disable_projection=true;
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    if(test_number==17 || test_number==18){
        LOG::SCOPE scope("restart hack","restart hack");
        ARRAY<TV> X_save(particles.X);
        ARRAY<T> mass_save(particles.mass);
        FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/deformable_body_collection_particles."+FILE_UTILITIES::Number_To_String(frame),particles);
        TETRAHEDRALIZED_VOLUME<T>& volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        int particles_per_torus=volume.mesh.number_nodes/tori_initial_states.m;
        if(volume.mesh.number_nodes%particles_per_torus!=0 || particles.Size()%particles_per_torus!=0) PHYSBAM_FATAL_ERROR();
        int number_of_tori_read=particles.Size()/particles_per_torus;
        int number_of_old_tori=number_of_tori_read;
        if(number_of_tori_read>tori_initial_states.m) PHYSBAM_FATAL_ERROR();
        particles.Add_Elements(particles_per_torus*(tori_initial_states.m-number_of_tori_read));
        T time=Time_At_Frame(frame);
        TV fall_X(0,-(T).5*(T)9.8*sqr(time),0);
        TV fall_V(0,-(T)9.8*time,0);
        LOG::cout<<"fall_X = "<<fall_X<<std::endl;
        LOG::cout<<"fall_V = "<<fall_V<<std::endl;
        LOG::Stat("old tori",number_of_old_tori);
        LOG::Stat("new tori",tori_initial_states.m-number_of_old_tori);
        LOG::Stat("total tori",tori_initial_states.m);
        for(int torus=number_of_old_tori+1;torus<=tori_initial_states.m;torus++){
            RIGID_BODY_STATE<TV>& state=tori_initial_states(torus);
            for(int i=0;i<particles_per_torus;i++){int p=i+particles_per_torus*(torus-1);
                particles.X(p)=state.frame.t+fall_X+ROTATION<TV>::From_Rotation_Vector(time*state.twist.angular).Rotate(X_save(p)-state.frame.t);
                particles.V(p)=fall_V+TV::Cross_Product(state.twist.angular,X_save(p)-state.frame.t);
                particles.mass(p)=mass_save(p);}}
        if(solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Check_For_Intersection(false,
                solids_parameters.triangle_collision_parameters.collisions_small_number))
            throw std::runtime_error("Tried to start with intersections Found");}
    else{
        BASE::Read_Output_Files_Solids(frame);
        solid_body_collection.Update_Simulated_Particles();}
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    if(BINDING_SPRINGS<TV>* binding_springs=solid_body_collection.template Find_Force<BINDING_SPRINGS<TV>*>()){
        T critical_distance=(T).01;
        T base_stiffness=binding_springs->constant_stiffness;
        ARRAY<T>& stiffness=binding_springs->stiffness;stiffness.Resize(soft_bindings.bindings.m);
        for(int a=0;a<soft_bindings.bindings.m;a++){VECTOR<int,2> binding=soft_bindings.bindings(a);
            T distance=(particles.X(binding.x)-particles.X(binding.y)).Magnitude();
            T ratio=distance/critical_distance;
            stiffness(a)=ratio>1?base_stiffness/ratio:base_stiffness;}
        binding_springs->Set_Overdamping_Fraction((T)1);
        LOG::cout<<"clamping binding stiffness: base = "<<base_stiffness<<", min = "<<stiffness.Min()<<std::endl;}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if((test_number==1 || test_number==2 || test_number==3 || test_number==4 || test_number==6) && id==int(1)) twist.linear=curve.Derivative(time);
    else if((test_number==1 || test_number==2 || test_number==3 || test_number==4 || test_number==6) && id==int(2)) twist.linear=curve2.Derivative(time);
    else if((test_number==7 || test_number==14 || test_number==15) && id==int(3)) twist.linear=time>cylinder_time?TV():cylinder_velocity;
    else return false;
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if((test_number==1 || test_number==2 || test_number==3 || test_number==4 || test_number==6) && id==int(1)) frame.t=curve.Value(time);
    else if((test_number==1 || test_number==2 || test_number==3 || test_number==4 || test_number==6) && id==int(2)) frame.t=curve2.Value(time);
    else if((test_number==7 || test_number==14 || test_number==15) && id==int(1)) frame.r=ROTATION<TV>(-roller_speed*time,TV(0,0,1))*roller_orientation;
    else if((test_number==7 || test_number==14 || test_number==15) && id==int(2)) frame.r=ROTATION<TV>(roller_speed*time,TV(0,0,1))*roller_orientation;
    else if((test_number==7 || test_number==14 || test_number==15) && id==int(3)) frame.t=time>cylinder_time?cylinder_start:cylinder_start+(time-cylinder_time)*cylinder_velocity;
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<fixed_particles.m;i++) V(fixed_particles(i))=TV();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<fixed_particles.m;i++) V(fixed_particles(i))=TV();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<fixed_particles.m;i++) X(fixed_particles(i))=fixed_particle_positions(i);
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<fixed_particles.m;i++) X(fixed_particles(i))=TV();
}
//#####################################################################
};
}
#endif

