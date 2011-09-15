//#####################################################################
// Copyright 2007-2008, Jon Gretarsson, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_WATER
//#####################################################################
//   1. Falling deformable sphere and smoke
//   2. Heavy rigid objects falling into water
//   3. Heavy block falling and compressing in column of water
//   4. Many rigid bodies falling into a pool with a fountain
//   6. thin shell rigid boat hit with one light sphere and one heavy sphere
//   7. cloth hit with a heavy rigid body
//   8. balloon filling with water
//   9. Fish flopping in water
//  10. buoyant balls floating in water
//  11. submerged bag being pulled out of water
//  12. sphere splashing into a pool of water
//  13. domain filling with water; many rigid objects being tossed around
//  14. cloth hit with heavy rigid sphere
//  15. Falling deformable torus in water
//  16. case 7 with no coupling
//  17. case 7, lighter
//  18. flag partially submerged in water
//  19. light flag falling into water
//  20. pulling cloth out of water
//  21. Latice of rigid bodies falling into a domain of water
//#####################################################################
#ifndef __STANDARD_TESTS_WATER__
#define __STANDARD_TESTS_WATER__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_WATER:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef GRID<TV> T_GRID;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::stream_type;using BASE::data_directory;using BASE::solid_body_collection;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::parse_args;using BASE::test_number;using BASE::mpi_world;using BASE::resolution;

    WATER_STANDARD_TESTS_3D<T_GRID > water_tests;
    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    T_GRID mattress_grid;
    int deformable_object_id;
    T solid_density;
    T stiffness_ratio;
    int boat_index,light_sphere_index,heavy_sphere_index;
    T light_sphere_initial_height,heavy_sphere_initial_height;
    T light_sphere_drop_time,heavy_sphere_drop_time;
    T balloon_initial_radius;
    T initial_fluid_height;
    T boat_mass;
    bool implicit_springs;
    ARRAY<CYLINDER<T> > fountain_source;
    ARRAY<BOX<TV> > fountain_source_boxes;
    ARRAY<TV> fountain_source_velocity;
    MATRIX<T,4> world_to_source;
    int bodies;
    int sub_test;
    TV left_corner,right_corner;
    INTERPOLATION_CURVE<T,T> bag_curve;

    INTERPOLATION_CURVE<T,TV> curve;

    ARRAY<PAIR<int,TV> > constrained_node_positions;
    ARRAY<TRIPLE<int,T,TV> > constrained_nodes;
    ARRAY<int> rigid_bodies_to_simulate;
    ARRAY<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*> deformable_objects_to_simulate;
    ARRAY<int> rigid_bodies_to_collide_against;
    ORIENTED_BOX<TV> fish_bounding_box;
    LEVELSET_IMPLICIT_OBJECT<TV>* fish_levelset;
    bool solid_node;
    bool mpi;

    STANDARD_TESTS_WATER(const STREAM_TYPE stream_type)
        :BASE(stream_type,solid_node?0:1,fluids_parameters.WATER),
        water_tests(*this,fluids_parameters,solid_body_collection.rigid_body_collection),
        solids_tests(*this,solid_body_collection),deformable_object_id(0),solid_density((T)2000),stiffness_ratio(0),light_sphere_index(0),heavy_sphere_index(0),
        light_sphere_initial_height((T)1.75),heavy_sphere_initial_height((T)1.75),light_sphere_drop_time((T)1),heavy_sphere_drop_time((T)1.5),balloon_initial_radius((T)0),
        initial_fluid_height((T)0),boat_mass((T)7),implicit_springs(0),world_to_source(MATRIX<T,4>::Identity_Matrix()),bodies(5),sub_test(0),fish_levelset(0)
    {
        solid_node=mpi_world->initialized && !mpi_world->rank;
        mpi=mpi_world->initialized;
    }

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(T_ARRAYS_BOOL*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-suboption",1);
    parse_args->Add_Double_Argument("-stiffness",(T)1);
    parse_args->Add_Option_Argument("-implicit");

    parse_args->Add_Integer_Argument("-bodies",5);
    parse_args->Add_Double_Argument("-mass",(T)1);
    parse_args->Add_Integer_Argument("-xprocs",0);
    parse_args->Add_Integer_Argument("-zprocs",0);
    parse_args->Add_Integer_Argument("-iterations",200);

    parse_args->Add_Double_Argument("-height",1.5);
    parse_args->Add_Double_Argument("-stop_time",1);
    parse_args->Add_Double_Argument("-radius",.2);

}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    water_tests.Initialize(Water_Test_Number(test_number),resolution);
    sub_test=parse_args->Get_Integer_Value("-suboption");
    stiffness_ratio=(T)parse_args->Get_Double_Value("-stiffness");
    implicit_springs=parse_args->Get_Option_Value("-implicit");

    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    last_frame=1000;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;

    fluids_parameters.solve_neumann_regions=false;
    fluids_parameters.incompressible_iterations=200;
    *fluids_parameters.grid=water_tests.grid;
    fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;
    fluids_parameters.second_order_cut_cell_method=true;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.verbose_dt=true;
    solid_body_collection.print_residuals=false;
    solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1;

    if(solid_node || !mpi) solids_parameters.use_rigid_deformable_contact=true;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests_Water/Test_%d_Resolution_%d_Stiffness_%d_Suboption_%d",test_number,resolution,stiffness_ratio,sub_test);
    if(parse_args->Is_Value_Set("-mass")){boat_mass=(T)parse_args->Get_Double_Value("-mass");output_directory+=STRING_UTILITIES::string_sprintf("_mass%f",boat_mass);solid_density=boat_mass;}
    if(parse_args->Is_Value_Set("-bodies")){bodies=parse_args->Get_Integer_Value("-bodies");output_directory+=STRING_UTILITIES::string_sprintf("_bodies%f",bodies);}
        
    T spout_stop_time=0,ball_initial_height=0,spout_radius=0;
    if(parse_args->Is_Value_Set("-height")){ball_initial_height=(T)parse_args->Get_Double_Value("-height");output_directory+=STRING_UTILITIES::string_sprintf("_height%.02f",ball_initial_height);}
    if(parse_args->Is_Value_Set("-stop_time")){spout_stop_time=(T)parse_args->Get_Double_Value("-stop_time");output_directory+=STRING_UTILITIES::string_sprintf("_stopat%.02f",spout_stop_time);}
    if(parse_args->Is_Value_Set("-radius")){spout_radius=(T)parse_args->Get_Double_Value("-radius");output_directory+=STRING_UTILITIES::string_sprintf("_radius%.02f",spout_radius);}

    fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[2][1]=false;
    fluids_parameters.density=(T)1000;
    solids_parameters.implicit_solve_parameters.cg_iterations=400;

    switch(test_number){
        case 15:
        case 1:
            last_frame=240;
            fluids_parameters.density=(T)1000;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            //mattress_grid=T_GRID(sub_test*20,sub_test*5,sub_test*20,(T).55,(T).95,(T).4,(T).52,(T).55,(T).95);
            mattress_grid=T_GRID(sub_test*20,sub_test*5,sub_test*20,(T).3,(T).7,(T).6,(T).72,(T).3,(T).7);
            (*fluids_parameters.grid).Initialize(10*resolution+1,10*resolution+1,10*resolution+1,(T)0,(T)1,(T)0,(T)1,(T)0,(T)1);
            //solids_parameters.implicit_solve_parameters.cg_iterations=1000;
            if(!parse_args->Is_Value_Set("-iterations")) solids_parameters.implicit_solve_parameters.cg_iterations=500;
            solid_density=(T)300;
            break;
        case 2:
            fluids_parameters.density=(T)10;
            last_frame=100;
            fluids_parameters.reseeding_frame_rate=5;
            fluids_parameters.domain_walls[2][2]=false;
            /* 
               case 1: heavy ball into bowl
               case 2: light ball into bowl
            */
            (*fluids_parameters.grid).Initialize(20*resolution+1,30*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)3,(T)-1,(T)1);
            left_corner=TV((T)-.5,(T)0,(T)-.5);
            right_corner=TV((T).5,(T)1,(T).5);
            break;
        case 3:
            fluids_parameters.density=(T)1000;
            solid_density=(T)30000;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            (*fluids_parameters.grid).Initialize(2*resolution+1,2*resolution+1,10*resolution+1,(T)0,(T).2,(T)0,(T).2,(T)0,(T)1);
            mattress_grid=T_GRID(20,20,20,(T).05,(T).15,(T).05,(T).15,(T).9,(T)1.05);
            //mattress_grid=T_GRID(20,20,(T).05,(T).15,(T).8,(T).95);
            break;
        case 4:
            solids_parameters.rigid_body_collision_parameters.use_push_out=true;
            //solids_parameters.enforce_rigid_rigid_contact_in_cg=true;
            if(solid_node) solids_parameters.use_rigid_deformable_contact=true;
            fluids_parameters.density=(T)1000;
            solid_density=(T)1000;
            fountain_source.Append(CYLINDER<T>(TV((T).7,(T).225,(T).35),TV((T).7,(T).25,(T).35),(T).05));
            fountain_source.Append(CYLINDER<T>(TV((T).48,(T).225,(T).48),TV((T).48,(T).25,(T).48),(T).04));
            fountain_source.Append(CYLINDER<T>(TV((T).35,(T).225,(T).5),TV((T).35,(T).25,(T).5),(T).05));
            fountain_source.Append(CYLINDER<T>(TV((T).6,(T).225,(T).6),TV((T).6,(T).25,(T).6),(T).045));
            fountain_source.Append(CYLINDER<T>(TV((T).375,(T).225,(T).325),TV((T).6,(T).25,(T).6),(T).055));
            fountain_source_velocity.Append(TV((T)0,(T)3,(T)0));
            fountain_source_velocity.Append(TV((T)0,(T)3.5,(T)0));
            fountain_source_velocity.Append(TV((T)0,(T)3.2,(T)0));
            fountain_source_velocity.Append(TV((T)0,(T)3,(T)0));
            fountain_source_velocity.Append(TV((T)0,(T)3.3,(T)0));
            fluids_parameters.reseeding_frame_rate=10;
            (*fluids_parameters.grid).Initialize(10*resolution+1,15*resolution+1,10*resolution+1,(T)0,(T)1,(T)0,(T)1.5,(T)0,(T)1);
            break;
        case 5:
            fluids_parameters.density=(T)1000;
            solid_density=(T)100000;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            (*fluids_parameters.grid).Initialize(10*resolution+1,20*resolution+1,10*resolution+1,(T)0,(T)1,(T)0,(T)2,(T)0,(T)1);
            break;
        case 6:
            fluids_parameters.density=(T)1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            if(solid_node) solids_parameters.use_rigid_deformable_contact=true;
            solids_parameters.rigid_body_collision_parameters.use_push_out=true;
            solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T).1;
            solid_density=(T)3000;
            fluids_parameters.grid->Initialize(20*resolution+1,15*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)1.5,(T)-1,(T)1);
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T).6*(T)1.5;;
            break;
        case 7:
            //if(solid_node) solids_parameters.use_rigid_deformable_contact=true;
            solids_parameters.use_rigid_deformable_contact=false;
            if(sub_test!=0) fluids_parameters.second_order_cut_cell_method=true; else fluids_parameters.second_order_cut_cell_method=false;
            solids_parameters.triangle_collision_parameters.allow_intersections=true;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            fluids_parameters.density=(T)1000;
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=0;
            solids_parameters.implicit_solve_parameters.cg_iterations=500;
            solid_density=(T)8000;
            fluids_parameters.grid->Initialize(20*resolution+1,25*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)2.5,(T)-1,(T)1);
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T).45*(T)2.5;
            break;
        case 16:
            //if(solid_node) solids_parameters.use_rigid_deformable_contact=true;
            fluids_parameters.simulate=false;
            solids_parameters.use_rigid_deformable_contact=false;
            fluids_parameters.fluid_affects_solid=false;fluids_parameters.solid_affects_fluid=false;
            if(sub_test!=0) fluids_parameters.second_order_cut_cell_method=true; else fluids_parameters.second_order_cut_cell_method=false;
            solids_parameters.triangle_collision_parameters.allow_intersections=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            fluids_parameters.number_of_regions=0;
            fluids_parameters.density=(T)1000;
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=0;
            solids_parameters.implicit_solve_parameters.cg_iterations=500;
            solid_density=(T)8000;
            fluids_parameters.grid->Initialize(20*resolution+1,25*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)2.5,(T)-1,(T)1);
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T).45*(T)2.5;
            break;
        case 8:
            fluids_parameters.reincorporate_removed_particle_velocity=true;
            solids_parameters.use_rigid_deformable_contact=false;
            fluids_parameters.density=(T)1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            fluids_parameters.reseeding_frame_rate=5;
            fluids_parameters.grid->Initialize(20*resolution+1,25*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)2.5,(T)-1,(T)1);
            fluids_parameters.domain_walls[2][2]=false;
            light_sphere_initial_height=(T)1.2;
            heavy_sphere_initial_height=(T)2.9;  // height of the funnel
            initial_fluid_height=(T)0;balloon_initial_radius=(T).325;

            // TODO Choose this parameter once drop time and spout size are chosen
            heavy_sphere_drop_time=(T)10; // this should be soon after steady state; a variable dependent on other parameter choices

            fountain_source.Append(CYLINDER<T>(TV((T)0,heavy_sphere_initial_height-(T).7,(T)0),TV((T)0,heavy_sphere_initial_height-(T).3,(T)0),spout_radius));
            fountain_source_velocity.Append(TV((T)0,ball_initial_height,(T)0));
            light_sphere_drop_time=spout_stop_time;
            break;
        case 9:
            fluids_parameters.reincorporate_removed_particle_velocity=true;
            fluids_parameters.removed_particle_mass_scaling=60;
            fluids_parameters.density=(T)1000;
            if(!parse_args->Is_Value_Set("-iterations")) solids_parameters.implicit_solve_parameters.cg_iterations=800;
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T)3.5;
            fluids_parameters.grid->Initialize(32*resolution+1,36*resolution+1,24*resolution+1,(T)-8,(T)8,(T)0,(T)18,(T)-6,(T)6);
            PHYSBAM_ASSERT(fluids_parameters.grid->dX.x==fluids_parameters.grid->dX.y && fluids_parameters.grid->dX.y==fluids_parameters.grid->dX.z);
            break;
        case 10:
            last_frame=500;
            fluids_parameters.density=(T)1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            fluids_parameters.grid->Initialize(40*resolution+1,15*resolution+1,10*resolution+1,(T)-2,(T)2,(T)0,(T)1.5,(T)-.5,(T).5);
            solid_density=(T)1; for(int i=1;i<sub_test;++i) solid_density*=(T)10;
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T)1.2;
            light_sphere_initial_height=ball_initial_height;
            break;
        case 11:
            fluids_parameters.reincorporate_removed_particle_velocity=true;
            fluids_parameters.removed_particle_mass_scaling=60;
            fluids_parameters.density=(T)1000;
            if(!parse_args->Is_Value_Set("-iterations")) solids_parameters.implicit_solve_parameters.cg_iterations=50;
            fluids_parameters.reseeding_frame_rate=5;
            fluids_parameters.grid->Initialize(20*resolution+1,30*resolution+1,20*resolution+1,(T)-1.5,(T)1.5,(T)0,(T)4,(T)-1.5,(T)1.5);
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T)1.35;
            balloon_initial_radius=(T).325;
            heavy_sphere_drop_time=(T)2.5;
            light_sphere_initial_height=(T).75;
            bag_curve.Add_Control_Point(0,(T)0);bag_curve.Add_Control_Point((T)1.79,(T)1.3425);bag_curve.Add_Control_Point((T)10.,(T)1.3425);
            break;
        case 12:
            frame_rate=96;
            fluids_parameters.reincorporate_removed_particle_velocity=true;
            fluids_parameters.removed_particle_mass_scaling=60;
            last_frame=2000;
            fluids_parameters.reseeding_frame_rate=3;
            solids_parameters.implicit_solve_parameters.cg_iterations=200;
            fluids_parameters.grid->Initialize(15*resolution+1,15*resolution+1,10*resolution+1,(T)0,(T)1.5,(T)0,(T)1.5,(T)0,(T)1);
            break;
        case 13:
            solids_parameters.rigid_body_collision_parameters.use_push_out=true;
            if(solid_node) solids_parameters.use_rigid_deformable_contact=true;
            fountain_source.Append(CYLINDER<T>(TV((T)-1,(T)1.7,(T)0),TV((T)-.8,(T)1.7,(T)0),(T).2));
            fountain_source_velocity.Append(TV((T)3,(T)0,(T)0));
            fluids_parameters.grid->Initialize(20*resolution+1,20*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)2,(T)-1,(T)1);
            break;
        case 14:
            solids_parameters.use_rigid_deformable_contact=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            fluids_parameters.density=(T)1000;
            solids_parameters.triangle_collision_parameters.check_mesh_for_self_intersection=false;
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=0;
            solids_parameters.implicit_solve_parameters.cg_iterations=500;
            solid_density=(T)1200;
            frame_rate=48;
            fluids_parameters.grid->Initialize(20*resolution+1,25*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)2.5,(T)-1,(T)1);
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T).45*(T)2.5;
            break;
        case 17:
            solids_parameters.use_rigid_deformable_contact=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            fluids_parameters.density=(T)1000;
            solids_parameters.triangle_collision_parameters.allow_intersections=true;
            solids_parameters.triangle_collision_parameters.check_mesh_for_self_intersection=false;
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=0;
            solids_parameters.implicit_solve_parameters.cg_iterations=500;
            //solid_density=(T)1200;
            solid_density=(T)1900;
            frame_rate=48;
            fluids_parameters.grid->Initialize(20*resolution+1,25*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)2.5,(T)-1,(T)1);
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T).45*(T)2.5;
            break;
        case 21:
            last_frame=500;
            initial_fluid_height=(T).4;
            light_sphere_initial_height=(T).8;
            fluids_parameters.grid->Initialize(15*resolution+1,20*resolution+1,10*resolution+1,(T)0,(T)1.5,(T)0,(T)2,(T)0,(T)1);
            break;
        case 18:
            frame_rate=48;
            last_frame=250;
            fluids_parameters.density=(T)100;
//                solid_density=(T)6;                // This appears to be close to buoyancy; a little on the heavy side...
            solid_density=(T)1;
            fluids_parameters.grid->Initialize(40*resolution+1,25*resolution+1,20*resolution+1,(T)-2,(T)2,(T)0,(T)2.5,(T)-1,(T)1);
            solids_parameters.use_rigid_deformable_contact=false;
            solids_parameters.triangle_collision_parameters.allow_intersections=true;
            solids_parameters.triangle_collision_parameters.check_mesh_for_self_intersection=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=1;
            solids_parameters.implicit_solve_parameters.cg_iterations=500;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=500;

            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=false;
            initial_fluid_height=(T)1.5;

            fountain_source_boxes.Append(BOX<TV>(TV((T)-2,(T)0,(T)-1),TV((T)-1.95,initial_fluid_height,(T)1)));
            fountain_source_velocity.Append(TV((T).5,(T)0,(T)0));
            fountain_source_boxes.Append(BOX<TV>(TV((T)1.95,(T)0,(T)-1),TV((T)2,initial_fluid_height,(T)1)));
            fountain_source_velocity.Append(TV((T).5,(T)0,(T)0));
            break;
        case 19:
            frame_rate=48;
            last_frame=250;
            solid_density=(T)35;                // This appears to be close to buoyancy; a little on the heavy side...
            fluids_parameters.grid->Initialize(20*resolution+1,20*resolution+1,20*resolution+1,(T)-1,(T)1,(T)0,(T)2,(T)-1,(T)1);
            solids_parameters.use_rigid_deformable_contact=false;
            solids_parameters.triangle_collision_parameters.allow_intersections=true;
            solids_parameters.triangle_collision_parameters.check_mesh_for_self_intersection=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=1;
            solids_parameters.implicit_solve_parameters.cg_iterations=500;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1e-3;
            initial_fluid_height=(T)1;

            fountain_source.Append(CYLINDER<T>(TV((T).6,(T).7,(T).2),TV((T).8,(T).7,(T).2),(T).25));
            fountain_source_velocity.Append(TV((T)-1,(T)0,(T)0));
            fountain_source.Append(CYLINDER<T>(TV((T)-.6,(T).7,(T)-.2),TV((T)-.8,(T).7,(T)-.2),(T).25));
            fountain_source_velocity.Append(TV((T)1,(T)0,(T)0));
            break;
        case 20:
            frame_rate=48;
            last_frame=250;
            solid_density=(T)35;                // This appears to be close to buoyancy; a little on the heavy side...
            fluids_parameters.grid->Initialize(20*resolution+1,20*resolution+1,10*resolution+1,(T)-1,(T)1,(T)0,(T)2,(T)-.5,(T).5);
            solids_parameters.use_rigid_deformable_contact=false;
            solids_parameters.triangle_collision_parameters.allow_intersections=true;
            solids_parameters.triangle_collision_parameters.check_mesh_for_self_intersection=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=1;
            solids_parameters.implicit_solve_parameters.cg_iterations=500;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1e-4;
            initial_fluid_height=(T)1;

            curve.Add_Control_Point((T)0,TV((T).7,(T).5,(T)0));
            curve.Add_Control_Point((T)2,TV((T)0,(T)1.5,(T)0));
            curve.Add_Control_Point((T)10,TV((T)0,(T)1.5,(T)0));
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);
    }

    if(parse_args->Is_Value_Set("-iterations")){
        solids_parameters.implicit_solve_parameters.cg_iterations=parse_args->Get_Integer_Value("-iterations");output_directory+=STRING_UTILITIES::string_sprintf("_%f_iterations",solids_parameters.implicit_solve_parameters.cg_iterations);}

    switch(test_number){
        case 9:THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this,(T).5,(T).5,&rigid_bodies_to_collide_against);break;
        case 12:THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this,(T).5,(T).4);break;
        default:THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);break;}

    fluids_parameters.domain_walls[2][1]=true;

    // give mon hints
    LOG::cout<<"MONITOR begin_frame="<<this->first_frame<<std::endl;
    LOG::cout<<"MONITOR output_directory="<<(FILE_UTILITIES::Get_Working_Directory()+"/"+output_directory)<<std::endl;
    LOG::cout<<"MONITOR end_frame="<<last_frame<<std::endl;
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
//void Clear_Constrained_Particles(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
//{
//    for(int i=1;i<=constrained_node_positions.m;i++) particle_is_simulated(constrained_node_positions(i).x)=false;
//}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    if(test_number==9){
        ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
        T desired_x=(T)two_pi/16;
        ROTATION<TV> desired_rotation=ROTATION<TV>(desired_x*sin(4*time),TV(0,1,0));
        for(int i=1;i<=arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(desired_rotation);}}
}
//#####################################################################
// Function Initialize_Joint_Between
//#####################################################################
void Initialize_Joint_Between(JOINT<TV>* joint,const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,TV up)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.joint_mesh.Add_Articulation(parent.particle_index,child.particle_index,joint);
    const FRAME<TV>& pf=parent.Frame();
    const FRAME<TV>& cf=child.Frame();
    TV x=(cf.t-pf.t).Normalized();
    TV y=up.Projected_Orthogonal_To_Unit_Direction(x).Normalized();
    TV z=TV::Cross_Product(x,y);
    FRAME<TV> J((T).5*(cf.t+pf.t),ROTATION<TV>(MATRIX<T,3>(x,y,z)));
    joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(pf));
    joint->Set_Child_To_Joint_Frame(J.Inverse_Times(cf));
}
//#####################################################################
// Function Floppy_Fish
//#####################################################################
void Floppy_Fish()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.enforce_poststabilization_in_cg=false;
    solids_parameters.use_post_cg_constraints=false;
    arb.Set_Iterative_Tolerance((T)1e-4);
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    arb.Set_Poststabilization_Iterations(5);
    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;

    T friction=(T).2;
    T scale=(T)2;
    T bone_density=200;
    T bone_unscaled_volume=.1875;
    RIGID_BODY_STATE<TV> fish_state(FRAME<TV>(TV(0,(T)3,0),ROTATION<TV>((T)half_pi,TV(1,0,0))));

    ARRAY<RIGID_BODY<TV>*> bones;
    T bone_scales[5]={(T)1,(T)1,(T)1,(T).7,(T).6};
    TV bone_positions[5]={TV(-scale*(T)2,(T)3,0),TV(-scale*(T).75,(T)3,0),TV(scale*(T).5,(T)3,0),TV(scale*(T)1.7,(T)3,0),TV(scale*(T)2.5,(T)3,0)};
    for(int i=0;i<5;i++){
        T bone_scale=scale*bone_scales[i];
        RIGID_BODY<TV>& bone=solids_tests.Add_Rigid_Body("miniplank25wide2",bone_scale,friction);
        bones.Append(&bone);
        bone.Set_Frame(FRAME<TV>(bone_positions[i]));
        bone.Set_Mass(bone_density*bone_unscaled_volume*std::pow(bone_scale,TV::dimension));}

    T joint_strengths[4]={(T)500,(T)500,(T)200,(T)100};
    for(int i=2;i<=bones.m;i++){
        JOINT<TV>* joint=new POINT_JOINT<TV>;
        Initialize_Joint_Between(joint,*bones(i-1),*bones(i),TV(0,0,1));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(joint_strengths[i-2]);joint_function->Set_Target_Angle(ROTATION<TV>());}

    // add the fish
    T flesh_density=(T)200;
    TETRAHEDRALIZED_VOLUME<T>* fish=0;
    if(!sub_test){
        BOX<TV> box(-TV((T)3.93278,(T)1.07277,(T)0.384066),TV((T)2.68344,(T)1.1747,(T)0.384353));box*=scale;
        VECTOR<int,3> counts(20,15,5);
        GRID<TV> fish_mattress_grid=GRID<TV>(counts,box);
        solids_tests.Create_Mattress(fish_mattress_grid,true,&fish_state,flesh_density);
        fish_bounding_box=ORIENTED_BOX<TV>(box,fish_state.frame);}
    else{
        fish=&solids_tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",fish_state,true,false,flesh_density,scale);}

    // binding the deformable particles to the rigid bodies
    for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++) solids_tests.Bind_Particles_In_Rigid_Body(rigid_body_collection.Rigid_Body(p));

    if(fish){
        fish->Update_Number_Nodes();fish->Initialize_Triangulated_Surface();
        TRIANGULATED_SURFACE<T>& triangulated_surface=*fish->triangulated_surface;
        triangulated_surface.Update_Triangle_List();triangulated_surface.Initialize_Hierarchy();
        fish_levelset=solids_tests.Read_Or_Initialize_Implicit_Surface(STRING_UTILITIES::string_sprintf("%s/fish_undeformed_levelset.phi",output_directory.c_str()),triangulated_surface);}
}
//#####################################################################
// Function Water_Test_Number
//#####################################################################
static int Water_Test_Number(const int test_number)
{
    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
        case 17:
        case 18:
        case 19:
        case 20:
        case 21:
            return 1;
        default:
            return 1;}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
T Initial_Phi(const TV& X) const
{
    T phi=(T)1;
    switch(test_number){
        case 15:
        case 1: phi=X.y-(T).4*(T)1;break;
        case 3: phi=X.y-(T).6*(T)0.5;break;
        case 2:{
                /* 
                   case 1: heavy ball into bowl
                   case 2: light ball into bowl
                 */
            phi=X.y-(T).45;T solid_phi=solid_body_collection.rigid_body_collection.Rigid_Body(light_sphere_index).Implicit_Geometry_Extended_Value(X);
            if(solid_phi<0) phi=-solid_phi; // TODO Leave water outside bowl?
            break;}
        case 4:
            phi=X.y-(T).5;
            break;
        case 5:
            phi=X.y-(T).8;
            break;
        case 9:{
            T phi_fish=sub_test?(*fish_levelset)(X):fish_bounding_box.Signed_Distance(X);
            phi=max(X.y-initial_fluid_height,-phi_fish);
            break;}
        case 6:
        case 7:
        case 16:
        case 8:
        case 10:
        case 11:
        case 14:
        case 17:
        case 18:
        case 19:
        case 20:
        case 21:
            phi=X.y-initial_fluid_height;break;
        case 12:
            phi=max(X.y-(T).400235234,-solid_body_collection.rigid_body_collection.Rigid_Body(heavy_sphere_index).Implicit_Geometry_Extended_Value(X));
            break;
        case 13:
            phi=X.y-(T).1;
            break;
        default:
            phi=water_tests.Initial_Phi(X);}

    for(int i=1;i<=rigid_bodies_to_simulate.m;++i)
        phi=max(phi,-solid_body_collection.rigid_body_collection.Rigid_Body(rigid_bodies_to_simulate(i)).Implicit_Geometry_Extended_Value(X));
    for(int i=1;i<=fountain_source.m;++i) phi=min(phi,fountain_source(i).Signed_Distance(X));
    return phi;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=Initial_Phi(iterator.Location());
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=water_tests.Initial_Velocity(iterator.Location())[iterator.Axis()];

    if(test_number==18) for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(iterator.Axis()==1?(T)1:(T)0);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    for(int i=1;i<=fountain_source.m;i++) BASE::Get_Source_Velocities(fountain_source(i),world_to_source,fountain_source_velocity(i));
    for(int i=1;i<=fountain_source_boxes.m;i++) BASE::Get_Source_Velocities(fountain_source_boxes(i),world_to_source,fountain_source_velocity(i+fountain_source.m));
}
//#####################################################################
// Function Find_Placement
//#####################################################################
FRAME<TV> Find_Placement(RANDOM_NUMBERS<T>& random,const BOX<TV>& bounding_box,ARRAY<ORIENTED_BOX<TV> >& bounding_boxes,const BOX<TV>& world,bool want_rotate)
{
    for(int i=1;i<=10000;i++){
        FRAME<TV> frame;
        if(want_rotate) frame.r=random.template Get_Rotation<TV>();
        ORIENTED_BOX<TV> oriented_box(bounding_box,frame.r);
        BOX<TV> new_box(oriented_box.Bounding_Box());
        frame.t=random.Get_Uniform_Vector(world.min_corner-new_box.min_corner,world.max_corner-new_box.max_corner);
        oriented_box.corner+=frame.t;
        bool okay=true;
        for(int j=1;j<=bounding_boxes.m;j++) if(oriented_box.Intersection(bounding_boxes(j))){okay=false;break;}
        if(okay){
            bounding_boxes.Append(oriented_box);
            return frame;}}
    PHYSBAM_FATAL_ERROR("Could not find suitable placement");
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    // Add rigid bodies and initialize deformable objects into solids_tests
    switch(test_number){
        case 1:
        case 3:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
            tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
            //tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
            tetrahedralized_volume.Check_Signed_Volumes_And_Make_Consistent(true);
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(tetrahedralized_volume,solid_density,true);
            solids_tests.Copy_And_Add_Structure(tetrahedralized_volume);
            //int frames=16;
            //T initial_velocity=-(T)9.8*(T)frames/24;
            //for(int i=1;i<=particles.array_collection->Size();i++) particles.V(i)=TV((T)0,initial_velocity,(T)0);
            break;}
        case 15:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solids_tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_44K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T).5,(T).75,(T).5),ROTATION<TV>(-(T)half_pi,TV::Axis_Vector(1)))),true,true,1000,(T).15);
            
//                mattress_grid=T_GRID(sub_test*20,sub_test*5,sub_test*20,(T).3,(T).7,(T).6,(T).72,(T).3,(T).7);
            tetrahedralized_volume.Check_Signed_Volumes_And_Make_Consistent(true);
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(tetrahedralized_volume,solid_density,true);
            break;}
        case 2:{
            T diameter=right_corner.x-left_corner.x;
            // Set up the thing we're dropping our object into
            RIGID_BODY<TV>& bowl=solids_tests.Add_Rigid_Body("bowl",diameter*(T).55,(T)0);
            bowl.X()=TV((T)0,(T).05,(T)0);
            bowl.Rotation()=ROTATION<TV>(-(T)pi/2,TV::Axis_Vector(1));
            bowl.is_static=true;
            bowl.Set_Mass((T)1e10);
            bowl.Set_Coefficient_Of_Restitution((T)0);bowl.Set_Coefficient_Of_Friction((T)1);
            light_sphere_index=bowl.particle_index;
            rigid_bodies_to_simulate.Append(bowl.particle_index);
            
            // Set up the thing we're dropping
            RIGID_BODY<TV>& dropped_object=solids_tests.Add_Rigid_Body("sphere",(T).8*(T).5*diameter,(T)0);
            dropped_object.X()=TV((T)0,(T)2,(T)0);
            dropped_object.Set_Name("falling body");
            dropped_object.Set_Mass(sub_test==1?(T)1e10:(T)1e-10);
            dropped_object.Set_Coefficient_Of_Restitution((T)0);dropped_object.Set_Coefficient_Of_Friction((T)1);
            rigid_bodies_to_simulate.Append(dropped_object.particle_index);
            break;}
        case 4:{
            // place a bunch of rigid bodies - blocks and spheres
            int num_spheres=bodies;int num_blocks=bodies;
            T scale=(T).04;
            T sphere_mass=(T)four_thirds_pi*cube(scale)*solid_density,block_mass=(T)cube(2*scale)*solid_density;
            ARRAY<ORIENTED_BOX<TV> > bounding_boxes;
            GRID<TV> grid;
            if(fluids_parameters.mpi_grid) grid=fluids_parameters.mpi_grid->global_grid;
            else grid=*fluids_parameters.grid;
            
            // TODO Rewrite to not depend on the grid dimensions; may differ in MPI
            BOX<TV> world(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y+(grid.domain.max_corner.y-grid.domain.min_corner.y)*(T).75,(T)1.5*grid.domain.max_corner.y,grid.domain.min_corner.z,grid.domain.max_corner.z);

            RANDOM_NUMBERS<T> random;

            for(int i=1;i<=num_spheres;i++){
                RIGID_BODY<TV>& rigid_body_sphere=solids_tests.Add_Rigid_Body("sphere",scale,(T)0,true,false);
                rigid_body_sphere.Update_Bounding_Box();
                rigid_body_sphere.Set_Frame(Find_Placement(random,rigid_body_sphere.axis_aligned_bounding_box,bounding_boxes,world,true));
                rigid_body_sphere.Set_Mass(sphere_mass);///random.Get_Uniform_Number((T).1,(T)2));
                rigid_bodies_to_simulate.Append(rigid_body_sphere.particle_index);}

            for(int i=1;i<=num_blocks;i++){
                RIGID_BODY<TV>& rigid_body_block=solids_tests.Add_Rigid_Body("subdivided_box",scale,(T)0);
                rigid_body_block.Update_Bounding_Box();
                rigid_body_block.Set_Frame(Find_Placement(random,rigid_body_block.axis_aligned_bounding_box,bounding_boxes,world,true));
                rigid_body_block.Set_Mass(block_mass);///random.Get_Uniform_Number((T).1,(T)2));
                rigid_bodies_to_simulate.Append(rigid_body_block.particle_index);}
            break;}
        case 5:{
            //std::string rigid_directory=data_directory+"/Rigid_Bodies"+(TV::m==3?"":"_2D");
            //int id=solid_body_collection.rigid_body_collection.Add_Rigid_Body(stream_type,rigid_directory+"/box",(T).4,true,true,false,false);
            //RIGID_BODY<TV>& rigid_body_square=solid_body_collection.rigid_body_collection.Rigid_Body(id);
            RIGID_BODY<TV>& rigid_body_square=solids_tests.Add_Rigid_Body("box",(T).4,(T)0);
            rigid_body_square.Set_Frame(FRAME<TV>(TV((T).4,(T)1.4,(T).4)));
            rigid_body_square.Set_Coefficient_Of_Restitution((T)0);
            rigid_body_square.Set_Name("box");
            T volume=cube((T).4);
            rigid_body_square.Set_Mass(volume*solid_density);
            rigid_bodies_to_simulate.Append(rigid_body_square.particle_index);
            break;}
        case 6:
            Boat();break;
        case 7:
        case 16:{
            // light cloth
            TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Cloth_Panel(75,(T)1,(T)1.333,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T).25,(T)1.15,(T)0))));
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_surface,(T)10,true);
            // rigid body
            RIGID_BODY<TV>& torus=solids_tests.Add_Rigid_Body("torus",(T).25,(T)0);
            torus.Set_Frame(FRAME<TV>(TV((T)0,(T)1.5,(T)0)));
            torus.Set_Name("torus");
            T volume=cube((T).15);
            torus.Set_Mass(volume*solid_density);
            rigid_bodies_to_simulate.Append(torus.particle_index);
            break;}
        case 17:{
            // light cloth
            TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Cloth_Panel(75,(T)1,(T)1.333,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T).25,(T)1.15,(T)0))));
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_surface,(T)10,true);
            // rigid body
            RIGID_BODY<TV>& torus=solids_tests.Add_Rigid_Body("torus",(T).25,(T).2);
            torus.Set_Frame(FRAME<TV>(TV((T)0,(T)2.15,(T)0)));
            torus.Set_Name("torus");
            T volume=torus.Volume();
            torus.Set_Mass(volume*solid_density);
            rigid_bodies_to_simulate.Append(torus.particle_index);
            break;}
        case 8:{
            Balloon();

//             RIGID_BODY<TV>& pipe = *new RIGID_BODY<TV>(rigid_body_collection);

//             CYLINDER<T> pipe_structure(TV((T)0,heavy_sphere_initial_height,(T)0),TV((T)0,heavy_sphere_initial_height+(T)10,(T)0),(T).2);
//             pipe.Add_Structure(*pipe_structure.Generate_Triangles(100,100));
//             DEBUG_UTILITIES::Debug_Breakpoint();
            RIGID_BODY<TV>& pipe=solids_tests.Add_Rigid_Body("funnel_revolve",(T).3,(T)0);
            pipe.Set_Name("funnel");
            pipe.is_static=true;pipe.thin_shell=true;
            pipe.Set_Frame(FRAME<TV>(TV((T)0,(T)heavy_sphere_initial_height,(T)0)));
            rigid_bodies_to_simulate.Append(pipe.particle_index);
            break;}
        case 9: Floppy_Fish();break;
        case 10:{
            T fluid_mass=(T)four_thirds_pi*cube((T).4)*fluids_parameters.density;
            // really light cork
            RIGID_BODY<TV>& light_cork=solids_tests.Add_Rigid_Body("sphere",(T).4,(T)0);
            light_cork.Set_Frame(FRAME<TV>(TV((T)-1.5,light_sphere_initial_height,(T)0)));
            light_cork.Set_Name("light_cork");light_cork.Set_Mass((T).1*fluid_mass);
            light_cork.Set_Coefficient_Of_Restitution((T)1);
            rigid_bodies_to_simulate.Append(light_cork.particle_index);

            // lightish cork
            RIGID_BODY<TV>& lightish_cork=solids_tests.Add_Rigid_Body("sphere",(T).4,(T)0);
            lightish_cork.Set_Frame(FRAME<TV>(TV((T)-.5,light_sphere_initial_height,(T)0)));
            lightish_cork.Set_Name("lightish_cork");lightish_cork.Set_Mass((T).5*fluid_mass);
            lightish_cork.Set_Coefficient_Of_Restitution((T)1);
            rigid_bodies_to_simulate.Append(lightish_cork.particle_index);

            // barely light cork
            RIGID_BODY<TV>& barely_cork=solids_tests.Add_Rigid_Body("sphere",(T).4,(T)0);
            barely_cork.Set_Frame(FRAME<TV>(TV((T).5,light_sphere_initial_height,(T)0)));
            barely_cork.Set_Name("heavy_cork");barely_cork.Set_Mass((T).9*fluid_mass);
            barely_cork.Set_Coefficient_Of_Restitution((T)1);
            rigid_bodies_to_simulate.Append(barely_cork.particle_index);

            // heavy cork
            RIGID_BODY<TV>& heavy_cork=solids_tests.Add_Rigid_Body("sphere",(T).4,(T)0);
            heavy_cork.Set_Frame(FRAME<TV>(TV((T)1.5,light_sphere_initial_height,(T)0)));
            heavy_cork.Set_Name("cork");heavy_cork.Set_Mass((T)10*fluid_mass);
            heavy_cork.Set_Coefficient_Of_Restitution((T)1);
            rigid_bodies_to_simulate.Append(heavy_cork.particle_index);
            break;}
        case 21:{
            T fluid_mass=(T)four_thirds_pi*cube((T).1)*fluids_parameters.density;
            for(int i=0;i<6;++i) for(int j=0;j<4;++j) {
                RIGID_BODY<TV>& sphere=solids_tests.Add_Rigid_Body("sphere",(T).1,(T)0);
                TV position=TV((T).125+(T)i*(T).25,light_sphere_initial_height,(T).125+(T)j*(T).25);
                sphere.Set_Frame(FRAME<TV>(position));
                sphere.Set_Name(STRING_UTILITIES::string_sprintf("cork_%i_%i",i,j));
                sphere.Set_Coefficient_Of_Restitution((T)1);
                switch(j) {
                    case 0: sphere.Set_Mass(fluid_mass*(T).1); break;
                    case 1: sphere.Set_Mass(fluid_mass*(T).5); break;
                    case 2: sphere.Set_Mass(fluid_mass*(T).9); break;
                    case 3: sphere.Set_Mass(fluid_mass*(T)10); break;}
                rigid_bodies_to_simulate.Append(sphere.particle_index);
                constrained_nodes.Append(TRIPLE<int,T,TV>(sphere.particle_index,(T).5*(T)j+(T)(i+1)*(T).225877/(T)6,position));}
                break;}
        case 11:
            Balloon();
            break;
        case 12:{
            T fraction=(T).2;
            T volume=(T)four_thirds_pi*cube((T).05);
            RIGID_BODY<TV>& sphere=solids_tests.Add_Rigid_Body("sphere",(T).1,(T)0);
            sphere.Set_Frame(FRAME<TV>(TV((T)1.25,(T).55,(T).5)));
            sphere.Twist().linear=TV((T)-6,(T)-7.23333333333,(T)0);
            sphere.Set_Mass(sub_test==1?(T)1e10:fraction*volume*fluids_parameters.density);
            sphere.Set_Coefficient_Of_Restitution((T)1);
            sphere.Set_Coefficient_Of_Friction((T).4);
            rigid_bodies_to_simulate.Append(sphere.particle_index);
            heavy_sphere_index=sphere.particle_index;
            break;}
        case 13:{
            // place a bunch of rigid bodies - blocks and spheres
            int num_spheres=bodies;int num_blocks=bodies;
            T scale=(T).04;
            T sphere_mass=(T)four_thirds_pi*cube(scale)*solid_density,block_mass=(T)cube(2*scale)*solid_density;
            ARRAY<ORIENTED_BOX<TV> > bounding_boxes;
            GRID<TV> grid;
            if(fluids_parameters.mpi_grid) grid=fluids_parameters.mpi_grid->global_grid;
            else grid=*fluids_parameters.grid;
            
            // TODO Rewrite to not depend on the grid dimensions; may differ in MPI
            BOX<TV> world(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y,grid.domain.min_corner.y+(T).04,grid.domain.min_corner.z,grid.domain.max_corner.z);

            // slide all the rigid body walls down
            for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) rigid_body_collection.rigid_body_particle.X(i).y-=(T).75;

            RANDOM_NUMBERS<T> random;

            for(int i=1;i<=num_spheres;i++){
                RIGID_BODY<TV>& rigid_body_sphere=solids_tests.Add_Rigid_Body("sphere",scale,(T).15,true,false);
                rigid_body_sphere.Update_Bounding_Box();
                rigid_body_sphere.Set_Frame(Find_Placement(random,rigid_body_sphere.axis_aligned_bounding_box,bounding_boxes,world,true));
                rigid_body_sphere.X().y=scale*(T)5;
                rigid_body_sphere.Set_Mass(sphere_mass);///random.Get_Uniform_Number((T).1,(T)2));
                rigid_bodies_to_simulate.Append(rigid_body_sphere.particle_index);}

            for(int i=1;i<=num_blocks;i++){
                RIGID_BODY<TV>& rigid_body_block=solids_tests.Add_Rigid_Body("subdivided_box",scale,(T).15);
                rigid_body_block.Update_Bounding_Box();
                rigid_body_block.Set_Frame(Find_Placement(random,rigid_body_block.axis_aligned_bounding_box,bounding_boxes,world,true));
                rigid_body_block.X().y=scale*(T)5;
                rigid_body_block.Set_Mass(block_mass);///random.Get_Uniform_Number((T).1,(T)2));
                rigid_bodies_to_simulate.Append(rigid_body_block.particle_index);}
            break;}
        case 14:{
            // light cloth
            TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Cloth_Panel(75,(T)1,(T)1.333,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)1.15,(T)0))));
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_surface,(T)10,true);
            // Set up the thing we're dropping
            RIGID_BODY<TV>& dropped_object=solids_tests.Add_Rigid_Body("sphere_66k",(T).125,(T)0,true,true);
            dropped_object.X()=TV((T)0,(T)1.5,(T)0);
            dropped_object.Set_Name("falling body");
            T volume=(T)four_thirds_pi*cube((T).25);
            dropped_object.Set_Mass(volume*solid_density);
            dropped_object.Set_Coefficient_Of_Restitution((T)0);
            rigid_bodies_to_simulate.Append(dropped_object.particle_index);
            break;}
        case 18:{
            T top_of_flag=(T)1.9;
            T pole_x_location=(T)-1.55;
            T flag_width=(T)1.4;

            // light cloth
            TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Cloth_Panel(70,flag_width,(T)2,
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(pole_x_location+flag_width,top_of_flag-flag_width/(T)2,(T)0),ROTATION<TV>((T)pi/2,TV::Axis_Vector(1)))));
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_surface,solid_density,true);

            BOX<TV> binding_box=BOX<TV>(TV((T)-2,(T)0,(T)-1),TV(pole_x_location+(T).01,(T)2,(T)1));
            for(int i=1;i<=triangulated_surface.particles.array_collection->Size();++i){TV& position=triangulated_surface.particles.X(i);
                if(binding_box.Lazy_Inside(position)){
                    constrained_node_positions.Append(PAIR<int,TV>(i,position));}}
            break;}
        case 19:{
            T flag_width=(T).75;
            // light cloth
            TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Cloth_Panel(70,flag_width,(T)2,
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)2,(T)0),ROTATION<TV>((T)pi/2,TV::Axis_Vector(3)))));
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_surface,solid_density,true);
            break;}
        case 20:{
            int resolution=60;
            T flag_width=(T).7;
            // light cloth
            TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Cloth_Panel(resolution,flag_width,(T)2,
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).5,(T)0))));
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_surface,solid_density,true);

            for(int i=2*resolution+1;i<=triangulated_surface.particles.array_collection->Size();i+=2*resolution+1)
                constrained_node_positions.Append(PAIR<int,TV>(i,triangulated_surface.particles.X(i)));

//             left_point=triangulated_surface.particles.array_collection->Size();
//             right_point=2*resolution+1;
            break;}
        default: LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.collision_body_list.Add_Bodies(*rigid_body_collection.rigid_geometry_collection.collision_body_list);

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    // Add deformable object collision structures and forces
    switch(test_number){
        case 1:
        case 15:
        case 3:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solids_tests.Initialize_Tetrahedron_Collisions(1,tetrahedralized_volume,solids_parameters.triangle_collision_parameters);
            //solid_body_collection.Add_Force(Create_Edge_Springs(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,tetrahedralized_volume,stiffness_ratio*20,(T)3));
            //solid_body_collection.Add_Force(Create_Tet_Springs(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,tetrahedralized_volume,stiffness_ratio*10,(T)3));
            bool limit_time_step_by_strain_rate=false;
            solid_body_collection.Add_Force(Create_Edge_Springs(tetrahedralized_volume,stiffness_ratio*20,(T)3,
                limit_time_step_by_strain_rate,(T).5,true,(T)0,true,implicit_springs));
            solid_body_collection.Add_Force(Create_Tet_Springs(tetrahedralized_volume,stiffness_ratio*10,(T)3,
                limit_time_step_by_strain_rate,(T).1,true,(T).5,true,(T)0,true,implicit_springs));
            deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(tetrahedralized_volume.Get_Boundary_Object()));
            break;}
        case 7:
        case 14:
        case 16:
        case 17:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,10/(1+sqrt((T)2)),(T)15,false,(T).1,true,(T)0,true,true));
            solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,10/(1+sqrt((T)2)),(T)15,false,(T).1,true,(T)0,true,true));
            deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_surface));
            break;}
        case 8:
        case 11:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,(T)6e3,(T)5,false,(T).1,true,(T)0,true,true));
            solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,10/(1+sqrt((T)2)),(T)15,false,(T).1,true,(T)0,true,true));
            deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_surface));
            break;}
        case 9:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            T stiffness=(T)2e5,damping=(T).03;
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1));
            deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(tetrahedralized_volume.Get_Boundary_Object()));
            break;}
        case 18:
        case 19:
        case 20:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,20/(1+sqrt((T)2)),(T)15,false));
            solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,20/(1+sqrt((T)2)),(T)8,false));

            deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_surface));
            break;}
        default:break;
    }

    // Add everything to the simulation
    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
    for(int i=1;i<=deformable_objects_to_simulate.m;++i){
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& collision_structure=*deformable_objects_to_simulate(i);
        collision_structure.object.Initialize_Hierarchy();
        Add_To_Fluid_Simulation(collision_structure);}
    for(int i=1;i<=rigid_bodies_to_simulate.m;++i){
        RIGID_BODY<TV>& rigid_body_to_add=rigid_body_collection.Rigid_Body(rigid_bodies_to_simulate(i));
        if(rigid_body_to_add.thin_shell) Add_Thin_Shell_To_Fluid_Simulation(rigid_body_to_add); else Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_to_add);}

    // add a ground and add it to the collision body list
    RIGID_BODY<TV>* ground=0;
    switch(test_number){
        case 10:
        case 21: ground=&solids_tests.Add_Ground((T).4,(T)0,(T)1);break;
        case 12: ground=&solids_tests.Add_Ground((T).4,(T)0,(T).15);break;
        default: ground=&solids_tests.Add_Ground();break;}
    if(ground){
        if(test_number==9) rigid_bodies_to_collide_against.Append(ground->particle_index);}

    // collide structures with the ground and walls only
    if(test_number==9){
        deformable_body_collection.collisions.Use_Structure_Collide_Collision_Body(true);
        for(int s=1;s<=deformable_body_collection.deformable_geometry.structures.m;s++) for(int r=1;r<=rigid_bodies_to_collide_against.m;r++)
            deformable_body_collection.collisions.structure_collide_collision_body(s).Set(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(rigid_bodies_to_collide_against(r)));
        for(int i=1;i<=solid_body_collection.solids_forces.m;i++)
            solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;}
}
//#####################################################################
// Function Boat
//#####################################################################
void Boat()
{
    // boat
    RIGID_BODY<TV>& boat=solids_tests.Add_Rigid_Body("Thin_Shells/canoe_high",(T).25,(T)1,false);
    boat_index=boat.particle_index;
    boat.Update_Bounding_Box();T dy=boat.X().y-boat.Axis_Aligned_Bounding_Box().min_corner.y;
    boat.X()=TV((T)0,initial_fluid_height+(T)1.05*dy,(T)0);
    boat.thin_shell=true;
    boat.Set_Mass(boat_mass);

    // light sphere
    RIGID_BODY<TV>& light_sphere=solids_tests.Add_Rigid_Body("sphere",(T).06,(T).4);
    light_sphere_index=light_sphere.particle_index;
    light_sphere.Set_Frame(FRAME<TV>(TV((T)0,light_sphere_initial_height,-(T).1)));
    T light_sphere_density=(T)50;
    light_sphere.Set_Mass(light_sphere.Volume()*light_sphere_density);

    // heavy sphere
    RIGID_BODY<TV>& heavy_sphere=solids_tests.Add_Rigid_Body("sphere",(T).06,(T).4);
    heavy_sphere_index=heavy_sphere.particle_index;
    heavy_sphere.Set_Frame(FRAME<TV>(TV((T)0,heavy_sphere_initial_height,(T).15)));
    heavy_sphere.Set_Mass(heavy_sphere.Volume()*(T)10000);

    rigid_bodies_to_simulate.Append(boat.particle_index);
    rigid_bodies_to_simulate.Append(light_sphere_index);
    rigid_bodies_to_simulate.Append(heavy_sphere_index);
}
//#####################################################################
// Function Balloon
//#####################################################################
void Balloon()
{
//    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Triangulated_Object(data_directory+"/Tetrahedralized_Volumes/sphere_66k.tri",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,light_sphere_initial_height,(T)0))),true,true,balloon_initial_radius*(T)2);
    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Triangulated_Object(data_directory+"/Rigid_Bodies/sphere.tri",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,light_sphere_initial_height,(T)0))),true,true,balloon_initial_radius*(T)2);
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_surface,(T)10,true);

    T cut_sphere_radius=(T).3;
    T seperation_distance=(T)1;

    SPHERE<TV> analytic_cut_sphere(TV((T)0,light_sphere_initial_height+seperation_distance,(T)0),cut_sphere_radius*(T)2);
    T a=((sqr(balloon_initial_radius)/(sqr(cut_sphere_radius)))*seperation_distance)/((T)1+sqr(balloon_initial_radius)/sqr(cut_sphere_radius));
    T height=light_sphere_initial_height+seperation_distance-a+(T).08;
    T dist=a*(balloon_initial_radius/cut_sphere_radius)-(T).2;

    ARRAY<int> deletion_list; // List of deleted triangles
    ARRAY<bool> is_constrained;
    is_constrained.Resize(triangulated_surface.particles.array_collection->Size());
    ARRAYS_COMPUTATIONS::Fill(is_constrained,false);

    for(int i=1;i<=triangulated_surface.mesh.elements.m;++i){
        const TRIANGLE_3D<T>& triangle=triangulated_surface.Get_Element(i);
        int node1,node2,node3;triangulated_surface.mesh.elements(i).Get(node1,node2,node3);
        if(analytic_cut_sphere.Lazy_Inside(triangle.x1) || analytic_cut_sphere.Lazy_Inside(triangle.x2) || analytic_cut_sphere.Lazy_Inside(triangle.x3)){
            if(analytic_cut_sphere.Lazy_Inside(triangle.x1) && analytic_cut_sphere.Lazy_Inside(triangle.x2) && analytic_cut_sphere.Lazy_Inside(triangle.x3))
                deletion_list.Append(i);
            else {
                if(analytic_cut_sphere.Lazy_Inside(triangle.x1)) is_constrained(node1)=true;
                if(analytic_cut_sphere.Lazy_Inside(triangle.x2)) is_constrained(node2)=true;
                if(analytic_cut_sphere.Lazy_Inside(triangle.x3)) is_constrained(node3)=true;}}}

    triangulated_surface.mesh.Delete_Elements(deletion_list);
    ARRAY<int> condensation_mapping;
    triangulated_surface.Discard_Valence_Zero_Particles_And_Renumber(condensation_mapping);

    for(int i=1;i<=is_constrained.m;++i)
        if(is_constrained(i)){
            TV& position=triangulated_surface.particles.X(condensation_mapping(i));
            T angle=atan2(position(3),position(1));
            position=TV(cos(angle)*dist,height,sin(angle)*dist);
            constrained_node_positions.Append(PAIR<int,TV>(condensation_mapping(i),triangulated_surface.particles.X(condensation_mapping(i))));}
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    for(int i=1;i<=fountain_source.m;i++) Adjust_Phi_With_Source(fountain_source(i),world_to_source);
    for(int i=1;i<=fountain_source_boxes.m;i++) Adjust_Phi_With_Source(fountain_source_boxes(i),world_to_source);
    return (fountain_source.m && fountain_source_boxes.m);
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==21)
        for(int i=1;i<=constrained_nodes.m;++i)
            if(time<constrained_nodes(i).y) X(constrained_nodes(i).x)=constrained_nodes(i).z;
    if(test_number==6){
        if(time<light_sphere_drop_time) X(light_sphere_index).y=light_sphere_initial_height;
        if(time<heavy_sphere_drop_time) X(heavy_sphere_index).y=heavy_sphere_initial_height;}
    if(test_number==4){
        // slide all the rigid body walls down
        for(int i=1;i<=4;i++) X(i).y=(T).75;
    }
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==8 && time<heavy_sphere_drop_time)
        for(int i=1;i<=constrained_node_positions.m;++i){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            X(node_pair.x)=node_pair.y;}
    if(test_number==11)
        for(int i=1;i<=constrained_node_positions.m;++i){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            X(node_pair.x)=node_pair.y;X(node_pair.x).y+=bag_curve.Value(time);}
    if(test_number==18)
        for(int i=1;i<=constrained_node_positions.m;++i){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            X(node_pair.x)=node_pair.y;}
    if(test_number==20){
        for(int i=1;i<=constrained_node_positions.m;++i){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            X(node_pair.x)=curve.Value(time);
            X(node_pair.x)(3)=node_pair.y(3);}}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==6){
        if(velocity_time<light_sphere_drop_time) twist(light_sphere_index).linear.y=(T)0;
        if(velocity_time<heavy_sphere_drop_time) twist(heavy_sphere_index).linear.y=(T)0;}
   if(test_number==21)
        for(int i=1;i<=constrained_nodes.m;++i) if(velocity_time<constrained_nodes(i).y) twist(constrained_nodes(i).x).linear=TV();
 }
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==8 && velocity_time<heavy_sphere_drop_time)
        for(int i=1;i<=constrained_node_positions.m;++i){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            V(node_pair.x)=TV();}
    if(test_number==11)
        for(int i=1;i<=constrained_node_positions.m;++i){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            V(node_pair.x)=TV((T)0,bag_curve.Derivative(velocity_time),(T)0);}
    if(test_number==18)
        for(int i=1;i<=constrained_node_positions.m;++i){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            V(node_pair.x)=TV();}
    if(test_number==20)
        for(int i=1;i<=constrained_node_positions.m;++i)
            V(constrained_node_positions(i).x)=curve.Derivative(velocity_time);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==18)
        for(int i=1;i<=constrained_node_positions.m;++i){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            V(node_pair.x)=TV();}
    if(test_number==20)
        for(int i=1;i<=constrained_node_positions.m;++i)
            V(constrained_node_positions(i).x)=TV();
}
//#####################################################################
};
}
#endif
