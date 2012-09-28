//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Falling deformable sphere and smoke
//   2. Floppy bar test
//   4. Balloon filling with smoke then shooting off
//   5. Balloon shooting off at angle, with ceiling (hole in middle)
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::data_directory;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::solids_evolution;
    using BASE::Add_To_Fluid_Simulation;using BASE::Add_Thin_Shell_To_Fluid_Simulation;

    SMOKE_STANDARD_TESTS_3D<GRID<TV> > smoke_tests;
    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    int deformable_circle_id;
    GRID<TV> mattress_grid;
    ARRAY<int> deformable_object_enslaved_nodes;
    int rigid_body_id;
    T velocity_multiplier;
    int number_side_panels;
    T aspect_ratio,side_length;
    T stiffness_multiplier,damping_multiplier,bending_stiffness_multiplier,bending_damping_multiplier;
    int heavy_sphere_index,light_sphere_index;
    T heavy_sphere_drop_time,light_sphere_drop_time,heavy_sphere_initial_height,light_sphere_initial_height;
    T balloon_initial_radius;
    CYLINDER<T> source_cylinder;
    TV source_velocity;

    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_each_matrix;
    bool output_iterators;

    ARRAY<PAIR<int,TV> > constrained_node_positions;
    ARRAY<int> rigid_bodies_to_simulate;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.SMOKE),
        smoke_tests(*this,fluids_parameters,fluid_collection.incompressible_fluid_collection,solid_body_collection.rigid_body_collection),
        solids_tests(*this,solid_body_collection),deformable_circle_id(0),rigid_body_id(0),velocity_multiplier((T)1),number_side_panels(40),aspect_ratio((T)1.1),side_length((T).5),
        stiffness_multiplier((T)1),damping_multiplier((T)1),bending_stiffness_multiplier((T)1),bending_damping_multiplier((T)1),heavy_sphere_drop_time((T)1.5),source_velocity((T)0,(T).5,(T)0),
        run_self_tests(false),print_poisson_matrix(false),print_index_map(false),print_matrix(false),print_each_matrix(false),output_iterators(false)
    {
    }

    // Unused callbacks
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add("-cg_iterations",&fluids_parameters.incompressible_iterations,"value","cg iterations");
    parse_args->Add("-slip",&fluids_parameters.use_slip,"use slip");
    parse_args->Add("-viscosity",&fluids_parameters.viscosity,"value","fluid viscosity");
    parse_args->Add("-test_system",&run_self_tests,"Run self tests");
    parse_args->Add("-print_poisson_matrix",&print_poisson_matrix,"print_poisson_matrix");
    parse_args->Add("-print_index_map",&print_index_map,"print_index_map");
    parse_args->Add("-print_matrix",&print_matrix,"print_matrix");
    parse_args->Add("-print_each_matrix",&print_each_matrix,"print_each_matrix");
    parse_args->Add("-output_iterators",&output_iterators,"output_iterators");
    parse_args->Add_Not("-no_preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"do not use preconditioner");
    parse_args->Add("-preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"preconditioner");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    smoke_tests.Initialize(Smoke_Test_Number(test_number),resolution);
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    last_frame=1000;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.use_post_cg_constraints=false; // TODO: if this isn't false rigids don't get affected
    solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;

    *fluids_parameters.grid=smoke_tests.grid;
    fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;
        
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1.1;

    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_preconditioner_for_slip_system=true;
    if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=true;
        
    //if(solid_node || !mpi) solids_parameters.use_rigid_deformable_contact=true;

    switch(test_number){
        case 1:
            last_frame=1000;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            fluids_parameters.density=(T)1000;
            fluids_parameters.gravity=(T)9.8;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
            source_velocity=TV((T)0,(T)1,(T)0);
            source_cylinder=CYLINDER<T>(TV((T).5,(T)0,(T).5),TV((T).5,(T).05,(T).5),(T).15);
            break;
        case 2:
            last_frame=1000;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.density=(T)1000;
            velocity_multiplier=(T)4;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
            break;
        case 3:
            last_frame=1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.density=(T)10;
            velocity_multiplier=(T)5;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
            break;
        case 4:
            last_frame=1000;
            //solids_parameters.implicit_solve_parameters.cg_iterations=400;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.density=(T)1;
            velocity_multiplier=(T)15;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
            balloon_initial_radius=(T).2;
            light_sphere_drop_time=(T).5;
            heavy_sphere_drop_time=(T)1.5;
            heavy_sphere_initial_height=(T).25;
            light_sphere_initial_height=(T).25;
            source_velocity=TV((T)0,(T)4,(T)0);
            source_cylinder=CYLINDER<T>(TV((T).5,(T)0,(T).5),TV((T).5,(T).05,(T).5),(T).05);
            solids_parameters.implicit_solve_parameters.cg_iterations=1500;
            break;
        case 5:
            last_frame=1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.density=(T)1;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
            balloon_initial_radius=(T).2;
            light_sphere_drop_time=(T).5;
            heavy_sphere_drop_time=(T)2.5;
            heavy_sphere_initial_height=(T).35;
            light_sphere_initial_height=(T).35;
            source_velocity=TV(2*(T)root_two,2*(T)root_two,(T)0);
            source_cylinder=CYLINDER<T>(TV((T).5,(T)0,(T).5),TV((T).5,(T).1,(T).5),(T).05);
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}
        
    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);

    switch(test_number){
        case 1:
            break;
        case 2:
            mattress_grid=GRID<TV>(8,2,2,(T)0,(T)4,(T)2,(T)3,(T)2,(T)3);
            //mattress_grid=GRID<TV>(2,2,2,(T).5,(T).6,(T).5,(T).6,(T).5,(T).6);
            solids_parameters.implicit_solve_parameters.cg_iterations=4000;
            break;
        case 3:
        case 4:
        case 5:
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}


    if(fluids_parameters.use_slip)
        output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d_Resolution_%d_slip",test_number,resolution);
    else
        output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d_Resolution_%d",test_number,resolution);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE
{}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==4 && time<heavy_sphere_drop_time)
        for(int i=0;i<constrained_node_positions.m;i++){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            X(node_pair.x)=node_pair.y;}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==4 && velocity_time<heavy_sphere_drop_time)
        for(int i=0;i<constrained_node_positions.m;i++){
            PAIR<int,TV>& node_pair=constrained_node_positions(i);
            V(node_pair.x)=TV();}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    //for(int i=0;i<deformable_object_enslaved_nodes.m;i++) V(deformable_object_enslaved_nodes(i))=TV();
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(fluids_parameters.use_slip){
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).run_self_tests=run_self_tests;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).output_iterators=output_iterators;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_matrix=print_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_each_matrix=print_each_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_poisson_matrix=print_poisson_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_index_map=print_index_map;}
}
//#####################################################################
// Function Smoke_Test_Number
//#####################################################################
static int Smoke_Test_Number(const int test_number)
{
    switch(test_number){
        case 1:
        case 2:
        case 3:
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
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=velocity_multiplier*smoke_tests.Initial_Velocity(iterator.Location())[iterator.Axis()];
    BASE::Initialize_Velocities();
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    if(smoke_tests.test_number==1) BASE::Adjust_Density_And_Temperature_With_Sources(source_cylinder,smoke_tests.world_to_source,smoke_tests.rho,fluids_parameters.temperature_products);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    if(smoke_tests.test_number==1 && time<heavy_sphere_drop_time) BASE::Get_Source_Velocities(source_cylinder,smoke_tests.world_to_source,source_velocity);
}
//#####################################################################
// Function Balloon
//#####################################################################
void Balloon()
{
    TV balloon_initial_position=TV((T).5,light_sphere_initial_height,(T).5);
    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_tests.Create_Triangulated_Object(data_directory+"/Rigid_Bodies/sphere.tri",RIGID_BODY_STATE<TV>(FRAME<TV>(balloon_initial_position)),true,true,balloon_initial_radius);
    //triangulated_surface.density=(T)10;
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_surface,(T)2,true);

    T cut_sphere_radius=(T).1;
    T seperation_distance=balloon_initial_radius+cut_sphere_radius-(T).05;

    SPHERE<TV> analytic_cut_sphere(balloon_initial_position-TV((T)0,seperation_distance,(T)0),cut_sphere_radius);
    T a=(sqr(seperation_distance)+sqr(balloon_initial_radius)-sqr(cut_sphere_radius))/(2*seperation_distance);
    T height=light_sphere_initial_height-balloon_initial_radius;//seperation_distance+a;//+(T).08;
    T dist=sqrt(sqr(balloon_initial_radius)-sqr(a));
    //RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("sphere",cut_sphere_radius,0);
    //rigid_body.Frame()=FRAME<TV>(TV((T)0,light_sphere_initial_height-seperation_distance,(T)0));

    ARRAY<int> deletion_list; // List of deleted triangles
    ARRAY<bool> is_constrained;
    is_constrained.Resize(triangulated_surface.particles.Size());
    is_constrained.Fill(false);

    for(int i=0;i<triangulated_surface.mesh.elements.m;i++){
        const TRIANGLE_3D<T>& triangle=triangulated_surface.Get_Element(i);
        int node1,node2,node3;triangulated_surface.mesh.elements(i).Get(node1,node2,node3);
        if(analytic_cut_sphere.Lazy_Inside(triangle.X.x) || analytic_cut_sphere.Lazy_Inside(triangle.X.y) || analytic_cut_sphere.Lazy_Inside(triangle.X.z)){
            if(analytic_cut_sphere.Lazy_Inside(triangle.X.x) && analytic_cut_sphere.Lazy_Inside(triangle.X.y) && analytic_cut_sphere.Lazy_Inside(triangle.X.z))
                deletion_list.Append(i);
            else {
                if(analytic_cut_sphere.Lazy_Inside(triangle.X.x)) is_constrained(node1)=true;
                if(analytic_cut_sphere.Lazy_Inside(triangle.X.y)) is_constrained(node2)=true;
                if(analytic_cut_sphere.Lazy_Inside(triangle.X.z)) is_constrained(node3)=true;}}}

    triangulated_surface.mesh.Delete_Elements(deletion_list);
    ARRAY<int> condensation_mapping;
    triangulated_surface.Discard_Valence_Zero_Particles_And_Renumber(condensation_mapping);

    // rotate all the points around (.5,0,.5) at a 45 degree angle
    FRAME<TV> frame;
    if(test_number==4) frame=FRAME<TV>(TV((T)0,(T)0,(T)0),ROTATION<TV>((T)0,TV((T)0,(T)0,(T)1)));
    else if(test_number==5) frame=FRAME<TV>(TV((T)0,(T)0,(T)0),ROTATION<TV>((T)pi/4,TV((T)0,(T)0,(T)1)));
    else PHYSBAM_NOT_IMPLEMENTED();
    for(int i=0;i<is_constrained.m;i++){
        if(is_constrained(i)){
            TV& position=triangulated_surface.particles.X(condensation_mapping(i));
            T angle=atan2(position(3)-balloon_initial_position.z,position(1)-balloon_initial_position.x);
            position=TV(balloon_initial_position.x+cos(angle)*dist,height,balloon_initial_position.z+sin(angle)*dist);
            position=frame*position;
            constrained_node_positions.Append(PAIR<int,TV>(condensation_mapping(i),triangulated_surface.particles.X(condensation_mapping(i))));}
        else if(condensation_mapping(i)){
            TV& position=triangulated_surface.particles.X(condensation_mapping(i));
            position=frame*position;}}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    solid_body_collection.Set_CFL_Number(10);

    if(test_number==1){
        RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("sphere",(T).1,(T)0);
        rigid_body.Frame()=FRAME<TV>(TV((T).5,(T).6,(T).5));
        rigid_body.Set_Coefficient_Of_Restitution((T)0);
        rigid_body.name="circle";
        rigid_body_id=rigid_body.particle_index;
        rigid_body.Set_Mass(4.3);
        TRIANGULATED_SURFACE<T>& surface=rigid_body.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
        surface.Initialize_Hierarchy();
        surface.mesh.Initialize_Adjacent_Elements();
    }
    else if(test_number==2){
        solids_tests.Create_Mattress(mattress_grid,true);
        deformable_object_enslaved_nodes.Resize(mattress_grid.counts.y*mattress_grid.counts.z);
        deformable_object_enslaved_nodes=IDENTITY_ARRAY<int>(mattress_grid.counts.y*mattress_grid.counts.z);
    }
    else if(test_number==3){
        RIGID_BODY_STATE<TV> state;
        state.frame.t=TV((T).5,(T).5,(T).5);
        solids_tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,&state);
    }
    else if(test_number==4){
        Balloon();}
    else if(test_number==5){
        // set up the ceiling
        RIGID_BODY<TV>& ceiling_one=solids_tests.Add_Rigid_Body("ground",(T).005,(T)0);
        ceiling_one.is_static=true;
        ceiling_one.Frame().t=TV((T)0,(T)1.5,(T).55);
        rigid_bodies_to_simulate.Append(ceiling_one.particle_index);
        RIGID_BODY<TV>& ceiling_two=solids_tests.Add_Rigid_Body("ground",(T).005,(T)0);
        ceiling_two.is_static=true;
        ceiling_two.Frame().t=TV((T)0,(T)1.5,(T)-.55);
        rigid_bodies_to_simulate.Append(ceiling_two.particle_index);
        RIGID_BODY<TV>& ceiling_three=solids_tests.Add_Rigid_Body("ground",(T).005,(T)0);
        ceiling_three.is_static=true;
        ceiling_three.Frame().t=TV((T).55,(T)1.5,(T)0);
        rigid_bodies_to_simulate.Append(ceiling_three.particle_index);
        RIGID_BODY<TV>& ceiling_four=solids_tests.Add_Rigid_Body("ground",(T).005,(T)0);
        ceiling_four.is_static=true;
        ceiling_four.Frame().t=TV((T)-.55,(T)1.5,(T)0);
        rigid_bodies_to_simulate.Append(ceiling_four.particle_index);
    }

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
    T solid_gravity=(T)9.8;
    if(test_number==1){
        solids_tests.Add_Ground();
        solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
        Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_id));}
    else if(test_number==2){
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(tetrahedralized_volume,2000);
        tetrahedralized_volume.Initialize_Hierarchy();
        solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)4e5,(T).45,(T).01)));
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(tetrahedralized_volume.Get_Boundary_Object());
        deformable_collisions.object.Initialize_Hierarchy();
        Add_To_Fluid_Simulation(deformable_collisions);
    }
    else if(test_number==3){
        solid_gravity=(T)9.8;
        solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
        solids_tests.Add_Ground();
        TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
        T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
        solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,linear_stiffness,linear_damping)); // were *2 and *10
        T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
        solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,bending_stiffness,bending_damping));
        triangulated_surface.Initialize_Hierarchy();
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_surface),0,true);}
    else if(test_number==4 || test_number==5){
        TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
        solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,(T)5,(T)3,false,(T).1,true,(T)0,true,true));
        solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,5/(1+sqrt((T)2)),(T)1,false,(T).1,true,(T)0,true,true));

        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_surface);
        deformable_collisions.object.Initialize_Hierarchy();
        Add_To_Fluid_Simulation(deformable_collisions);
        for(int i=0;i<rigid_bodies_to_simulate.m;i++){
            RIGID_BODY<TV>& rigid_body_to_add=rigid_body_collection.Rigid_Body(rigid_bodies_to_simulate(i));
            if(rigid_body_to_add.thin_shell) Add_Thin_Shell_To_Fluid_Simulation(rigid_body_to_add); else Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_to_add);}
        //Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(heavy_sphere_index));
        //solids_tests.Add_Ground();
    }

    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++) solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->compute_half_forces=true;
}
//#####################################################################
};
}
#endif
