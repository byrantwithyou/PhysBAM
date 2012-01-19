//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Rigid circle
//   2. Boundary condition test
//   3. Boundary condition test (second order spatial convergence)
//   4. Boundary condition test (second order spatial convergence, unaligned wall)
//   5. Boundary condition test (second order spatial convergence, angled wall)
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_WENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_FACE_INDEX.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{
const void* global_face_samples;
const void* global_cell_samples;
int global_resolution;
template<class T,int d> void Print_Samples(const char* label,const ARRAY<T,VECTOR<int,d> >& array)
{
    if(!global_resolution) return;
    LOG::cout<<"W@Z-"<<label<<" "<<global_resolution<<" "<<array.Subset(*(const ARRAY<VECTOR<int,d> >*)global_cell_samples)<<std::endl;
}
template<class T,int d> void Print_Samples(const char* label,const ARRAY<T,FACE_INDEX<d> >& array)
{
    if(!global_resolution) return;
    LOG::cout<<"W@Z-"<<label<<" "<<global_resolution<<" "<<array.Subset(*(const ARRAY<FACE_INDEX<d> >*)global_face_samples)<<std::endl;
}

template void Print_Samples<double,1>(char const*,ARRAY<double,FACE_INDEX<1> > const&);
template void Print_Samples<double,1>(char const*,ARRAY<double,VECTOR<int,1> > const&);
template void Print_Samples<double,2>(char const*,ARRAY<double,FACE_INDEX<2> > const&);
template void Print_Samples<double,2>(char const*,ARRAY<double,VECTOR<int,2> > const&);
template void Print_Samples<double,3>(char const*,ARRAY<double,FACE_INDEX<3> > const&);
template void Print_Samples<double,3>(char const*,ARRAY<double,VECTOR<int,3> > const&);
template void Print_Samples<float,1>(char const*,ARRAY<float,FACE_INDEX<1> > const&);
template void Print_Samples<float,1>(char const*,ARRAY<float,VECTOR<int,1> > const&);
template void Print_Samples<float,2>(char const*,ARRAY<float,FACE_INDEX<2> > const&);
template void Print_Samples<float,2>(char const*,ARRAY<float,VECTOR<int,2> > const&);
template void Print_Samples<float,3>(char const*,ARRAY<float,FACE_INDEX<3> > const&);
template void Print_Samples<float,3>(char const*,ARRAY<float,VECTOR<int,3> > const&);

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef BOUNDARY_CONDITIONS_CALLBACKS<TV> VBC;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::solids_fluids_parameters;using BASE::output_directory;using BASE::last_frame;
    using BASE::frame_rate;using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::solid_body_collection;using BASE::solids_evolution;
    using BASE::parse_args;using BASE::test_number;using BASE::resolution;using BASE::data_directory;using BASE::stream_type;

    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    T velocity_multiplier;
    T mass_multiplier;
    T stiffness_multiplier,damping_multiplier,bending_stiffness_multiplier,bending_damping_multiplier;
    int parameter;
    T boundary_offset;
    int extra_cells;
    int base_resolution;

    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_each_matrix;
    bool output_iterators;
    bool use_eno_advection;
    RANGE<TV> domain;
    LINE_2D<T> planes[4];
    typename VBC::RAY_TYPE plane_types[4];
    T angle;
    T outside_tolerance;
    bool viscosity_only;
    ARRAY<FACE_INDEX<TV::m> > face_samples;
    ARRAY<TV_INT> cell_samples;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.SMOKE),
        solids_tests(*this,solid_body_collection),velocity_multiplier(1),mass_multiplier(0),stiffness_multiplier((T)1),damping_multiplier((T)1),
        bending_stiffness_multiplier((T)1),bending_damping_multiplier((T)1),parameter(0),boundary_offset(0),extra_cells(0),base_resolution(4),run_self_tests(false),
        print_poisson_matrix(false),print_index_map(false),print_matrix(false),print_each_matrix(false),output_iterators(false),
        use_eno_advection(false),angle(0),outside_tolerance(0),viscosity_only(false)
    {
    }

    virtual ~STANDARD_TESTS()
    {
        if(test_number==3 || test_number==4 || test_number==5){
            LOG::cout<<std::setprecision(16)<<std::endl;
            
            LOG::cout<<"LOCATIONS ";
            for(int i=1;i<face_samples.m;i++) LOG::cout<<fluid_collection.grid.Axis_X_Face(face_samples(i))<<" ";
            LOG::cout<<std::endl;
            Print_Samples("CONVERGENCE",fluid_collection.incompressible_fluid_collection.face_velocities);
            LOG::cout<<"LOCATIONS ";
            for(int i=1;i<cell_samples.m;i++) LOG::cout<<fluid_collection.grid.X(cell_samples(i))<<" ";
            LOG::cout<<std::endl;
            //Print_Samples("PRESSURE",dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).pressure);
        }
    }

    // Unused callbacks
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Double_Argument("-mass",(T)1);
    parse_args->Add_Integer_Argument("-cg_iterations",3000);
    parse_args->Add_Double_Argument("-viscosity",(T)0);
    parse_args->Add_Option_Argument("-test_system");
    parse_args->Add_Option_Argument("-print_poisson_matrix");
    parse_args->Add_Option_Argument("-print_index_map");
    parse_args->Add_Option_Argument("-print_matrix");
    parse_args->Add_Option_Argument("-print_each_matrix");
    parse_args->Add_Option_Argument("-output_iterators");
    parse_args->Add_Option_Argument("-no_preconditioner");
    parse_args->Add_Option_Argument("-preconditioner");
    parse_args->Add_Option_Argument("-leakproof");
    parse_args->Add_Integer_Argument("-parameter",0);
    parse_args->Add_Option_Argument("-maccormack");
    parse_args->Add_Double_Argument("-boundary_offset",(T)0);
    parse_args->Add_Option_Argument("-levelset_viscosity");
    parse_args->Add_Option_Argument("-print_viscosity_matrix");
    parse_args->Add_Double_Argument("-angle",(T)0);
    parse_args->Add_Integer_Argument("-base_resolution",4);
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    last_frame=100;
    mass_multiplier=(T)parse_args->Get_Double_Value("-mass");
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;

    if(test_number==2 || test_number==3 || test_number==4) fluids_parameters.use_psi_R=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.use_post_cg_constraints=false;
    solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;

    fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;
    fluids_parameters.incompressible_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    solids_parameters.implicit_solve_parameters.cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    parameter=parse_args->Get_Integer_Value("-parameter");

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=1000;
    solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;

    fluids_parameters.use_slip=true;
    run_self_tests=parse_args->Is_Value_Set("-test_system");
    print_poisson_matrix=parse_args->Is_Value_Set("-print_poisson_matrix");
    print_index_map=parse_args->Is_Value_Set("-print_index_map");
    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    print_each_matrix=parse_args->Is_Value_Set("-print_each_matrix");
    output_iterators=parse_args->Is_Value_Set("-output_iterators");
    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-no_preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=false;

    solids_fluids_parameters.use_leakproof_solve=false;
    if(parse_args->Is_Value_Set("-leakproof")) solids_fluids_parameters.use_leakproof_solve=true;

    fluids_parameters.viscosity=parse_args->Get_Double_Value("-viscosity");
    if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=true;
    fluids_parameters.use_levelset_viscosity=parse_args->Is_Value_Set("-levelset_viscosity");
    fluids_parameters.print_viscosity_matrix=parse_args->Is_Value_Set("-print_viscosity_matrix");
    angle=parse_args->Get_Double_Value("-angle")*pi/180;

    fluids_parameters.use_removed_positive_particles=true;
    fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.write_removed_positive_particles=true;
    fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.store_particle_ids=true;
    fluids_parameters.removed_positive_particle_buoyancy_constant=0;
    //solid_body_collection.print_residuals=true;

    fluids_parameters.gravity=0;
    fluids_parameters.density=(T)100;
    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.use_body_force=false;

    boundary_offset=(T)parse_args->Get_Double_Value("-boundary_offset");
    base_resolution=parse_args->Get_Integer_Value("-base_resolution");
    
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d_Resolution_%d",test_number,resolution);
    if(parse_args->Is_Value_Set("-maccormack")){
        fluids_parameters.use_maccormack_semi_lagrangian_advection=true;
        fluids_parameters.use_maccormack_for_incompressible=true;}

    fluids_parameters.use_second_order_pressure=true;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
    if(ADVECTION_MACCORMACK_UNIFORM<GRID<TV>,T,ADVECTION<GRID<TV>,T> >* mac=
        dynamic_cast<ADVECTION_MACCORMACK_UNIFORM<GRID<TV>,T,ADVECTION<GRID<TV>,T> >*>(fluids_parameters.incompressible->advection)){
        mac->ensure_second_order=true;
        mac->clamp_extrema=false;}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time)
{
    BASE::Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
T Get_Source_Velocity(const FACE_INDEX<2>& f) const
{
    return 0;
    if(test_number==2) return velocity_multiplier;
    if(test_number==3 || test_number==4){
//        if(!viscosity_only && f.axis==2) return 0; // Only force vertically in viscosity-only case.
        bool top_d=plane_types[3]==VBC::noslip || (plane_types[3]==VBC::slip && f.axis==2);
        bool bot_d=plane_types[1]==VBC::noslip || (plane_types[1]==VBC::slip && f.axis==2);
        T y=(fluids_parameters.grid->Axis_X_Face(f).y-domain.min_corner.y)/(domain.max_corner.y-domain.min_corner.y);
        if(top_d && bot_d) return y*(1-y)*(2*y+5);
        if(!top_d && bot_d) return y*((16*y-39)*y+30)/6;
        if(top_d && !bot_d) return (1-y)*((32*y+5)*y+5)/6;
        return ((6-4*y)*y*y+5)/6;}
    return 0;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
T Get_Source_Velocity(const TV_INT& c) const
{
    return 0;
    if(test_number==2) return velocity_multiplier;
    if(test_number==3 || test_number==4){
//        if(!viscosity_only && f.axis==2) return 0; // Only force vertically in viscosity-only case.
        bool top_d=plane_types[3]==VBC::noslip;
        bool bot_d=plane_types[1]==VBC::noslip;
        T y=(fluids_parameters.grid->X(c).y-domain.min_corner.y)/(domain.max_corner.y-domain.min_corner.y);
        if(top_d && bot_d) return y*(1-y)*(2*y+5);
        if(!top_d && bot_d) return y*((16*y-39)*y+30)/6;
        if(top_d && !bot_d) return (1-y)*((32*y+5)*y+5)/6;
        return ((6-4*y)*y*y+5)/6;}
    return 0;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<2> >& face_velocities,ARRAY<bool,FACE_INDEX<2> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    return;
    if(test_number==2){
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,1);it.Valid();it.Next()){
            face_velocities(it.Full_Index())=velocity_multiplier;psi_N(it.Full_Index())=true;}}
    if(test_number==3 || test_number==4){
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,1);it.Valid();it.Next()){
            face_velocities(it.Full_Index())=Get_Source_Velocity(it.Full_Index());psi_N(it.Full_Index())=true;}}
}
//#####################################################################
// Function Get_Reflection_Conditions
//#####################################################################
void Get_Reflection_Conditions(ARRAY<T,FACE_INDEX<2> >& psi_R,const T time)
{
    if(test_number==2 || test_number==3 || test_number==4){
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,3);it.Valid();it.Next())
            psi_R(it.Full_Index())=mass_multiplier;
        psi_R(FACE_INDEX<2>(2,TV_INT(fluids_parameters.grid->counts.x+1,1)))=mass_multiplier;}
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(fluids_parameters.use_slip){dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).run_self_tests=run_self_tests;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).output_iterators=output_iterators;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_matrix=print_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_each_matrix=print_each_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_poisson_matrix=print_poisson_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_index_map=print_index_map;}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    switch(test_number){
        case 1: Rigid_Circle();break;
        case 2: case 3: Boundary_Conditions();break;
        case 4: Boundary_Conditions_Offset();break;
        case 5: Angled_Boundary_Conditions();break;
        default:break;}

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(deformable_body_collection.soft_bindings);
    deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);

    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
    for(int k=0;k<deformable_body_collection.deformables_forces.m;k++) deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
    for(int i=0;i<rigid_body_collection.rigids_forces.m;i++) rigid_body_collection.rigids_forces(i)->compute_half_forces=true;
}
//#####################################################################
// Function Mark_Outside
//#####################################################################
void Mark_Outside(ARRAY<bool,FACE_INDEX<TV::m> >& outside) PHYSBAM_OVERRIDE
{
    if(test_number==2 || test_number==3 || test_number==4){
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION);it.Valid();it.Next())
            outside(it.Full_Index())=true;}
    if(test_number==5){
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){TV X=it.Location();
            for(int i=0;i<4;i++) if(planes[i].Signed_Distance(X)>=-outside_tolerance) outside(it.Full_Index())=true;}}
}
//#####################################################################
// Function Mark_Outside
//#####################################################################
void Mark_Outside(ARRAY<bool,TV_INT>& outside) PHYSBAM_OVERRIDE
{
    if(test_number==5){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){TV X=it.Location();
            for(int i=0;i<4;i++) if(planes[i].Signed_Distance(X)>=-outside_tolerance) outside(it.index)=true;}}
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    if(test_number==4)
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){TV X=it.Location()-(T).5;
            TV V=X.Orthogonal_Vector()+X;
            fluid_collection.incompressible_fluid_collection.face_velocities(it.Full_Index())=V(it.Axis());}
}
//#####################################################################
// Function Get_Viscosity_Boundary_Along_Ray
//#####################################################################
typename VBC::RAY_TYPE Get_Boundary_Along_Ray(const FACE_INDEX<TV::m>& f1,const FACE_INDEX<TV::m>& f2,T& theta,T& value) PHYSBAM_OVERRIDE
{
    VECTOR<RANGE<TV_INT>,TV::dimension> fi=fluids_parameters.grid->Face_Indices(0);
    TV X1=fluids_parameters.grid->Axis_X_Face(f1);
    TV X2=fluids_parameters.grid->Axis_X_Face(f2);

    if(test_number==2 || test_number==3 || test_number==4){
        if(f2.index.x==fi(f2.axis).min_corner.x-(f2.axis!=1)){ // left
            theta=(domain.min_corner.x-X1.x)/(X2.x-X1.x);
//            if(test_number==2) value=velocity_multiplier;
//            if((test_number==3 || test_number==4)) value=Get_Source_Velocity(f2);
            if(plane_types[0]!=VBC::slip) return plane_types[0];
            return f2.axis==1?VBC::free:VBC::noslip;}
//            return VBC::noslip;}

        if(f2.index.y==fi(f2.axis).min_corner.y-(f2.axis!=2)){ // bottom
            theta=(domain.min_corner.y-X1.y)/(X2.y-X1.y);
            if(plane_types[1]!=VBC::slip) return plane_types[1];
            return f2.axis==1?VBC::free:VBC::noslip;}

        if(f2.index.x==fi(f2.axis).max_corner.x+(f2.axis!=1)){ // right
            theta=(domain.max_corner.x-X1.x)/(X2.x-X1.x);
            if(plane_types[2]!=VBC::slip) return plane_types[2];
            return f2.axis==2?VBC::free:VBC::noslip;}

        if(f2.index.y==fi(f2.axis).max_corner.y+(f2.axis!=2)){ // top
            theta=(domain.max_corner.y-X1.y)/(X2.y-X1.y);
            if(plane_types[3]!=VBC::slip) return plane_types[3];
            return f2.axis==1?VBC::free:VBC::noslip;}}

    if(test_number==5){
        for(int i=0;i<4;i++) if(planes[i].Signed_Distance(X2)>=-outside_tolerance){
                if(i==0 && f1.axis==1){value=velocity_multiplier*sqr(1-fluids_parameters.grid->Axis_X_Face(f2).y);}
                planes[i].Segment_Line_Intersection(X1,X2,theta);
                theta=clamp(theta,(T)0,(T)1);
                return plane_types[i];}}

    LOG::cout<<f1<<"  "<<f2<<"  "<<fi<<std::endl;

    return VBC::unused;
}
//#####################################################################
// Function Get_Viscosity_Boundary_Along_Ray
//#####################################################################
typename VBC::RAY_TYPE Get_Boundary_Along_Ray(const TV_INT& c1,const TV_INT& c2,T& theta,T& value) PHYSBAM_OVERRIDE
{
    RANGE<TV_INT> di=fluids_parameters.grid->Domain_Indices();
    TV X1=fluids_parameters.grid->X(c1);
    TV X2=fluids_parameters.grid->X(c2);

    if(test_number==2 || test_number==3 || test_number==4){
        if(c2.x<di.min_corner.x){ // left
            theta=(domain.min_corner.x-X1.x)/(X2.x-X1.x);
//            if(test_number==2) value=velocity_multiplier;
//            if((test_number==3 || test_number==4)) value=Get_Source_Velocity(c2);
            if(plane_types[0]!=VBC::slip) return plane_types[0];
            return VBC::noslip;}

        if(c2.y<di.min_corner.y){ // bottom
            theta=(domain.min_corner.y-X1.y)/(X2.y-X1.y);
            if(plane_types[1]!=VBC::slip) return plane_types[1];
            return VBC::noslip;}

        if(c2.x>di.max_corner.x){ // right
            theta=(domain.max_corner.x-X1.x)/(X2.x-X1.x);
            if(plane_types[2]!=VBC::slip) return plane_types[2];
            return VBC::noslip;}

        if(c2.y>di.max_corner.y){ // top
            theta=(domain.max_corner.y-X1.y)/(X2.y-X1.y);
            if(plane_types[3]!=VBC::slip) return plane_types[3];
            return VBC::noslip;}}

    if(test_number==5){
        for(int i=0;i<4;i++) if(planes[i].Signed_Distance(X2)>=-outside_tolerance){
                if(i==0){value=velocity_multiplier*sqr(1-fluids_parameters.grid->X(c2).y);}
                planes[i].Segment_Line_Intersection(X1,X2,theta);
                theta=clamp(theta,(T)0,(T)1);
                return plane_types[i];}}

    LOG::cout<<c1<<"  "<<c2<<"  "<<domain<<std::endl;

    return VBC::unused;
}
//#####################################################################
// Function Rigid_Circle
//#####################################################################
void Rigid_Circle()
{
    T radius=.3;
    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",(T)radius,(T)0);
    rigid_body.thin_shell=false;
    rigid_body.Set_Frame(FRAME<TV>(TV((T).5,(T).5)));
    rigid_body.Set_Coefficient_Of_Restitution((T)0);
    rigid_body.Set_Name("circle");
    T density=100;
    rigid_body.Set_Mass((T)pi*sqr(radius)*(T)density*mass_multiplier);
    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body);

    fluids_parameters.grid->Initialize(resolution+1,(int)(1.5*resolution)+1,(T)0,(T)1,(T)0,(T)1.5);
    fluids_parameters.gravity=(T)9.8;
    fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.domain_walls[1][2]=true;
    fluids_parameters.domain_walls[2][1]=true;
    fluids_parameters.domain_walls[2][2]=false;

    solids_tests.Add_Ground();
    solids_tests.Add_Gravity();
}
//#####################################################################
// Function Boundary_Conditions
//#####################################################################
void Boundary_Conditions()
{
    fluids_parameters.grid->Initialize(resolution+1,resolution+1,(T)0,(T)1,(T)0,(T)1);
    domain=fluids_parameters.grid->domain;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.domain_walls[2][1]=false;
    fluids_parameters.domain_walls[1][2]=false;
    fluids_parameters.domain_walls[2][2]=false;
    plane_types[0]=VBC::noslip;
    plane_types[1]=VBC::noslip;
    plane_types[2]=VBC::noslip;
    plane_types[3]=VBC::noslip;
    Compute_Sample_Points();
}
//#####################################################################
// Function Boundary_Conditions
//#####################################################################
void Boundary_Conditions_Offset()
{
    last_frame=5;
    T extra=boundary_offset*resolution+(T).5;
    extra_cells=(int)extra;
    T left=extra-extra_cells;
    mass_multiplier=(left-1)/left;

    fluids_parameters.grid->Initialize(resolution+1,resolution+1+extra_cells,(T)0,(T)1,-(T)extra_cells/resolution,(T)1);
    domain=fluids_parameters.grid->domain;
    domain.min_corner.y=-boundary_offset;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.domain_walls[2][1]=false;
    fluids_parameters.domain_walls[1][2]=false;
    fluids_parameters.domain_walls[2][2]=false;
    plane_types[0]=VBC::noslip;
    plane_types[1]=VBC::noslip;
    plane_types[2]=VBC::noslip;
    plane_types[3]=VBC::noslip;
    Compute_Sample_Points();
}
//#####################################################################
// Function Boundary_Conditions
//#####################################################################
void Bounding_Edge_From_Endpoints(int i,const TV& p1,const TV& p2,typename VBC::RAY_TYPE t)
{
    TV d=p2-p1;
    T len=d.Normalize();
    planes[i].x1=(p1+p2)/2;
    planes[i].normal=TV(d.y,-d.x);
    plane_types[i]=t;

    int id=solid_body_collection.rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies_2D/ground",len/200);
    solid_body_collection.rigid_body_collection.rigid_body_particle.X(id)=planes[i].x1;
    solid_body_collection.rigid_body_collection.rigid_body_particle.rotation(id)=ROTATION<VECTOR<T,2> >::From_Rotated_Vector(TV(0,1),planes[i].normal);
    solid_body_collection.rigid_body_collection.Rigid_Body(id).is_static=true;
}
//#####################################################################
// Function Boundary_Conditions
//#####################################################################
void Bounding_Quad_From_Corners(const TV& p1,const TV& p2,const TV& p3,const TV& p4,typename VBC::RAY_TYPE t12,typename VBC::RAY_TYPE t23,typename VBC::RAY_TYPE t34,typename VBC::RAY_TYPE t41)
{
    Bounding_Edge_From_Endpoints(0,p1,p2,t12);
    Bounding_Edge_From_Endpoints(1,p2,p3,t23);
    Bounding_Edge_From_Endpoints(2,p3,p4,t34);
    Bounding_Edge_From_Endpoints(3,p4,p1,t41);
}
//#####################################################################
// Function Compute_Sample_Points
//#####################################################################
void Compute_Sample_Points()
{
    PHYSBAM_ASSERT(resolution%base_resolution==0 && resolution/base_resolution%2==1);
    int m=resolution/base_resolution,n=m/2;
    for(int y=0;y<base_resolution;y++)
        for(int x=2; x<=base_resolution; x++)
            face_samples.Append(FACE_INDEX<TV::m>(1,TV_INT((x-1)*(m-1)+x,y*m-n+extra_cells)));
    for(int y=2; y<=base_resolution; y++)
        for(int x=0;x<base_resolution;x++)
            face_samples.Append(FACE_INDEX<TV::m>(2,TV_INT(x*m-n,(y-1)*(m-1)+y+extra_cells)));
    for(int y=0;y<base_resolution;y++)
        for(int x=0;x<base_resolution;x++)
            cell_samples.Append(TV_INT(x*m-n,y*m-n+extra_cells));

    global_face_samples=&face_samples;
    global_cell_samples=&cell_samples;
    global_resolution=resolution;
}
//#####################################################################
// Function Boundary_Conditions
//#####################################################################
void Angled_Boundary_Conditions()
{
    last_frame=5;
    TV x(1,tan(angle)),y(0,1);
    Bounding_Quad_From_Corners(y,TV(),x,x+y,VBC::noslip,VBC::noslip,VBC::free,VBC::free);

    for(int i=0;i<4;i++) PHYSBAM_ASSERT(planes[i].Signed_Distance((x+y)/2)<0);

    int extra=(int)ceil(((x+y).y-1)*resolution);

    fluids_parameters.grid->Initialize(resolution+1,resolution+1+extra,0,1,0,1+(T)extra/resolution);
    fluids_parameters.use_levelset_viscosity=true;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.domain_walls[1][2]=false;
    fluids_parameters.domain_walls[2][1]=false;
    fluids_parameters.domain_walls[2][2]=false;
    outside_tolerance=(T)1e-5;
    Compute_Sample_Points();
}
//#####################################################################
};
}
#endif
