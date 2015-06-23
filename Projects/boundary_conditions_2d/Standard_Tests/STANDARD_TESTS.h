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

#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_WENO.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
#include <Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>

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
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef BOUNDARY_CONDITIONS_CALLBACKS<TV> VBC;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::solids_fluids_parameters;using BASE::output_directory;using BASE::last_frame;
    using BASE::frame_rate;using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::solid_body_collection;using BASE::solids_evolution;
    using BASE::test_number;using BASE::resolution;using BASE::data_directory;using BASE::stream_type;

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
    bool use_maccormack;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args,0,fluids_parameters.SMOKE),
        solids_tests(stream_type_input,data_directory,solid_body_collection),velocity_multiplier(1),mass_multiplier(1),stiffness_multiplier((T)1),damping_multiplier((T)1),
        bending_stiffness_multiplier((T)1),bending_damping_multiplier((T)1),parameter(0),boundary_offset(0),extra_cells(0),base_resolution(4),run_self_tests(false),
        print_poisson_matrix(false),print_index_map(false),print_matrix(false),print_each_matrix(false),output_iterators(false),
        use_eno_advection(false),angle(0),outside_tolerance(0),viscosity_only(false),use_maccormack(false)
    {
        parse_args.Add("-mass",&mass_multiplier,"scale","scale mass");
        parse_args.Add("-cg_iterations",&fluids_parameters.incompressible_iterations,"iterations","solver iterations");
        parse_args.Add("-viscosity",&fluids_parameters.viscosity,"value","viscosity");
        parse_args.Add("-test_system",&run_self_tests,"run self tests");
        parse_args.Add("-print_poisson_matrix",&print_poisson_matrix,"print poisson matrix");
        parse_args.Add("-print_index_map",&print_index_map,"print index map");
        parse_args.Add("-print_matrix",&print_matrix,"print matrix");
        parse_args.Add("-print_each_matrix",&print_each_matrix,"print each matrix");
        parse_args.Add("-output_iterators",&output_iterators,"output iterators");
        parse_args.Add_Not("-no_preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"do not use preconditioner");
        parse_args.Add("-preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"use preconditioner");
        parse_args.Add("-leakproof",&solids_fluids_parameters.use_leakproof_solve,"use leakproof solve");
        parse_args.Add("-maccormack",&use_maccormack,"use maccormack");
        parse_args.Add("-boundary_offset",&boundary_offset,"value","boundary offset");
        parse_args.Add("-levelset_viscosity",&fluids_parameters.use_levelset_viscosity,"use levelset viscosity");
        parse_args.Add("-print_viscosity_matrix",&fluids_parameters.print_viscosity_matrix,"print viscosity matrix");
        parse_args.Add("-angle",&angle,"angle","angle");
        parse_args.Add("-base_resolution",&base_resolution,"res","base resolution");
        parse_args.Parse();

        solids_tests.data_directory=data_directory;
        last_frame=100;
        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;

        angle*=pi/180;
        solids_parameters.implicit_solve_parameters.cg_iterations=fluids_parameters.incompressible_iterations;
        if(test_number==2 || test_number==3 || test_number==4) fluids_parameters.use_psi_R=true;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.rigid_body_collision_parameters.use_push_out=false;
        solids_parameters.use_post_cg_constraints=false;
        solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;

        fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;

        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.use_trapezoidal_rule_for_velocities=false;
        solids_parameters.implicit_solve_parameters.cg_restart_iterations=1000;
        solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;

        fluids_parameters.use_slip=true;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_preconditioner_for_slip_system=true;

        solids_fluids_parameters.use_leakproof_solve=false;

        if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=true;

        fluids_parameters.use_removed_positive_particles=true;
        fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.write_removed_positive_particles=true;
        fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.removed_positive_particle_buoyancy_constant=0;
        //solid_body_collection.print_residuals=true;

        fluids_parameters.gravity=TV();
        fluids_parameters.density=(T)100;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_density=fluids_parameters.use_temperature=false;
        fluids_parameters.use_body_force=false;

        output_directory=LOG::sprintf("Standard_Tests/Test_%d_Resolution_%d",test_number,resolution);
        if(use_maccormack){
            fluids_parameters.use_maccormack_semi_lagrangian_advection=true;
            fluids_parameters.use_maccormack_for_incompressible=true;}

        fluids_parameters.use_second_order_pressure=true;
    }

    virtual ~STANDARD_TESTS()
    {
        if(test_number==3 || test_number==4 || test_number==5){
            LOG::cout<<std::setprecision(16)<<std::endl;
            
            LOG::cout<<"LOCATIONS ";
            for(int i=1;i<face_samples.m;i++) LOG::cout<<fluid_collection.grid.Face(face_samples(i))<<" ";
            LOG::cout<<std::endl;
            Print_Samples("CONVERGENCE",fluid_collection.incompressible_fluid_collection.face_velocities);
            LOG::cout<<"LOCATIONS ";
            for(int i=1;i<cell_samples.m;i++) LOG::cout<<fluid_collection.grid.X(cell_samples(i))<<" ";
            LOG::cout<<std::endl;
            //Print_Samples("PRESSURE",dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).pressure);
        }
    }

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override {}
    void Update_Time_Varying_Material_Properties(const T time) override {}
    void Initialize_Euler_State() override {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) override {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
    void Update_Solids_Parameters(const T time) override {}
    void Preprocess_Solids_Substep(const T time,const int substep) override {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) override {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
    void Limit_Solids_Dt(T& dt,const T time) override {}
    void Adjust_Density_And_Temperature_With_Sources(const T time) override {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) override {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) override {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) override {}
    void Post_Initialization() override {}
    void Preprocess_Substep(const T dt,const T time) override {}
    void Postprocess_Substep(const T dt,const T time) override {}
    void Limit_Dt(T& dt,const T time) override {}
    void Postprocess_Frame(const int frame) override {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
    if(ADVECTION_MACCORMACK_UNIFORM<TV,T,ADVECTION<GRID<TV>,T> >* mac=
        dynamic_cast<ADVECTION_MACCORMACK_UNIFORM<TV,T,ADVECTION<GRID<TV>,T> >*>(fluids_parameters.incompressible->advection)){
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
        T y=(fluids_parameters.grid->Face(f).y-domain.min_corner.y)/(domain.max_corner.y-domain.min_corner.y);
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
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<2> >& face_velocities,ARRAY<bool,FACE_INDEX<2> >& psi_N,const T time) override
{
    return;
    if(test_number==2){
        for(FACE_ITERATOR<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);it.Valid();it.Next()){
            face_velocities(it.Full_Index())=velocity_multiplier;psi_N(it.Full_Index())=true;}}
    if(test_number==3 || test_number==4){
        for(FACE_ITERATOR<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);it.Valid();it.Next()){
            face_velocities(it.Full_Index())=Get_Source_Velocity(it.Full_Index());psi_N(it.Full_Index())=true;}}
}
//#####################################################################
// Function Get_Reflection_Conditions
//#####################################################################
void Get_Reflection_Conditions(ARRAY<T,FACE_INDEX<2> >& psi_R,const T time)
{
    if(test_number==2 || test_number==3 || test_number==4){
        for(FACE_ITERATOR<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,2);it.Valid();it.Next())
            psi_R(it.Full_Index())=mass_multiplier;
        psi_R(FACE_INDEX<2>(2,TV_INT(fluids_parameters.grid->counts.x+1,1)))=mass_multiplier;}
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) override
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
void Initialize_Bodies() override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    switch(test_number){
        case 1: Rigid_Circle();break;
        case 2: case 3: Boundary_Conditions();break;
        case 4: Boundary_Conditions_Offset();break;
        case 5: Angled_Boundary_Conditions();break;
        default:break;}

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

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
void Mark_Outside(ARRAY<bool,FACE_INDEX<TV::m> >& outside) override
{
    if(test_number==2 || test_number==3 || test_number==4){
        for(FACE_ITERATOR<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION);it.Valid();it.Next())
            outside(it.Full_Index())=true;}
    if(test_number==5){
        for(FACE_ITERATOR<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){TV X=it.Location();
            for(int i=0;i<4;i++) if(planes[i].Signed_Distance(X)>=-outside_tolerance) outside(it.Full_Index())=true;}}
}
//#####################################################################
// Function Mark_Outside
//#####################################################################
void Mark_Outside(ARRAY<bool,TV_INT>& outside) override
{
    if(test_number==5){
        for(CELL_ITERATOR<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){TV X=it.Location();
            for(int i=0;i<4;i++) if(planes[i].Signed_Distance(X)>=-outside_tolerance) outside(it.index)=true;}}
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() override
{
    if(test_number==4)
        for(FACE_ITERATOR<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){TV X=it.Location()-(T).5;
            TV V=X.Orthogonal_Vector()+X;
            fluid_collection.incompressible_fluid_collection.face_velocities(it.Full_Index())=V(it.Axis());}
}
//#####################################################################
// Function Get_Viscosity_Boundary_Along_Ray
//#####################################################################
typename VBC::RAY_TYPE Get_Boundary_Along_Ray(const FACE_INDEX<TV::m>& f1,const FACE_INDEX<TV::m>& f2,T& theta,T& value) override
{
    VECTOR<RANGE<TV_INT>,TV::dimension> fi=fluids_parameters.grid->Face_Indices(0);
    TV X0=fluids_parameters.grid->Face(f1);
    TV X1=fluids_parameters.grid->Face(f2);

    if(test_number==2 || test_number==3 || test_number==4){
        if(f2.index.x==fi(f2.axis).min_corner.x-(f2.axis!=1)){ // left
            theta=(domain.min_corner.x-X0.x)/(X1.x-X0.x);
//            if(test_number==2) value=velocity_multiplier;
//            if((test_number==3 || test_number==4)) value=Get_Source_Velocity(f2);
            if(plane_types[0]!=VBC::slip) return plane_types[0];
            return f2.axis==1?VBC::free:VBC::noslip;}
//            return VBC::noslip;}

        if(f2.index.y==fi(f2.axis).min_corner.y-(f2.axis!=2)){ // bottom
            theta=(domain.min_corner.y-X0.y)/(X1.y-X0.y);
            if(plane_types[1]!=VBC::slip) return plane_types[1];
            return f2.axis==1?VBC::free:VBC::noslip;}

        if(f2.index.x==fi(f2.axis).max_corner.x+(f2.axis!=1)){ // right
            theta=(domain.max_corner.x-X0.x)/(X1.x-X0.x);
            if(plane_types[2]!=VBC::slip) return plane_types[2];
            return f2.axis==2?VBC::free:VBC::noslip;}

        if(f2.index.y==fi(f2.axis).max_corner.y+(f2.axis!=2)){ // top
            theta=(domain.max_corner.y-X0.y)/(X1.y-X0.y);
            if(plane_types[3]!=VBC::slip) return plane_types[3];
            return f2.axis==1?VBC::free:VBC::noslip;}}

    if(test_number==5){
        for(int i=0;i<4;i++) if(planes[i].Signed_Distance(X1)>=-outside_tolerance){
                if(i==0 && f1.axis==1){value=velocity_multiplier*sqr(1-fluids_parameters.grid->Face(f2).y);}
                planes[i].Segment_Line_Intersection(X0,X1,theta);
                theta=clamp(theta,(T)0,(T)1);
                return plane_types[i];}}

    LOG::cout<<f1<<"  "<<f2<<"  "<<fi<<std::endl;

    return VBC::unused;
}
//#####################################################################
// Function Get_Viscosity_Boundary_Along_Ray
//#####################################################################
typename VBC::RAY_TYPE Get_Boundary_Along_Ray(const TV_INT& c1,const TV_INT& c2,T& theta,T& value) override
{
    RANGE<TV_INT> di=fluids_parameters.grid->Domain_Indices();
    TV X0=fluids_parameters.grid->X(c1);
    TV X1=fluids_parameters.grid->X(c2);

    if(test_number==2 || test_number==3 || test_number==4){
        if(c2.x<di.min_corner.x){ // left
            theta=(domain.min_corner.x-X0.x)/(X1.x-X0.x);
//            if(test_number==2) value=velocity_multiplier;
//            if((test_number==3 || test_number==4)) value=Get_Source_Velocity(c2);
            if(plane_types[0]!=VBC::slip) return plane_types[0];
            return VBC::noslip;}

        if(c2.y<di.min_corner.y){ // bottom
            theta=(domain.min_corner.y-X0.y)/(X1.y-X0.y);
            if(plane_types[1]!=VBC::slip) return plane_types[1];
            return VBC::noslip;}

        if(c2.x>di.max_corner.x){ // right
            theta=(domain.max_corner.x-X0.x)/(X1.x-X0.x);
            if(plane_types[2]!=VBC::slip) return plane_types[2];
            return VBC::noslip;}

        if(c2.y>di.max_corner.y){ // top
            theta=(domain.max_corner.y-X0.y)/(X1.y-X0.y);
            if(plane_types[3]!=VBC::slip) return plane_types[3];
            return VBC::noslip;}}

    if(test_number==5){
        for(int i=0;i<4;i++) if(planes[i].Signed_Distance(X1)>=-outside_tolerance){
                if(i==0){value=velocity_multiplier*sqr(1-fluids_parameters.grid->X(c2).y);}
                planes[i].Segment_Line_Intersection(X0,X1,theta);
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
    rigid_body.Frame()=FRAME<TV>(TV((T).5,(T).5));
    rigid_body.Set_Coefficient_Of_Restitution((T)0);
    rigid_body.name="circle";
    T density=100;
    rigid_body.Set_Mass((T)pi*sqr(radius)*(T)density*mass_multiplier);
    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body);

    fluids_parameters.grid->Initialize(TV_INT(resolution,(T)1.5*resolution)+1,RANGE<TV>(TV(),TV(1,(T)1.5)));
    fluids_parameters.gravity.y=-(T)9.8;
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
    fluids_parameters.grid->Initialize(TV_INT()+resolution+1,RANGE<TV>::Unit_Box());
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

    fluids_parameters.grid->Initialize(TV_INT(resolution,resolution+extra_cells)+1,RANGE<TV>(TV(0,-(T)extra_cells/resolution),TV(1,1)));
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
    planes[i].x0=(p1+p2)/2;
    planes[i].normal=TV(d.y,-d.x);
    plane_types[i]=t;

    int id=solid_body_collection.rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies_2D/ground",len/200);
    solid_body_collection.rigid_body_collection.rigid_body_particles.frame(id).t=planes[i].x0;
    solid_body_collection.rigid_body_collection.rigid_body_particles.frame(id).r=ROTATION<VECTOR<T,2> >::From_Rotated_Vector(TV(0,1),planes[i].normal);
    solid_body_collection.rigid_body_collection.Rigid_Body(id).is_static=true;
}
//#####################################################################
// Function Boundary_Conditions
//#####################################################################
void Bounding_Quad_From_Corners(const TV& p1,const TV& p2,const TV& p3,const TV& p4,typename VBC::RAY_TYPE t01,typename VBC::RAY_TYPE t12,typename VBC::RAY_TYPE t23,typename VBC::RAY_TYPE t30)
{
    Bounding_Edge_From_Endpoints(0,p1,p2,t01);
    Bounding_Edge_From_Endpoints(1,p2,p3,t12);
    Bounding_Edge_From_Endpoints(2,p3,p4,t23);
    Bounding_Edge_From_Endpoints(3,p4,p1,t30);
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

    fluids_parameters.grid->Initialize(TV_INT(resolution,resolution+extra)+1,RANGE<TV>(TV(),TV(1,1+(T)extra/resolution)));
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
