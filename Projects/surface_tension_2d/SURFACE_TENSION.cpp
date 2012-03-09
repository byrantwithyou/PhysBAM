#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLECTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_CUT.h>
#include "SURFACE_TENSION.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4355) // 'this' : used in base member initializer list
#endif
template<class T> SURFACE_TENSION<T>::
SURFACE_TENSION(const STREAM_TYPE stream_type)
    :BASE(stream_type,1),solids_tests(*this,solid_body_collection),run_self_tests(false),print_poisson_matrix(false),print_index_map(false),
    print_matrix(false),print_each_matrix(false),use_full_ic(false),output_iterators(false),use_decoupled_viscosity(false),max_dt(0),exact_dt(0),current_dt(0),implicit_solid(true),use_cut_volume(false),use_low_order_advection(false),
    front_tracked_structure(0),rebuild_curve(0),deformable_collisions(0),fsi(0),number_surface_particles(0),rebuild_surface(false),free_particles(0),psi_D(0),
    circle_radius(0),circle_perturbation(0),oscillation_mode(0),use_massless_structure(false),coupled_particles(0),make_ellipse(false),use_phi(false),remesh(false),m(1),s(1),kg(1),solid_refinement(0),
    solid_density(0),solid_width(0),analytic_solution(0),omega(0),laplace_number(0),use_T_nu(false)
{
    LOG::cout<<std::setprecision(16);
    debug_particles.array_collection->template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    debug_particles.Store_Velocity(true);
    Store_Debug_Particles(&debug_particles);
}
#ifdef _MSC_VER
#pragma warning(pop)
#endif
//#####################################################################
// Destructor
//#####################################################################
template<class T> SURFACE_TENSION<T>::
~SURFACE_TENSION()
{
    delete deformable_collisions;
    if(!use_low_order_advection) delete fluids_parameters.incompressible->advection;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Register_Options()
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-cg_iterations",3000);
    parse_args->Add_Double_Argument("-viscosity",(T)0);
    parse_args->Add_Option_Argument("-test_system");
    parse_args->Add_Option_Argument("-print_poisson_matrix");
    parse_args->Add_Option_Argument("-print_index_map");
    parse_args->Add_Option_Argument("-print_matrix");
    parse_args->Add_Option_Argument("-print_each_matrix");
    parse_args->Add_Option_Argument("-use_full_ic");
    parse_args->Add_Option_Argument("-output_iterators");
    parse_args->Add_Option_Argument("-no_preconditioner");
    parse_args->Add_Option_Argument("-preconditioner");
    parse_args->Add_Option_Argument("-leakproof");
    parse_args->Add_Option_Argument("-use_decoupled_viscosity");
    parse_args->Add_Double_Argument("-max_dt",0);
    parse_args->Add_Double_Argument("-dt",0);
    parse_args->Add_Option_Argument("-explicit_solid");
    parse_args->Add_Option_Argument("-build_surface");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Double_Argument("-surface_tension",.0728);
    parse_args->Add_Option_Argument("-cut_cell");
    parse_args->Add_Integer_Argument("-oscillation_mode",2);
    parse_args->Add_Integer_Argument("-solid_refinement",7);
    parse_args->Add_Option_Argument("-make_ellipse");
    parse_args->Add_Double_Argument("-epsilon",.05);
    parse_args->Add_Double_Argument("-m",1,"length unit");
    parse_args->Add_Double_Argument("-s",1,"time unit");
    parse_args->Add_Double_Argument("-kg",1,"mass unit");
    parse_args->Add_Double_Argument("-linear_force",0);
    parse_args->Add_Double_Argument("-rand",0);
    parse_args->Add_Option_Argument("-phi");
    parse_args->Add_Option_Argument("-two_phase");
    parse_args->Add_Option_Argument("-cut_mass");
    parse_args->Add_Integer_Argument("-cut_order",4);
    parse_args->Add_Option_Argument("-remesh");
    parse_args->Add_Option_Argument("-low_order_advection");
    parse_args->Add_Double_Argument("-la",120);
    parse_args->Add_Option_Argument("-tnu");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Parse_Options()
{
    BASE::Parse_Options();
    last_frame=100;
    frame_rate=24;

    kg=(T)parse_args->Get_Double_Value("-kg");
    m=(T)parse_args->Get_Double_Value("-m");
    s=(T)parse_args->Get_Double_Value("-s");
    frame_rate/=s;
    laplace_number=(T)parse_args->Get_Double_Value("-la");
    use_T_nu=parse_args->Is_Value_Set("-tnu");

    fluids_parameters.cfl=(T).9;
    fluids_parameters.confinement_parameter=(T).04;
    fluids_parameters.rho_bottom=1;
    fluids_parameters.rho_top=(T).65;
    fluids_parameters.density_buoyancy_constant=fluids_parameters.temperature_buoyancy_constant=0;
    fluids_parameters.temperature_container.Set_Cooling_Constant(0);
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;

    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;
    fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;
    fluids_parameters.incompressible_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    solids_parameters.implicit_solve_parameters.cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=200;
    solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    //solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
    fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;

    fluids_parameters.use_slip=true;
//    fluids_parameters.use_slip=parse_args->Is_Value_Set("-slip");
    run_self_tests=parse_args->Is_Value_Set("-test_system");
    print_poisson_matrix=parse_args->Is_Value_Set("-print_poisson_matrix");
    print_index_map=parse_args->Is_Value_Set("-print_index_map");
    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    print_each_matrix=parse_args->Is_Value_Set("-print_each_matrix");
    use_full_ic=parse_args->Is_Value_Set("-use_full_ic");
    output_iterators=parse_args->Is_Value_Set("-output_iterators");
    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-no_preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=false;

    solids_fluids_parameters.use_leakproof_solve=false;
    if(parse_args->Is_Value_Set("-leakproof")) solids_fluids_parameters.use_leakproof_solve=true;

    fluids_parameters.viscosity=(T)(parse_args->Get_Double_Value("-viscosity")*kg/s);
    if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=true;
    fluids_parameters.surface_tension=(T)(parse_args->Get_Double_Value("-surface_tension")*kg*m/(s*s));
    solid_refinement=parse_args->Get_Integer_Value("-solid_refinement");

    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.store_particle_ids=true;
    // T default_removed_positive_particle_buoyancy_constant=fluids_parameters.removed_positive_particle_buoyancy_constant;
    fluids_parameters.removed_positive_particle_buoyancy_constant=0;
    //solid_body_collection.print_residuals=true;

    use_decoupled_viscosity=parse_args->Is_Value_Set("-use_decoupled_viscosity");
    fluids_parameters.use_coupled_implicit_viscosity=!use_decoupled_viscosity;
    max_dt=(T)(parse_args->Get_Double_Value("-max_dt")*s);
    exact_dt=(T)(parse_args->Get_Double_Value("-dt")*s);
    if(parse_args->Is_Value_Set("-explicit_solid")) implicit_solid=false;
    if(parse_args->Is_Value_Set("-build_surface")) rebuild_surface=true;
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    if(parse_args->Is_Value_Set("-cut_cell")) fluids_parameters.second_order_cut_cell_method=true;
    if(parse_args->Is_Value_Set("-make_ellipse")) make_ellipse=true;
    if(parse_args->Is_Value_Set("-phi")) use_phi=true;
    solids_parameters.use_post_cg_constraints=true;
    if(parse_args->Is_Value_Set("-two_phase")) two_phase=true;
    if(parse_args->Is_Value_Set("-cut_mass")) use_cut_volume=true;
    if(parse_args->Is_Value_Set("-remesh")) remesh=true;
    if(parse_args->Is_Value_Set("-low_order_advection")) use_low_order_advection=true;
    if(use_phi || remesh) solids_parameters.write_static_variables_every_frame=true;

    fluids_parameters.gravity=0;
    fluids_parameters.use_levelset_viscosity=false;
    fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.number_particles_per_cell=16;
    fluids_parameters.write_particles=true;
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;

    switch(test_number){
        case 1:
        case 3:
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,0*m,(T).04*m,0*m,(T).04*m);
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
            fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            break;
        case 2:
            fluids_parameters.grid->Initialize(5,5,10*m,11*m,10*m,11*m);
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
            fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            break;
        case 4:
        case 5:
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,0*m,1*m,0*m,1*m);
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
            fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            break;
        case 6:
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            (*fluids_parameters.grid).Initialize(resolution+1,resolution+1,0,1,(T)0,m);
            if(!parse_args->Is_Value_Set("-viscosity")) fluids_parameters.viscosity=100;
            solid_density=150;
            solid_width=(T)1/3;
            break;
        case 7:
            fluids_parameters.grid->Initialize(resolution+1,resolution*2+1,0*m,(T)1*m,0*m,(T)2*m);
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        case 8:
            fluids_parameters.grid->Initialize(2*resolution+1,resolution+1,0*m,(T)8*m,-(T)2*m,(T)2*m);
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        case 9:
        case 10:
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,0*m,(T)1*m,0*m,(T)1*m);
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
            fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            break;
        case 11:
            fluids_parameters.grid->Initialize(resolution+1,resolution/2*3+1,(T)-.01*m,(T).01*m,(T)-.01*m,(T).02*m);
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        case 12:
            resolution=64;
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,(T)-.5*m,(T).5*m,(T)-.5*m,(T).5*m);
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        case 13:
            fluids_parameters.grid->Initialize(resolution+1,resolution/2*3+1,(T)-.01*m,(T).01*m,(T)-.01*m,(T).02*m);
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

    Add_Rigid_Body_Walls();
    output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number,resolution);
}
//#####################################################################
// Function Add_Rigid_Body_Walls
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Add_Rigid_Body_Walls(const T coefficient_of_restitution,const T coefficient_of_friction,ARRAY<int>* walls_added)
{
    RIGID_BODY_COLLECTION<VECTOR<T,2> >& rigid_body_collection=solid_body_collection.rigid_body_collection;
    VECTOR<T,2> center=fluids_parameters.grid->domain.Center(),size=fluids_parameters.grid->domain.Edge_Lengths();
    int id;

    if(fluids_parameters.domain_walls(0)(0)){
        id=rigid_body_collection.Add_Rigid_Body(this->stream_type,data_directory+"/Rigid_Bodies_2D/ground",size.y*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Friction(coefficient_of_friction);
        rigid_body_collection.rigid_body_particle.X(id)=VECTOR<T,2>(fluids_parameters.grid->domain.min_corner.x,center.y);
        rigid_body_collection.rigid_body_particle.rotation(id)=ROTATION<VECTOR<T,2> >::From_Angle(-(T)pi/2);
        rigid_body_collection.Rigid_Body(id).Set_Name("left wall");
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(fluids_parameters.domain_walls(0)(1)){
        id=rigid_body_collection.Add_Rigid_Body(this->stream_type,data_directory+"/Rigid_Bodies_2D/ground",size.y*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Friction(coefficient_of_friction);
        rigid_body_collection.rigid_body_particle.X(id)=VECTOR<T,2>(fluids_parameters.grid->domain.max_corner.x,center.y);
        rigid_body_collection.rigid_body_particle.rotation(id)=ROTATION<VECTOR<T,2> >::From_Angle((T)pi/2);
        rigid_body_collection.Rigid_Body(id).Set_Name("right wall");
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(fluids_parameters.domain_walls(1)(0)){
        id=rigid_body_collection.Add_Rigid_Body(this->stream_type,data_directory+"/Rigid_Bodies_2D/ground",size.x*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Friction(coefficient_of_friction);
        rigid_body_collection.rigid_body_particle.X(id)=VECTOR<T,2>(center.x,fluids_parameters.grid->domain.min_corner.y);
        rigid_body_collection.Rigid_Body(id).Set_Name("bottom wall");
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(fluids_parameters.domain_walls(1)(1)){
        id=rigid_body_collection.Add_Rigid_Body(this->stream_type,data_directory+"/Rigid_Bodies_2D/ground",size.x*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Friction(coefficient_of_friction);
        rigid_body_collection.rigid_body_particle.X(id)=VECTOR<T,2>(center.x,fluids_parameters.grid->domain.max_corner.y);
        rigid_body_collection.rigid_body_particle.rotation(id)=ROTATION<VECTOR<T,2> >::From_Angle((T)pi);
        rigid_body_collection.Rigid_Body(id).Set_Name("top wall");
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}
}
template<class T> void SURFACE_TENSION<T>::
Parse_Late_Options() {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Initialize_Advection()
{
    PHYSBAM_ASSERT(!fluids_parameters.use_conservative_advection,"cannot do fluid coupling with conservative advection yet");
    fluids_parameters.particle_levelset_evolution->Levelset_Advection(0).
        Use_Semi_Lagrangian_Collidable_Advection(*fluids_parameters.collision_bodies_affecting_fluid,fluids_parameters.collidable_phi_replacement_value,
            fluids_parameters.incompressible->valid_mask);
    fluids_parameters.particle_levelset_evolution->Levelset(0).Set_Collision_Body_List(*fluids_parameters.collision_bodies_affecting_fluid);
    fluids_parameters.incompressible->Use_Semi_Lagrangian_Collidable_Advection(*fluids_parameters.collision_bodies_affecting_fluid);
    fluids_parameters.incompressible->collision_body_list=fluids_parameters.collision_bodies_affecting_fluid;

    fluids_parameters.particle_half_bandwidth=16;
    if(!use_low_order_advection){
        fluids_parameters.particle_levelset_evolution->runge_kutta_order_particles=3;
        fluids_parameters.particle_levelset_evolution->runge_kutta_order_levelset=3;
        fluids_parameters.particle_levelset_evolution->Use_Hamilton_Jacobi_Weno_Advection();
        fluids_parameters.incompressible->advection=new ADVECTION_HAMILTON_JACOBI_ENO<GRID<TV>,T>;
        convection_order=3;
        fluids_parameters.particle_levelset_evolution->Use_Reinitialization();}
    else{
        fluids_parameters.particle_levelset_evolution->runge_kutta_order_particles=1;
        fluids_parameters.particle_levelset_evolution->runge_kutta_order_levelset=1;
        convection_order=1;}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Initialize_Phi()
{
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    if(test_number==1 || test_number==3){
        SPHERE<TV> object(TV((T).02*m,(T).02*m),(T).01*m);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            TV X=it.Location();
            if(make_ellipse) X*=TV((T)1.1,(T).9);
            phi(it.index)=object.Signed_Distance(X);}}
    else if(test_number==2) phi.Fill(1);
    else if(test_number==4 || test_number==5){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            TV dX=it.Location()-TV((T)(.5*m),(T)(.5*m));
            T distance=dX.Magnitude();
            T angle=atan2(dX.y,dX.x);
            T radius=circle_radius+circle_perturbation*cos(oscillation_mode*angle);
            phi(it.index)=distance-radius;}}
    else if(test_number==6) phi.Fill(-1);
    else if(test_number==7){
        SPHERE<TV> object(TV((T).5*m,(T).5*m),(T).2*m);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next())
            phi(it.index)=-object.Signed_Distance(it.Location());}
    else if(test_number==10 || test_number==9){
        SPHERE<TV> object(TV((T).5*m,(T).5*m),(T).25*m);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            phi(it.index)=object.Signed_Distance(it.Location());}}
    else if(test_number==8) Initialize_Sine_Phi();
    else if(test_number==11 || test_number==13){
        SPHERE<TV> object(TV(),(T)1/300*m);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next())
            phi(it.index)=-object.Signed_Distance(it.Location());}
    else if(test_number==12){
        SPHERE<TV> object(TV(),(T).4*m);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next())
            phi(it.index)=-object.Signed_Distance(it.Location());}
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Preprocess_Substep(const T dt,const T time)
{
    current_dt=dt;
    if(SURFACE_TENSION_FORCE<TV>* force=solid_body_collection.deformable_body_collection.template Find_Force<SURFACE_TENSION_FORCE<TV>*>()) force->dt=dt;
    if(LINEAR_POINT_ATTRACTION<TV>* force=solid_body_collection.deformable_body_collection.template Find_Force<LINEAR_POINT_ATTRACTION<TV>*>()) force->dt=dt;

    if(test_number==4 || test_number==5){
        Test_Analytic_Velocity(time);
        Test_Analytic_Pressure(time);}
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Preprocess_Frame(const int frame)
{
    if(fluids_parameters.use_slip){
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).run_self_tests=run_self_tests;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).output_iterators=output_iterators;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_matrix=print_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_each_matrix=print_each_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).use_full_ic=use_full_ic;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_poisson_matrix=print_poisson_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).use_viscous_forces=!use_decoupled_viscosity;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_index_map=print_index_map;}
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Initialize_Velocities()
{
    if(fluids_parameters.use_coupled_implicit_viscosity){
        fluid_collection.incompressible_fluid_collection.viscosity.Resize(fluids_parameters.grid->Domain_Indices(1));
        fluid_collection.incompressible_fluid_collection.viscosity.Fill(fluids_parameters.viscosity);
        fluids_parameters.implicit_viscosity=false;}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    BASE::Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_BOOL& psi_N,const T time)
{
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Initialize_Bodies()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
//    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    switch(test_number){
        case 1: Kang_Circle(false);break;
        case 2: Solid_Circle();break;
        case 3: Kang_Circle(true);break;
        case 4: Oscillating_Circle(false);break;
        case 5: Oscillating_Circle(true);break;
        case 6: FSI_Analytic_Test();break;
        case 7: Kang_Circle(true);break;
        case 8: Sine_Wave();break;
        case 9: Kang_Circle(false);break;
        case 10: Kang_Circle(true);break;
        case 11: Kang_Circle(true);break;
        case 12: Kang_Circle(true);break;
        case 13: Kang_Circle(false);break;
        default: LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    if(test_number!=3) deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    use_pls_evolution_for_structure=false;//use_massless_structure;

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    SEGMENTED_CURVE_2D<T>* curve=rebuild_curve;
    if(!curve) curve=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE_2D<T>*>();
    if(curve){
        SURFACE_TENSION_FORCE<TV>* stf=new SURFACE_TENSION_FORCE<TV>(*curve,fluids_parameters.surface_tension);
        solid_body_collection.Add_Force(stf);
        stf->apply_explicit_forces=true;
        stf->apply_implicit_forces=implicit_solid;
        fluids_parameters.surface_tension=0;}

    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++) solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->compute_half_forces=true;
}
//#####################################################################
// Function Kang_Circle
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Kang_Circle(bool use_surface)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    if(test_number==7){
        fluids_parameters.gravity=(T)(8e-4)*m/(s*s);
        fluids_parameters.density=(T)10000*kg/(m*m);
        fluids_parameters.outside_density=(T)1000*kg/(m*m);
        fluids_parameters.viscosity=1*kg/s;
        fluids_parameters.surface_tension=(T).5;}
    else if(test_number==1 || test_number==3){
        fluids_parameters.gravity=(T)0*m/(s*s);
        fluids_parameters.density=(T)1000*kg/(m*m);
        fluids_parameters.outside_density=(T)1*kg/(m*m);}
    else if(test_number==10 || test_number==9){
        fluids_parameters.gravity=(T)0*m/(s*s);
        fluids_parameters.density=(T)10000*kg/(m*m);
        if(two_phase) fluids_parameters.outside_density=(T)1*kg/(m*m);
        else fluids_parameters.outside_density=(T)0;
        fluids_parameters.viscosity=(T)1*kg/s;
        fluids_parameters.surface_tension=(T)1;}
    else if(test_number==11 || test_number==13){
        two_phase=true;
        if(test_number==13) fluids_parameters.second_order_cut_cell_method=true;
        fluids_parameters.gravity=(T)9.8*m/(s*s);
        fluids_parameters.density=(T)1000*kg/(m*m);
        fluids_parameters.outside_density=(T)1.226*kg/(m*m);
        fluids_parameters.viscosity=(T)1.78e-5*kg/s;
        fluids_parameters.surface_tension=(T).0728;}
    else if(test_number==12){
        /* Fig. 7 rescaled by T_nu & U_sig; maximum time length: T_nu
           Fig. 8 rescaled by T_sig & U_sig; maximum time length: 3*T_sig
        */
        last_frame=1000;
        surface_tension=fluids_parameters.surface_tension;
        T density=1000,D=(T).8;
        T T_sig=D*sqrt(density*D/surface_tension),
            T_nu=T_sig*sqrt(laplace_number);
        LOG::cout<<"T_sig="<<T_sig<<", T_nu="<<T_nu<<std::endl;
        if(use_T_nu) frame_rate=(T)last_frame/T_nu;
        else frame_rate=(T)last_frame/(3*T_sig);
        two_phase=true;
        fluids_parameters.gravity=(T)0*m/(s*s);
        fluids_parameters.density=density*kg/(m*m);
        fluids_parameters.outside_density=density*kg/(m*m);
        // this is dynamic viscosity=kinematic viscosity * density
        if(laplace_number) fluids_parameters.viscosity=sqrt(surface_tension*D*density/laplace_number)*kg/s;
        // this is for infinite La
        else fluids_parameters.viscosity=0*kg/s;}
    if(parse_args->Is_Value_Set("-viscosity")) fluids_parameters.viscosity=(T)(parse_args->Get_Double_Value("-viscosity")*kg/s);
    if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=true;
    fluids_parameters.cfl=(T).9;
//    solids_parameters.write_static_variables_every_frame=true;
    use_massless_structure=use_surface;

    solid_body_collection.Set_CFL_Number(10);

//    fluids_parameters.second_order_cut_cell_method=false;
    fluids_parameters.use_particle_levelset=true;

    if(!rebuild_surface && use_surface){
        SPHERE<TV> object;
        if(test_number==3) object=SPHERE<TV>(TV((T).02*m,(T).02*m),(T).01*m);
        else if(test_number==10) object=SPHERE<TV>(TV((T).5*m,(T).5*m),(T).25*m);
        else if(test_number==7) object=SPHERE<TV>(TV((T).5*m,(T).5*m),(T).2*m);
        else if(test_number==11) object=SPHERE<TV>(TV(),(T)1./300*m);
        else if(test_number==12) object=SPHERE<TV>(TV((T)0*m,(T)0*m),(T).4*m);
        else PHYSBAM_FATAL_ERROR();
        SEGMENTED_CURVE_2D<T>* temp_structure=TESSELLATION::Tessellate_Boundary(object,solid_refinement);
        if(test_number==7 || test_number==11) 
            for(int i=0;i<temp_structure->particles.X.m;i++)
                temp_structure->particles.X(i).x=-temp_structure->particles.X(i).x+2*object.center.x;
        front_tracked_structure=&solids_tests.Copy_And_Add_Structure(*temp_structure);
        solid_body_collection.deformable_body_collection.particles.array_collection->Add_Elements(100*resolution);
        particle_segments.Resize(front_tracked_structure->mesh.elements.Flattened().Max());
        for(int i=0;i<front_tracked_structure->mesh.elements.m;i++){
            particle_segments(front_tracked_structure->mesh.elements(i).x).y=i;
            particle_segments(front_tracked_structure->mesh.elements(i).y).x=i;}
        if(make_ellipse) for(int i=0;i<particle_segments.m;i++) front_tracked_structure->particles.X(i)/=TV((T)1.1,(T).9);
        saved_tracked_particles_X=front_tracked_structure->particles.X.Prefix(particle_segments.m);
        particles.mass.Fill(2*(T)pi*object.radius*fluids_parameters.grid->dX.Max()*fluids_parameters.density/particles.mass.m/100);

        if(1){
            deformable_collisions=new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*front_tracked_structure);//.Get_Boundary_Object());
            deformable_collisions->object.Initialize_Hierarchy();
            Add_To_Fluid_Simulation(*deformable_collisions,false,true);}}

    if(rebuild_surface) Initialize_Surface_Particles(20*resolution);
}
//#####################################################################
// Function Oscillating_Circle
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Oscillating_Circle(bool use_surface)
{
    (*fluids_parameters.grid).Initialize(resolution,resolution,(T)0*m,(T)1*m,(T)0*m,(T)1*m,true);
    fluids_parameters.cfl=(T).9;
//    solids_parameters.write_static_variables_every_frame=true;
    use_massless_structure=use_surface;

    solid_body_collection.Set_CFL_Number(10);

    fluids_parameters.use_particle_levelset=true;

    circle_radius=(T)1/3*m;
    circle_perturbation=(T)(circle_radius*parse_args->Get_Double_Value("-epsilon"));
    oscillation_mode=parse_args->Get_Integer_Value("-oscillation_mode");
    fluids_parameters.surface_tension=(T)2/3*kg*m/(s*s);
    fluids_parameters.density=27*kg/(m*m);
    omega=sqrt(oscillation_mode*(oscillation_mode*oscillation_mode-1)*fluids_parameters.surface_tension/(fluids_parameters.density*circle_radius*circle_radius*circle_radius));
    T period=(T)(2*pi/omega);
    LOG::cout<<"Analytic period of oscillation: "<<period<<std::endl;

    // Forces here
    if(T linear_force=(T)(parse_args->Get_Double_Value("-linear_force"))){
        SPHERE<TV> circle;
        TRIANGULATED_AREA<T>* area=&solids_tests.Copy_And_Add_Structure(*TESSELLATION::Generate_Triangles(circle,solid_refinement));
        LINEAR_POINT_ATTRACTION<TV>* stf=new LINEAR_POINT_ATTRACTION<TV>(*area,TV(.5,.5),linear_force);
        for(int i=0;i<stf->referenced_particles.m;i++){int p=stf->referenced_particles(i);
            TV& X=area->particles.X(p);
            if(X==TV()){X=TV(.5,.5);continue;}
            T angle=atan2(X.y,X.x);
            T scale=circle_radius+circle_perturbation*cos(oscillation_mode*angle);
            X=TV(.5,.5)+X*scale;}
        solid_body_collection.deformable_body_collection.particles.mass.Fill((T)1);
        solid_body_collection.Add_Force(stf);
        stf->apply_explicit_forces=true;
        stf->apply_implicit_forces=implicit_solid;
        particle_segments.Resize(area->particles.X.m);
        if(1){
            deformable_collisions=new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*area);//.Get_Boundary_Object());
            deformable_collisions->object.Initialize_Hierarchy();
            Add_To_Fluid_Simulation(*deformable_collisions,false,true);}}
    else if(use_surface){
        SPHERE<TV> object;
        front_tracked_structure=&solids_tests.Copy_And_Add_Structure(*TESSELLATION::Tessellate_Boundary(object,solid_refinement));
        for(int i=0;i<front_tracked_structure->particles.X.m;i++){
            T angle=(T)(i*2*pi/front_tracked_structure->particles.X.m);
            T radius=circle_radius+circle_perturbation*cos(oscillation_mode*angle);
            front_tracked_structure->particles.X(i)=TV((T).5*m+radius*cos(angle),(T).5*m+radius*sin(angle));}
        solid_body_collection.deformable_body_collection.particles.mass.Fill((T)1);
        solid_body_collection.deformable_body_collection.particles.array_collection->Add_Elements(20*resolution);
        particle_segments.Resize(front_tracked_structure->mesh.elements.Flattened().Max());
        for(int i=0;i<front_tracked_structure->mesh.elements.m;i++){
            particle_segments(front_tracked_structure->mesh.elements(i).x).y=i;
            particle_segments(front_tracked_structure->mesh.elements(i).y).x=i;}
        saved_tracked_particles_X=front_tracked_structure->particles.X.Prefix(particle_segments.m);

        if(1){
            deformable_collisions=new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*front_tracked_structure);//.Get_Boundary_Object());
            deformable_collisions->object.Initialize_Hierarchy();
            Add_To_Fluid_Simulation(*deformable_collisions,false,true);}}

    if(rebuild_surface) Initialize_Surface_Particles(20*resolution);
}
//#####################################################################
// Function Sine_Wave
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Sine_Wave()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    use_massless_structure=true;
    solid_body_collection.Set_CFL_Number(10);
    fluids_parameters.use_particle_levelset=true;

    SPHERE<TV> object(TV((T).02*m,(T).02*m),(T).01*m);
    front_tracked_structure=&solids_tests.Copy_And_Add_Structure(*TESSELLATION::Tessellate_Boundary(object,solid_refinement));
    solid_body_collection.deformable_body_collection.particles.array_collection->Add_Elements(100*resolution);
    particle_segments.Resize(front_tracked_structure->mesh.elements.Flattened().Max());
    for(int i=0;i<front_tracked_structure->mesh.elements.m;i++){
        particle_segments(front_tracked_structure->mesh.elements(i).x).y=i;
        particle_segments(front_tracked_structure->mesh.elements(i).y).x=i;}
    if(make_ellipse) for(int i=0;i<particle_segments.m;i++) front_tracked_structure->particles.X(i)/=TV((T)1.1,(T).9);
    saved_tracked_particles_X=front_tracked_structure->particles.X.Prefix(particle_segments.m);
    particles.mass.Fill(2*(T)pi*object.radius*fluids_parameters.grid->dX.Max()*fluids_parameters.density/particles.mass.m/100);

    deformable_collisions=new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*front_tracked_structure);//.Get_Boundary_Object());
    deformable_collisions->object.Initialize_Hierarchy();
    Add_To_Fluid_Simulation(*deformable_collisions,false,true);
}
template<class T>
struct SINE_DIST:public NONLINEAR_FUNCTION<T(T)>
{
    T X,Y;
    virtual ~SINE_DIST(){}
    virtual T operator()(const T x) const
    {
        return (x-X)+(sin(x)-Y)*cos(x);
    }
};
//#####################################################################
// Function Initialize_Sine_Phi
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Initialize_Sine_Phi()
{
    const GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    SINE_DIST<T> sd;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        TV X=it.Location();
        sd.X=X.x;
        sd.Y=X.y;
        T mx=ITERATIVE_SOLVER<T>().Bisection_Secant_Root(sd,-2,10);
        T dist=(X-TV(mx,sin(mx))).Magnitude();
        T sy=sin(X.x);
        if(X.y<sy) phi(it.index)=dist;
        else{
            T l=grid.domain.min_corner.x-X.x+grid.dX.x*2;
            T r=X.x-grid.domain.max_corner.x+grid.dX.x*2;
            T t=X.y-grid.domain.max_corner.y+grid.dX.y*2;
            T b=-dist;
            phi(it.index)=std::max(std::max(l,r),std::max(t,b));}}
}
//#####################################################################
// Function Test_Analytic_Velocity
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Test_Analytic_Velocity(T time)
{
    T rh=circle_radius;
    T epsilon=circle_perturbation/circle_radius;

    T max_error=0;
    T L1_error=0;
    int cnt=0;
    ARRAY<T,FACE_INDEX<TV::m> > u2(fluid_collection.incompressible_fluid_collection.face_velocities);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
        if(fluids_parameters.particle_levelset_evolution->Levelset(0).Phi(it.Location())>fluids_parameters.grid->dX.Min()*0){
            fluid_collection.incompressible_fluid_collection.face_velocities(it.Full_Index())=0;
            continue;}
        TV dX=it.Location()-TV((T)(.5*m),(T)(.5*m));
        T r=dX.Magnitude();
        T theta=atan2(dX.y,dX.x);
        TV u=-epsilon*omega*rh*pow(r/rh,oscillation_mode-1)*sin(omega*time)*TV(cos((1-oscillation_mode)*theta),sin((1-oscillation_mode)*theta));
        T error=fabs(fluid_collection.incompressible_fluid_collection.face_velocities(it.Full_Index())-u(it.Axis()));
        if(error>max_error) max_error=error;
        L1_error+=error;
        cnt++;
        fluid_collection.incompressible_fluid_collection.face_velocities(it.Full_Index())-=u(it.Axis());}
    T max_u=fabs(epsilon*omega*rh*sin(omega*time));
    if(max_u<1e-10) max_u=(T)1e-10;

    LOG::cout<<"error velocity "<<max_error<<"        "<<L1_error/cnt<<"         "<<max_error/max_u<<"        "<<L1_error/cnt/max_u<<std::endl;
    
    PHYSBAM_DEBUG_WRITE_SUBSTEP("velocity error",0,0);
    fluid_collection.incompressible_fluid_collection.face_velocities=u2;
}
//#####################################################################
// Function Test_Analytic_Pressure
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Test_Analytic_Pressure(T time)
{
    T rh=circle_radius;
    T epsilon=circle_perturbation/circle_radius;

    T max_error=0;
    T L1_error=0;
    int cnt=0;
    ARRAY<T,TV_INT>& p=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).pressure;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
        if(fluids_parameters.particle_levelset_evolution->Levelset(0).Phi(it.Location())>fluids_parameters.grid->dX.Min()*0) continue;
        TV dX=it.Location()-TV((T)(.5*m),(T)(.5*m));
        T r=dX.Magnitude();
        T theta=atan2(dX.y,dX.x);
        T q=fluids_parameters.surface_tension/rh+epsilon/oscillation_mode*omega*omega*rh*rh*pow(r/rh,oscillation_mode)*cos(omega*time)*cos(oscillation_mode*theta);
        T error=fabs(p(it.index)-q);
        if(error>max_error) max_error=error;
        L1_error+=error;
        cnt++;}
    T max_p=fluids_parameters.surface_tension/rh+epsilon/oscillation_mode*omega*omega*rh*rh*fabs(cos(omega*time));
    if(max_p<1e-10) max_p=(T)1e-10;

    LOG::cout<<"error pressure "<<max_error<<"        "<<L1_error/cnt<<"      "<<max_error/max_p<<"        "<<L1_error/cnt/max_p<<std::endl;
}
//#####################################################################
// Function Solid_Circle
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Solid_Circle()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SPHERE<TV> object(TV(),1);
    solids_tests.Copy_And_Add_Structure(*TESSELLATION::Tessellate_Boundary(object,solid_refinement));
    particles.mass+=(T)1/particles.mass.m;
    if(make_ellipse) for(int i=0;i<particles.X.m;i++) particles.X(i)*=TV((T).9,(T)1.1);
    if(T rand=fluids_parameters.viscosity=(T)(parse_args->Get_Double_Value("-rand"))){
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1234);
        ARRAY<TV> dX(particles.X.m);
        random.Fill_Uniform(dX,-rand*m,rand*m);
        particles.X+=dX;}
}
//#####################################################################
// Function Adjust_Phi_With_Objects
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Adjust_Phi_With_Objects(const T time)
{
    return;
    if(front_tracked_structure){
        front_tracked_structure->particles.X.Prefix(particle_segments.m)=saved_tracked_particles_X;
        if(!use_massless_structure) Copy_Front_Tracked_Velocity_From_Fluid();
        front_tracked_structure->particles.X.Prefix(particle_segments.m)+=current_dt*front_tracked_structure->particles.V.Prefix(particle_segments.m);
        Sync_Front_Tracked_Particles_To_Level_Set();
        Sync_Front_Tracked_Particles_To_Level_Set();
        Sync_Front_Tracked_Particles_To_Level_Set();
        saved_tracked_particles_X=front_tracked_structure->particles.X.Prefix(particle_segments.m);}
}
//#####################################################################
// Function Sync_Particle_To_Level_Set
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Sync_Particle_To_Level_Set(int p)
{
    const T_LEVELSET& levelset=fluids_parameters.particle_levelset_evolution->Levelset(0);
    TV& X=front_tracked_structure->particles.X(p);
    X-=levelset.Extended_Phi(X)*levelset.Extended_Normal(X);
}
//#####################################################################
// Function Sync_Front_Tracked_Particles_To_Level_Set
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Sync_Front_Tracked_Particles_To_Level_Set()
{
    return;
    for(int i=0;i<particle_segments.m;i++) Sync_Particle_To_Level_Set(i);
//    for(int i=0;i<front_tracked_structure->mesh.elements.m;i++) Divide_Segment(i);
//    for(int i=particle_segments.m;i>=1;i--) Remove_Particle(i);
    for(int i=0;i<particle_segments.m;i++) solid_body_collection.deformable_body_collection.particles.mass(i)=Compute_New_Mass(i);
}
//#####################################################################
// Function Divide_Segment
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Divide_Segment(int e)
{
    VECTOR<int,2> seg=front_tracked_structure->mesh.elements(e);
    TV XL=front_tracked_structure->particles.X(seg.x),XR=front_tracked_structure->particles.X(seg.y);
    if((XL-XR).Magnitude_Squared()<sqr(fluids_parameters.grid->dX.Min())) return;

    LOG::cout<<"a e "<<front_tracked_structure->mesh.elements<<std::endl;
    LOG::cout<<"a p "<<particle_segments<<std::endl;
    int ne=front_tracked_structure->mesh.elements.m+1;
    int p=particle_segments.Append(VECTOR<int,2>(e,ne));
    LOG::cout<<p<<"  "<<front_tracked_structure->particles.array_collection->Size()<<std::endl;
    PHYSBAM_ASSERT(p<front_tracked_structure->particles.array_collection->Size());
    front_tracked_structure->mesh.elements(e).y=p;
    front_tracked_structure->mesh.elements.Append(VECTOR<int,2>(p,seg.y));
    particle_segments(seg.y).x=ne;

    front_tracked_structure->particles.X(p)=(T).5*(XL+XR);
    TV VL=front_tracked_structure->particles.V(seg.x),VR=front_tracked_structure->particles.V(seg.y);
    front_tracked_structure->particles.V(p)=(T).5*(VL+VR);
    Sync_Particle_To_Level_Set(p);
    LOG::cout<<"b e "<<front_tracked_structure->mesh.elements<<std::endl;
    LOG::cout<<"b p "<<particle_segments<<std::endl;
}
//#####################################################################
// Function Swap_Particles
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Swap_Particles(int p,int r)
{
    exchange(front_tracked_structure->mesh.elements(particle_segments(p).x).y,front_tracked_structure->mesh.elements(particle_segments(r).x).y);
    exchange(front_tracked_structure->mesh.elements(particle_segments(p).y).x,front_tracked_structure->mesh.elements(particle_segments(r).y).x);
    exchange(particle_segments(p),particle_segments(r));
    exchange(front_tracked_structure->particles.X(p),front_tracked_structure->particles.X(r));
    exchange(front_tracked_structure->particles.V(p),front_tracked_structure->particles.V(r));
}
//#####################################################################
// Function Swap_Segments
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Swap_Segments(int e,int f)
{
    exchange(particle_segments(front_tracked_structure->mesh.elements(e).x).y,particle_segments(front_tracked_structure->mesh.elements(f).x).y);
    exchange(particle_segments(front_tracked_structure->mesh.elements(e).y).x,particle_segments(front_tracked_structure->mesh.elements(f).y).x);
    exchange(front_tracked_structure->mesh.elements(e),front_tracked_structure->mesh.elements(f));
}
//#####################################################################
// Function Remove_Particle
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Remove_Particle(int p)
{
    VECTOR<int,2> segs=particle_segments(p);
    TV XL=front_tracked_structure->particles.X(front_tracked_structure->mesh.elements(particle_segments(p).x).x);
    TV XC=front_tracked_structure->particles.X(p);
    TV XR=front_tracked_structure->particles.X(front_tracked_structure->mesh.elements(particle_segments(p).y).y);
    if((XL-XC).Magnitude()+(XC-XR).Magnitude()>(T).5*fluids_parameters.grid->dX.Min()) return;

    LOG::cout<<"r e "<<front_tracked_structure->mesh.elements<<std::endl;
    LOG::cout<<"r p "<<particle_segments<<std::endl;
    Swap_Particles(p,particle_segments.m);
    LOG::cout<<"s e "<<front_tracked_structure->mesh.elements<<std::endl;
    LOG::cout<<"s p "<<particle_segments<<std::endl;
    p=particle_segments.m;
    segs=particle_segments(p);
    Swap_Segments(segs.y,front_tracked_structure->mesh.elements.m);
    LOG::cout<<"t e "<<front_tracked_structure->mesh.elements<<std::endl;
    LOG::cout<<"t p "<<particle_segments<<std::endl;
    segs.y=front_tracked_structure->mesh.elements.m;
    front_tracked_structure->mesh.elements(segs.x).y=front_tracked_structure->mesh.elements(segs.y).y;
    particle_segments(front_tracked_structure->mesh.elements(segs.x).y).x=segs.x;
    particle_segments.Pop();
    front_tracked_structure->mesh.elements.Pop();
    LOG::cout<<"u e "<<front_tracked_structure->mesh.elements<<std::endl;
    LOG::cout<<"u p "<<particle_segments<<std::endl;
}
//#####################################################################
// Function Compute_New_Mass
//#####################################################################
template<class T> T SURFACE_TENSION<T>::
Compute_New_Mass(int p)
{
    VECTOR<int,2> segs=particle_segments(p);
    TV XL=front_tracked_structure->particles.X(front_tracked_structure->mesh.elements(particle_segments(p).x).x);
    TV XC=front_tracked_structure->particles.X(p);
    TV XR=front_tracked_structure->particles.X(front_tracked_structure->mesh.elements(particle_segments(p).y).y);
    return (T).5*((XL-XC).Magnitude()+(XC-XR).Magnitude())*fluids_parameters.grid->dX.Min()*fluids_parameters.density;
}
//#####################################################################
// Function Copy_Front_Tracked_Velocity_From_Fluid
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Copy_Front_Tracked_Velocity_From_Fluid()
{
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    LINEAR_INTERPOLATION_MAC<TV,T> interp(*fluids_parameters.grid);

    ARRAY_VIEW<TV>& X=front_tracked_structure->particles.X,&V=front_tracked_structure->particles.V;
    for(int i=0;i<particle_segments.m;i++) V(i)=interp.Clamped_To_Array(face_velocities,X(i));
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Limit_Dt(T& dt,const T time)
{
    if(max_dt && dt>max_dt) dt=max_dt;
    if(exact_dt) dt=exact_dt;
    if(test_number==12){
        T dx=fluids_parameters.grid->dX.Min();
        T dt_surface=(T)(sqrt(dx*fluids_parameters.density/(pi*surface_tension))*dx);
        dt=min(dt,dt_surface);}
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Limit_Solids_Dt(T& dt,const T time)
{
    if(max_dt && dt>max_dt) dt=max_dt;
    if(exact_dt) dt=exact_dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);
    FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),frame));
    FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),frame),debug_particles);
    debug_particles.array_collection->Delete_All_Elements();
}
//#####################################################################
// Function Initialize_Surface_Particles
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Initialize_Surface_Particles(int number)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    number_surface_particles=number;
    PHYSBAM_ASSERT(particles.array_collection->Size()==0);
    particles.array_collection->Add_Elements(number_surface_particles);
    rebuild_curve=SEGMENTED_CURVE_2D<T>::Create(particles);
    free_particles=FREE_PARTICLES<TV>::Create(particles);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(free_particles);
    free_particles->nodes=IDENTITY_ARRAY<>(number_surface_particles);
}
//#####################################################################
// Function Rebuild_Surface
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Rebuild_Surface()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    for(int i=0;i<number_surface_particles;i++)
        particles.X(i)=TV();

    fluid_interpolation_entries.Remove_All();
    solid_interpolation_entries.Remove_All();
    rebuild_curve->mesh.elements.Remove_All();
    HASHTABLE<TV_INT,ARRAY<int> > used_cells;
    int next=1;
    const ARRAY<T,TV_INT>& phi=fluids_parameters.particle_levelset_evolution->Levelset(0).phi;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face=it.Full_Index();
        TV_INT cell1=face.First_Cell_Index(),cell2=face.Second_Cell_Index();
        T phi1=phi(cell1),phi2=phi(cell2);
        if((phi1>0)==(phi2>0)) continue;
        T theta=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
        TV X1=fluids_parameters.grid->X(cell1),X2=fluids_parameters.grid->X(cell2);
        TV X=(1-theta)*X1+theta*X2;

        typename MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::ENTRY entry={face,phi1<=0?1:2,X};
        fluid_interpolation_entries.Append(entry);

        PHYSBAM_ASSERT(next<number_surface_particles+1);
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));
        particles.X(next)=X;
        solid_interpolation_entries.Append(next);
        used_cells.Get_Or_Insert(cell2).Append(next);
        cell2(1-face.axis)++;
        used_cells.Get_Or_Insert(cell2).Append(next);
        next++;}

    for(typename HASHTABLE<TV_INT,ARRAY<int> >::ITERATOR it(used_cells);it.Valid();it.Next()){
        ARRAY<int>& array=it.Data();
        PHYSBAM_ASSERT(array.m==2);
        rebuild_curve->mesh.elements.Append(VECTOR<int,2>(array(0),array(1)));
        Add_Debug_Particle((rebuild_curve->particles.X(array(0))+rebuild_curve->particles.X(array(1)))/2,VECTOR<T,3>(1,0,0));}
    PHYSBAM_ASSERT(rebuild_curve->particles.X.Get_Array_Pointer()==particles.X.Get_Array_Pointer());
}
//#####################################################################
// Function Substitute_Coupling_Matrices
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Substitute_Coupling_Matrices(KRYLOV_SYSTEM_BASE<T>& coupled_system,T dt,T current_velocity_time,T current_position_time,bool velocity_update,bool leakproof_solve)
{
    if(use_massless_structure){
        SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>& system=dynamic_cast<SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>&>(coupled_system);
        if(!use_phi){
            if(remesh){
                T dx=fluids_parameters.grid->dX.Min();
                ARRAY_VIEW<TV>& X(front_tracked_structure->particles.X);
                ARRAY<TV> copy_X;
                copy_X=front_tracked_structure->particles.X.Prefix(front_tracked_structure->mesh.elements.Flattened().Max());
                X.Fill(TV());
                front_tracked_structure->mesh.elements.Remove_All();
                X(0)=copy_X(0);
                for(int i=0,j=0;j<copy_X.m;j++){
                    int j1=j%copy_X.m+1;
                    T full_length=(X(i)-copy_X(j1)).Magnitude();
                    if(full_length>(T)1.5*dx){
                        int j0=j-1,j2=j1%copy_X.m+1;if(j0==0) j0=copy_X.m;
                        X(++i)=(-copy_X(j0)+(T)9*copy_X(j)+(T)9*copy_X(j1)-copy_X(j2))/(T)16;
                        front_tracked_structure->mesh.elements.Append(TV_INT(i-1,i));
                        if(j1!=1){
                            X(++i)=copy_X(j1);
                            front_tracked_structure->mesh.elements.Append(TV_INT(i-1,i));}
                        else front_tracked_structure->mesh.elements.Append(TV_INT(i,1));}
                    else if(full_length<(T).5*dx){
                        if(j1==1){
                            front_tracked_structure->mesh.elements.Last()=TV_INT(i-1,1);
                            X(i)=TV();}}
                    else{
                        if(j1!=1){
                            X(++i)=copy_X(j1);
                            front_tracked_structure->mesh.elements.Append(TV_INT(i-1,i));}
                        else front_tracked_structure->mesh.elements.Append(TV_INT(i,1));}}}
            fsi=new FLUID_TO_SOLID_INTERPOLATION_CUT<TV>(system.index_map,*front_tracked_structure,fluids_parameters.density);}
        else{
            FLUID_TO_SOLID_INTERPOLATION_PHI<TV>* local_fsi=new FLUID_TO_SOLID_INTERPOLATION_PHI<TV>(system.index_map,fluids_parameters.particle_levelset_evolution->Levelset(0).phi,*front_tracked_structure,fluids_parameters.density);
            local_fsi->cut_order=parse_args->Get_Integer_Value("-cut_order");
            local_fsi->Setup_Mesh();
            fsi=dynamic_cast<FLUID_TO_SOLID_INTERPOLATION_CUT<TV>*>(local_fsi);}
        system.fluid_to_solid_interpolation=fsi;
        fsi->outside_density=fluids_parameters.outside_density;
        fsi->use_cut_volume=use_cut_volume;
        MATRIX_FLUID_GRADIENT_CUT<TV>* gradient=new MATRIX_FLUID_GRADIENT_CUT<TV>(system.index_map);
        system.fluid_gradient=gradient;    
        fsi->fluid_mass=&system.fluid_mass;
        fsi->gradient=gradient;
        SURFACE_TENSION_FORCE<TV>* force=solid_body_collection.deformable_body_collection.template Find_Force<SURFACE_TENSION_FORCE<TV>*>();
        if(!force) return;
        force->Update_Position_Based_State(current_position_time,true);}

    if(rebuild_surface){
        SURFACE_TENSION_FORCE<TV>* force=solid_body_collection.deformable_body_collection.template Find_Force<SURFACE_TENSION_FORCE<TV>*>();
        if(!force) return;

        SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>& system=dynamic_cast<SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>&>(coupled_system);
        MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>* fi=new MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>(system.index_map);
        MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>* si=new MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>(system.index_map.iterator_info);
        system.fluid_interpolation=fi;
        system.solid_interpolation=si;
        fi->entries=fluid_interpolation_entries;
        si->entries=solid_interpolation_entries;

        force->Update_Position_Based_State(current_position_time,true);

        PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle rebuild",0,0);}
}
//#####################################################################
// Function Advance_One_Time_Step_Begin_Callback
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Advance_One_Time_Step_Begin_Callback(const T dt,const T time)
{
    if(rebuild_curve) Rebuild_Surface();
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Update_Time_Varying_Material_Properties(const T time)
{
    if(!time) return;
    if(rebuild_curve) Rebuild_Surface();
    if((test_number==3 || test_number==5 || test_number==10 || test_number==7 || test_number==11 || test_number==12) && !use_phi){
        ARRAY<T,TV_INT>& phi=fluids_parameters.particle_levelset_evolution->Levelset(0).phi;
        LEVELSET_MAKER_UNIFORM_2D<T>::Compute_Level_Set(*front_tracked_structure,*fluids_parameters.grid,3,phi);}
}
//#####################################################################
// Function Store_Debug_Particles
//#####################################################################
template<class T> GEOMETRY_PARTICLES<VECTOR<T,2> >* SURFACE_TENSION<T>::
Store_Debug_Particles(GEOMETRY_PARTICLES<TV>* particle)
{
    static GEOMETRY_PARTICLES<TV>* stored_particles=0;
    GEOMETRY_PARTICLES<TV>* tmp=stored_particles;
    if(particle) stored_particles=particle;
    return tmp;
}
//#####################################################################
// Function Add_Debug_Particle
//#####################################################################
template<class TV> void PhysBAM::
Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color)
{
    typedef typename TV::SCALAR T;
    GEOMETRY_PARTICLES<TV>* particles=(GEOMETRY_PARTICLES<TV>*)SURFACE_TENSION<T>::Store_Debug_Particles();
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles->array_collection->template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    int p=particles->array_collection->Add_Element();
    particles->X(p)=X;
    (*color_attribute)(p)=color;
}
//#####################################################################
// Function FSI_Analytic_Test
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
FSI_Analytic_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    fluids_parameters.collision_bodies_affecting_fluid->use_collision_face_neighbors=true;
    T solid_gravity=(T)9.8*m/(s*s);
    fluids_parameters.surface_tension=0;
    debug_particles.array_collection->template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);

    fluids_parameters.gravity=(T)9.8*m;
    fluids_parameters.density=(T)100/(m*m);
    fluids_parameters.domain_walls[0][0]=true;
    fluids_parameters.domain_walls[0][1]=true;
    fluids_parameters.domain_walls[1][0]=false;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.second_order_cut_cell_method=true;
    fluids_parameters.use_levelset_viscosity=true;
    PHYSBAM_ASSERT(fluids_parameters.use_slip);
    solid_body_collection.Set_CFL_Number(10);

    for(int b=0;b<2;b++){
        Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(b));
        rigid_body_collection.Rigid_Body(b).is_static=true;}

    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Analytic_Box(m*TV(solid_width,(T)1.2),5);
    rigid_body.Set_Frame(FRAME<TV>(m*TV((T).5,.5)));
    rigid_body.Set_Coefficient_Of_Restitution((T)0);
    rigid_body.Set_Mass(solid_density*solid_width);

    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity*m));
    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body);

    T solid_mass=solid_body_collection.rigid_body_collection.rigid_body_particle.mass(2);
    T rho=fluids_parameters.density;
    TV size=fluids_parameters.grid->domain.Edge_Lengths();
    size.x=(size.x-solid_width)/2;
    analytic_solution=-(solid_mass*solid_gravity+rho*size.x*size.y*fluids_parameters.gravity)*size.x/(2*size.y*fluids_parameters.viscosity);
    LOG::cout<<"analytic_solution "<<analytic_solution<<std::endl;

    FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),0));
    FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),0),debug_particles);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Postprocess_Frame(const int frame)
{
    if(test_number==6){
        T v=fluid_collection.incompressible_fluid_collection.face_velocities(FACE_INDEX<2>(2,fluids_parameters.grid->counts/2));
        if(solid_body_collection.rigid_body_collection.rigid_body_particle.X.m>=3) v=solid_body_collection.rigid_body_collection.rigid_body_particle.twist(2).linear.y;
        LOG::cout<<"middle-velocity "<<v<<"   error from analytic solution "<<(v/analytic_solution-1)<<std::endl;}
    if(test_number==10 || test_number==9){ 
        ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
        int count=0;
        T mx=0,sum=0;
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            TV_INT cell1=it.First_Cell_Index(),cell2=it.Second_Cell_Index();
            T phi1=phi(cell1),phi2=phi(cell2);
            if(phi1<=0 || phi2<=0){
                T v=fabs(fluid_collection.incompressible_fluid_collection.face_velocities(it.Full_Index()));
                mx=std::max(v,mx);
                sum+=v;
                count++;}}
        LOG::cout<<"Parasitic current error:  L1="<<sum/count<<"  Linf="<<mx<<std::endl;}
    if(test_number==12){
        T D=(T).8,sum=0;
        T U_sig=sqrt(surface_tension/(fluids_parameters.density*D));
        int num=0;
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            T v=fluid_collection.incompressible_fluid_collection.face_velocities(it.Full_Index());
            sum+=sqr(v);
            num++;}
        LOG::cout<<"RMS velocity="<<sqrt(sum/(T)num)/U_sig<<std::endl;}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T> void SURFACE_TENSION<T>::
Postprocess_Substep(const T dt,const T time)
{
    if(test_number==6 && solid_body_collection.rigid_body_collection.rigid_body_particle.X.m>=3) solid_body_collection.rigid_body_collection.rigid_body_particle.X(2)=TV(.5,.5)*m;
    if(test_number==4 || test_number==5){
        // location test
        TV p((T)((.5+circle_radius)*m),(T)(.5*m));
        T value=fluids_parameters.particle_levelset_evolution->Levelset(0).Phi(p);
        for(int i=0;i<5;i++){
            p(0)-=value;
            value=fluids_parameters.particle_levelset_evolution->Levelset(0).Phi(p);}
        LOG::cout<<"Interface Location: "<<p<<std::endl;
        LOG::cout<<"Analytic Interface Location: "<<TV(.5,.5)*m+circle_radius+circle_perturbation*cos(omega*time)*TV(1,0)<<std::endl;

        // velocity test
        p=TV(.75,.5)*m;
        LINEAR_INTERPOLATION_MAC<TV,T> interp(*fluids_parameters.grid);
        LOG::cout<<"Vel at [.75,.5]: "<<interp.Clamped_To_Array(fluid_collection.incompressible_fluid_collection.face_velocities,p)<<std::endl;
        LOG::cout<<"Analytic Vel at [.75,.5]: "<<-omega*circle_perturbation*pow(T(.25)/circle_radius,oscillation_mode-1)*sin(omega*time)*TV(1,0)<<std::endl;
    }
}
//#####################################################################
// Function Debug_Particle_Set_Attribute
//#####################################################################
template<class TV,class ATTR> void PhysBAM::
Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr)
{
    typedef typename TV::SCALAR T;
    GEOMETRY_PARTICLES<TV>* particles=(GEOMETRY_PARTICLES<TV>*)SURFACE_TENSION<T>::Store_Debug_Particles();
    ARRAY_VIEW<ATTR>* attribute=particles->array_collection->template Get_Array<ATTR>(id);
    attribute->Last()=attr;
}
template class SURFACE_TENSION<float>;
template void PhysBAM::Add_Debug_Particle<VECTOR<float,1> >(VECTOR<float,1> const&,VECTOR<float,3> const&);
template void PhysBAM::Add_Debug_Particle<VECTOR<float,2> >(VECTOR<float,2> const&,VECTOR<float,3> const&);
template void PhysBAM::Add_Debug_Particle<VECTOR<float,3> >(VECTOR<float,3> const&,VECTOR<float,3> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<float,1>,VECTOR<float,1> >(ATTRIBUTE_ID,VECTOR<float,1> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<float,2>,VECTOR<float,2> >(ATTRIBUTE_ID,VECTOR<float,2> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<float,3>,VECTOR<float,3> >(ATTRIBUTE_ID,VECTOR<float,3> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SURFACE_TENSION<double>;
template void PhysBAM::Add_Debug_Particle<VECTOR<double,1> >(VECTOR<double,1> const&,VECTOR<double,3> const&);
template void PhysBAM::Add_Debug_Particle<VECTOR<double,2> >(VECTOR<double,2> const&,VECTOR<double,3> const&);
template void PhysBAM::Add_Debug_Particle<VECTOR<double,3> >(VECTOR<double,3> const&,VECTOR<double,3> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<double,1>,VECTOR<double,1> >(ATTRIBUTE_ID,VECTOR<double,1> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<double,2>,VECTOR<double,2> >(ATTRIBUTE_ID,VECTOR<double,2> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<double,3>,VECTOR<double,3> >(ATTRIBUTE_ID,VECTOR<double,3> const&);
#endif
