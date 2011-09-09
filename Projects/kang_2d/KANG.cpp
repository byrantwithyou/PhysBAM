#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLECTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_CUT.h>
#include "KANG.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4355) // 'this' : used in base member initializer list
#endif
template<class T> KANG<T>::
KANG(const STREAM_TYPE stream_type)
    :BASE(stream_type,1),solids_tests(*this,solid_body_collection),output_iterators(false),max_dt(0),exact_dt(0),
    circle_radius(0),circle_perturbation(0),oscillation_mode(0),make_ellipse(false),m(1),s(1),kg(1),
    omega(0),laplace_number(0)
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
template<class T> KANG<T>::
~KANG()
{
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T> void KANG<T>::
Register_Options()
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-cg_iterations",3000);
    parse_args->Add_Double_Argument("-solve_tolerance",1e-14);
    parse_args->Add_Double_Argument("-viscosity",(T)0);
    parse_args->Add_Option_Argument("-test_system");
    parse_args->Add_Option_Argument("-print_matrix");
    parse_args->Add_Option_Argument("-output_iterators");
    parse_args->Add_Option_Argument("-no_preconditioner");
    parse_args->Add_Option_Argument("-preconditioner");
    parse_args->Add_Double_Argument("-max_dt",0);
    parse_args->Add_Double_Argument("-dt",0);
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Double_Argument("-surface_tension",.0728);
    parse_args->Add_Integer_Argument("-oscillation_mode",2);
    parse_args->Add_Option_Argument("-make_ellipse");
    parse_args->Add_Double_Argument("-epsilon",.05);
    parse_args->Add_Double_Argument("-m",1,"length unit");
    parse_args->Add_Double_Argument("-s",1,"time unit");
    parse_args->Add_Double_Argument("-kg",1,"mass unit");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T> void KANG<T>::
Parse_Options()
{
    BASE::Parse_Options();
    last_frame=100;
    frame_rate=24;

    kg=(T)parse_args->Get_Double_Value("-kg");
    m=(T)parse_args->Get_Double_Value("-m");
    s=(T)parse_args->Get_Double_Value("-s");
    frame_rate/=s;

    fluids_parameters.cfl=(T).9;
    fluids_parameters.confinement_parameter=(T).04;
    fluids_parameters.rho_bottom=1;
    fluids_parameters.rho_top=(T).65;
    fluids_parameters.density_buoyancy_constant=fluids_parameters.temperature_buoyancy_constant=0;
    fluids_parameters.temperature_container.Set_Cooling_Constant(0);
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    use_kang=true;

    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;
    fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;
    fluids_parameters.incompressible_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    fluids_parameters.incompressible_tolerance=parse_args->Get_Double_Value("-solve_tolerance");

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=200;
    solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;
    fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=true;

    fluids_parameters.use_slip=true;
    test_system=parse_args->Is_Value_Set("-test_system");
    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    output_iterators=parse_args->Is_Value_Set("-output_iterators");
    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-no_preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=false;
    two_phase=true;

    fluids_parameters.viscosity=(T)(parse_args->Get_Double_Value("-viscosity")*kg/s);
    fluids_parameters.surface_tension=(T)(parse_args->Get_Double_Value("-surface_tension")*kg*m/(s*s));

    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.store_particle_ids=true;
    fluids_parameters.removed_positive_particle_buoyancy_constant=0;

    max_dt=(T)(parse_args->Get_Double_Value("-max_dt")*s);
    exact_dt=(T)(parse_args->Get_Double_Value("-dt")*s);
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    if(parse_args->Is_Value_Set("-make_ellipse")) make_ellipse=true;
    solids_parameters.use_post_cg_constraints=true;

    fluids_parameters.gravity=0;
    fluids_parameters.use_levelset_viscosity=false;
    fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.number_particles_per_cell=16;
    fluids_parameters.write_particles=true;
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;

    switch(test_number){
        case 1:
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,0*m,(T).04*m,0*m,(T).04*m);
            fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=false;
            fluids_parameters.domain_walls[2][1]=false;fluids_parameters.domain_walls[2][2]=false;
            break;
        case 2:
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,0*m,1*m,0*m,1*m);
            fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=false;
            fluids_parameters.domain_walls[2][1]=false;fluids_parameters.domain_walls[2][2]=false;
            break;
        case 3:
            fluids_parameters.grid->Initialize(resolution+1,resolution*2+1,0*m,(T)1*m,0*m,(T)2*m);
            fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;
            fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=true;
            break;
        case 4:
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,-m,m,-m,m);
            fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=true;
            fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=true;
            break;
        case 5:
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,-m,m,-m,m);
            fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;
            fluids_parameters.domain_walls[2][1]=false;fluids_parameters.domain_walls[2][2]=false;
            break;
        case 6:
            fluids_parameters.grid->Initialize(resolution+1,resolution+1,-m,m,-m,m);
            fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;
            fluids_parameters.domain_walls[2][1]=false;fluids_parameters.domain_walls[2][2]=false;
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

    output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number,resolution);
}
template<class T> void KANG<T>::
Parse_Late_Options() {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
template<class T> void KANG<T>::
Initialize_Advection()
{
    PHYSBAM_ASSERT(!fluids_parameters.use_conservative_advection,"cannot do fluid coupling with conservative advection yet");
    fluids_parameters.particle_levelset_evolution->Levelset_Advection(1).
        Use_Semi_Lagrangian_Collidable_Advection(*fluids_parameters.collision_bodies_affecting_fluid,fluids_parameters.collidable_phi_replacement_value,
            fluids_parameters.incompressible->valid_mask);
    fluids_parameters.particle_levelset_evolution->Levelset(1).Set_Collision_Body_List(*fluids_parameters.collision_bodies_affecting_fluid);
    fluids_parameters.incompressible->Use_Semi_Lagrangian_Collidable_Advection(*fluids_parameters.collision_bodies_affecting_fluid);
    fluids_parameters.incompressible->collision_body_list=fluids_parameters.collision_bodies_affecting_fluid;

    fluids_parameters.particle_half_bandwidth=16;
    fluids_parameters.particle_levelset_evolution->runge_kutta_order_particles=3;
    fluids_parameters.particle_levelset_evolution->runge_kutta_order_levelset=3;
    fluids_parameters.particle_levelset_evolution->Use_Hamilton_Jacobi_Weno_Advection();
    fluids_parameters.incompressible->advection=new ADVECTION_HAMILTON_JACOBI_ENO<GRID<TV>,T>;
    convection_order=3;
    fluids_parameters.particle_levelset_evolution->Use_Reinitialization();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class T> void KANG<T>::
Initialize_Phi()
{
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    if(test_number==1 || test_number==4){
        SPHERE<TV> object(TV((T).02*m,(T).02*m),(T).01*m);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            TV X=it.Location();
            if(make_ellipse) X*=TV((T)1.1,(T).9);
            phi(it.index)=object.Signed_Distance(X);}}
    else if(test_number==2){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            TV dX=it.Location()-TV((T)(.5*m),(T)(.5*m));
            T distance=dX.Magnitude();
            T angle=atan2(dX.y,dX.x);
            T radius=circle_radius+circle_perturbation*cos(oscillation_mode*angle);
            phi(it.index)=distance-radius;}}
    else if(test_number==3){
        SPHERE<TV> object(TV((T).5*m,(T).5*m),(T).2*m);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next())
            phi(it.index)=-object.Signed_Distance(it.Location());}
    else if(test_number==5 || test_number==6){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next())
            phi(it.index)=it.Location().x;}
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class T> void KANG<T>::
Preprocess_Substep(const T dt,const T time)
{
    if(test_number==2){
        Test_Analytic_Velocity(time);
        Test_Analytic_Pressure(time);}
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class T> void KANG<T>::
Preprocess_Frame(const int frame)
{
    if(test_number==6 && frame==1){
        ARRAY<T,FACE_INDEX<TV::m> >& u=fluid_collection.incompressible_fluid_collection.face_velocities;
        T mu0=fluids_parameters.viscosity;
        T mu1=fluids_parameters.outside_viscosity;
        T v0=(mu1-mu0)/(mu1+mu0);
        T dx=fluids_parameters.grid->dX(1);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,0,GRID<TV>::WHOLE_REGION,0,2);it.Valid();it.Next()){
            TV x=it.Location();
            if(x.x>0) u(it.Full_Index())=(1-v0)/(1+dx/2)*x.x+v0;
            else u(it.Full_Index())=(1+v0)/(1+dx/2)*x.x+v0;}}
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
template<class T> void KANG<T>::
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
template<class T> void KANG<T>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    BASE::Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T> void KANG<T>::
Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_BOOL& psi_N,const T time)
{
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void KANG<T>::
Initialize_Bodies()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    switch(test_number){
        case 1: Kang_Circle();break;
        case 2: Oscillating_Circle();break;
        case 3: Kang_Circle();break;
        case 4: Poisson_Test();break;
        case 5: Poiseuille_Flow_Test();break;
        case 6: Couette_Flow_Test();break;
        default: LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=false;

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    use_pls_evolution_for_structure=false;//use_massless_structure;

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    for(int i=1;i<=solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
    for(int k=1;k<=solid_body_collection.deformable_body_collection.deformables_forces.m;k++)
        solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
    for(int i=1;i<=solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->compute_half_forces=true;
}
//#####################################################################
// Function Kang_Circle
//#####################################################################
template<class T> void KANG<T>::
Kang_Circle()
{
    if(test_number==3){
        fluids_parameters.gravity=(T)(8e-4)*m/(s*s);
        fluids_parameters.density=(T)10000*kg/(m*m);
        fluids_parameters.outside_density=(T)1000*kg/(m*m);
        fluids_parameters.viscosity=1*kg/s;
        fluids_parameters.surface_tension=(T).5;}
    else if(test_number==1){
        fluids_parameters.gravity=(T)0*m/(s*s);
        fluids_parameters.density=(T)1000*kg/(m*m);
        fluids_parameters.outside_density=(T)1*kg/(m*m);}
    if(parse_args->Is_Value_Set("-viscosity")) fluids_parameters.viscosity=(T)(parse_args->Get_Double_Value("-viscosity")*kg/s);
    fluids_parameters.cfl=(T).9;

    solid_body_collection.Set_CFL_Number(10);

    fluids_parameters.use_particle_levelset=true;
}
//#####################################################################
// Function Poisson_Test
//#####################################################################
template<class T> void KANG<T>::
Poisson_Test()
{
    fluids_parameters.gravity=0;
    fluids_parameters.density=1;
    fluids_parameters.outside_density=1;
    fluids_parameters.viscosity=0;
    fluids_parameters.surface_tension=0;
    fluids_parameters.use_particle_levelset=true;
}
//#####################################################################
// Function Poiseuille_Flow_Test
//#####################################################################
template<class T> void KANG<T>::
Poiseuille_Flow_Test()
{
    fluids_parameters.gravity=(T)9.8*m/(s*s);
    fluids_parameters.density=(T)1*kg/(m*m);
    fluids_parameters.outside_density=(T)2*kg/(m*m);
    fluids_parameters.viscosity=(T)1*kg/s;
    fluids_parameters.outside_viscosity=(T)3*kg/s;
    fluids_parameters.surface_tension=0;
    fluids_parameters.use_particle_levelset=true;
}
//#####################################################################
// Function Couette_Flow_Test
//#####################################################################
template<class T> void KANG<T>::
Couette_Flow_Test()
{
    fluids_parameters.gravity=(T)0*m/(s*s);
    fluids_parameters.density=(T)1*kg/(m*m);
    fluids_parameters.outside_density=(T)1*kg/(m*m);
    fluids_parameters.viscosity=(T)1*kg/s;
    fluids_parameters.outside_viscosity=(T)2*kg/s;
    fluids_parameters.surface_tension=0;
    fluids_parameters.use_particle_levelset=true;
}
//#####################################################################
// Function Oscillating_Circle
//#####################################################################
template<class T> void KANG<T>::
Oscillating_Circle()
{
    (*fluids_parameters.grid).Initialize(resolution,resolution,(T)0*m,(T)1*m,(T)0*m,(T)1*m,true);
    fluids_parameters.cfl=(T).9;

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
}
//#####################################################################
// Function Test_Analytic_Velocity
//#####################################################################
template<class T> void KANG<T>::
Test_Analytic_Velocity(T time)
{
    T rh=circle_radius;
    T epsilon=circle_perturbation/circle_radius;

    T max_error=0;
    T L1_error=0;
    int cnt=0;
    ARRAY<T,FACE_INDEX<TV::m> > u2(fluid_collection.incompressible_fluid_collection.face_velocities);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
        if(fluids_parameters.particle_levelset_evolution->Levelset(1).Phi(it.Location())>fluids_parameters.grid->dX.Min()*0){
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
template<class T> void KANG<T>::
Test_Analytic_Pressure(T time)
{
#if 0
    T rh=circle_radius;
    T epsilon=circle_perturbation/circle_radius;

    T max_error=0;
    T L1_error=0;
    int cnt=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
        if(fluids_parameters.particle_levelset_evolution->Levelset(1).Phi(it.Location())>fluids_parameters.grid->dX.Min()*0) continue;
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
#endif
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
template<class T> void KANG<T>::
Limit_Dt(T& dt,const T time)
{
    if(max_dt && dt>max_dt) dt=max_dt;
    if(exact_dt) dt=exact_dt;
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
template<class T> void KANG<T>::
Limit_Solids_Dt(T& dt,const T time)
{
    if(max_dt && dt>max_dt) dt=max_dt;
    if(exact_dt) dt=exact_dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void KANG<T>::
Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);
    FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),frame));
    FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),frame),debug_particles);
    debug_particles.array_collection->Delete_All_Elements();
}
//#####################################################################
// Function Advance_One_Time_Step_Begin_Callback
//#####################################################################
template<class T> void KANG<T>::
Advance_One_Time_Step_Begin_Callback(const T dt,const T time)
{
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class T> void KANG<T>::
Update_Time_Varying_Material_Properties(const T time)
{
}
//#####################################################################
// Function Store_Debug_Particles
//#####################################################################
template<class T> GEOMETRY_PARTICLES<VECTOR<T,2> >* KANG<T>::
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
    GEOMETRY_PARTICLES<TV>* particles=(GEOMETRY_PARTICLES<TV>*)KANG<T>::Store_Debug_Particles();
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles->array_collection->template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    int p=particles->array_collection->Add_Element();
    particles->X(p)=X;
    (*color_attribute)(p)=color;
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class T> void KANG<T>::
Postprocess_Frame(const int frame)
{
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T> void KANG<T>::
Postprocess_Substep(const T dt,const T time)
{
    if(test_number==2){
        // location test
        TV p((T)((.5+circle_radius)*m),(T)(.5*m));
        T value=fluids_parameters.particle_levelset_evolution->Levelset(1).Phi(p);
        for(int i=0;i<5;i++){
            p(1)-=value;
            value=fluids_parameters.particle_levelset_evolution->Levelset(1).Phi(p);}
        LOG::cout<<"Interface Location: "<<p<<std::endl;
        LOG::cout<<"Analytic Interface Location: "<<TV(.5,.5)*m+circle_radius+circle_perturbation*cos(omega*time)*TV(1,0)<<std::endl;

        // velocity test
        p=TV(.75,.5)*m;
        LINEAR_INTERPOLATION_MAC<TV,T> interp(*fluids_parameters.grid);
        LOG::cout<<"Vel at [.75,.5]: "<<interp.Clamped_To_Array(fluid_collection.incompressible_fluid_collection.face_velocities,p)<<std::endl;
        LOG::cout<<"Analytic Vel at [.75,.5]: "<<-omega*circle_perturbation*pow(T(.25)/circle_radius,oscillation_mode-1)*sin(omega*time)*TV(1,0)<<std::endl;
    }
    if(test_number==5 || test_number==6) Initialize_Phi();
}
//#####################################################################
// Function Debug_Particle_Set_Attribute
//#####################################################################
template<class TV,class ATTR> void PhysBAM::
Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr)
{
    typedef typename TV::SCALAR T;
    GEOMETRY_PARTICLES<TV>* particles=(GEOMETRY_PARTICLES<TV>*)KANG<T>::Store_Debug_Particles();
    ARRAY_VIEW<ATTR>* attribute=particles->array_collection->template Get_Array<ATTR>(id);
    attribute->Last()=attr;
}
//#####################################################################
// Function Set_Boundary_Conditions_Callback
//#####################################################################
template<class T> void KANG<T>::
Set_Boundary_Conditions_Callback(ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& psi_D_value,
    ARRAY<T,FACE_INDEX<TV::dimension> >& psi_N_value) const
{
    if(test_number==4){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid,3,GRID<TV>::GHOST_REGION);it.Valid();it.Next()){
            TV x=it.Location();
            psi_D_value(it.index)=sqr(x.x)-sqr(x.y);}
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION);it.Valid();it.Next()){
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));
            TV x=it.Location();
            psi_N_value(it.Full_Index())=it.Axis()==1?2*x.x:-2*x.y;}}
    if(test_number==5){
        for(int s=1;s<=2;s++)
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,1,GRID<TV>::GHOST_REGION,s);it.Valid();it.Next()){
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));
                psi_N(it.Full_Index())=true;
                psi_N_value(it.Full_Index())=0;}}
    if(test_number==6){
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,1,GRID<TV>::GHOST_REGION,1);it.Valid();it.Next()){
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));
            psi_N(it.Full_Index())=true;
            psi_N_value(it.Full_Index())=-1;}
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid,1,GRID<TV>::GHOST_REGION,2);it.Valid();it.Next()){
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));
            psi_N(it.Full_Index())=true;
            psi_N_value(it.Full_Index())=1;}}
}
template class KANG<float>;
template void PhysBAM::Add_Debug_Particle<VECTOR<float,1> >(VECTOR<float,1> const&,VECTOR<float,3> const&);
template void PhysBAM::Add_Debug_Particle<VECTOR<float,2> >(VECTOR<float,2> const&,VECTOR<float,3> const&);
template void PhysBAM::Add_Debug_Particle<VECTOR<float,3> >(VECTOR<float,3> const&,VECTOR<float,3> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<float,1>,VECTOR<float,1> >(ATTRIBUTE_ID,VECTOR<float,1> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<float,2>,VECTOR<float,2> >(ATTRIBUTE_ID,VECTOR<float,2> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<float,3>,VECTOR<float,3> >(ATTRIBUTE_ID,VECTOR<float,3> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KANG<double>;
template void PhysBAM::Add_Debug_Particle<VECTOR<double,1> >(VECTOR<double,1> const&,VECTOR<double,3> const&);
template void PhysBAM::Add_Debug_Particle<VECTOR<double,2> >(VECTOR<double,2> const&,VECTOR<double,3> const&);
template void PhysBAM::Add_Debug_Particle<VECTOR<double,3> >(VECTOR<double,3> const&,VECTOR<double,3> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<double,1>,VECTOR<double,1> >(ATTRIBUTE_ID,VECTOR<double,1> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<double,2>,VECTOR<double,2> >(ATTRIBUTE_ID,VECTOR<double,2> const&);
template void PhysBAM::Debug_Particle_Set_Attribute<VECTOR<double,3>,VECTOR<double,3> >(ATTRIBUTE_ID,VECTOR<double,3> const&);
#endif
