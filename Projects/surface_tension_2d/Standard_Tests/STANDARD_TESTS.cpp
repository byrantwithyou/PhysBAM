#include <Grid_PDE/Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_CUT.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_CUT.h>
#include "STANDARD_TESTS.h"
using namespace PhysBAM;
namespace PhysBAM{
template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
template<class TV> void Add_Debug_Particle(const TV& X){Add_Debug_Particle(X,VECTOR<typename TV::SCALAR,3>(1,0,0));}
}
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4355) // 'this' : used in base member initializer list
#endif
template<class T> STANDARD_TESTS<T>::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :BASE(stream_type_input,parse_args,1,fluids_parameters.WATER),solids_tests(stream_type_input,data_directory,solid_body_collection),run_self_tests(false),print_poisson_matrix(false),print_index_map(false),
    print_matrix(false),print_each_matrix(false),use_full_ic(false),output_iterators(false),use_viscous_forces(false),max_dt(0),exact_dt(0),current_dt(0),implicit_solid(false),
    front_tracked_structure(0),rebuild_curve(0),deformable_collisions(0),fsi(0),number_surface_particles(0),rebuild_surface(false),free_particles(0),psi_D(0),
    circle_radius(0),circle_perturbation((T).05),oscillation_mode(2),use_massless_structure(false),coupled_particles(0),make_ellipse(false),m(1),s(1),kg(1),solid_refinement(7),
    solid_density(0),solid_width(0),analytic_solution(0),linear_force(0),rand(0),use_viscosity(false)
{
    LOG::cout<<std::setprecision(16);
    debug_particles.template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    debug_particles.Store_Velocity(true);
    Store_Debug_Particles(&debug_particles);

    fluids_parameters.use_preconditioner_for_slip_system=false;
    solids_parameters.implicit_solve_parameters.cg_iterations=3000;
    solids_fluids_parameters.use_leakproof_solve=false;
    fluids_parameters.surface_tension=(T).0728;

    parse_args.Add("-use_viscous_forces",&use_viscous_forces,"use_viscous_forces");
    parse_args.Add("-implicit_solid",&implicit_solid,"implicit_solid");
    parse_args.Add("-cg_iterations",&solids_parameters.implicit_solve_parameters.cg_iterations,"iterations","cg iterations");
    parse_args.Add("-viscosity",&fluids_parameters.viscosity,&use_viscosity,"viscosity","Viscosity");
    parse_args.Add("-test_system",&run_self_tests,"run_self_tests");
    parse_args.Add("-print_poisson_matrix",&print_poisson_matrix,"print_poisson_matrix");
    parse_args.Add("-print_index_map",&print_index_map,"print_index_map");
    parse_args.Add("-print_matrix",&print_matrix,"print_matrix");
    parse_args.Add("-print_each_matrix",&print_each_matrix,"print_each_matrix");
    parse_args.Add("-use_full_ic",&use_full_ic,"use_full_ic");
    parse_args.Add("-output_iterators",&output_iterators,"output_iterators");
    parse_args.Add_Not("-no_preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"do not use preconditioner");
    parse_args.Add("-preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"use preconditioner");
    parse_args.Add("-leakproof",&solids_fluids_parameters.use_leakproof_solve,"use leakproof solve");
    parse_args.Add("-max_dt",&max_dt,"dt","maximum dt");
    parse_args.Add("-dt",&exact_dt,"dt","use this constand dt");
    parse_args.Add("-build_surface",&rebuild_surface,"rebuild surface");
    parse_args.Add("-print_energy",&solid_body_collection.print_energy,"print energy statistics");
    parse_args.Add("-surface_tension",&fluids_parameters.surface_tension,"coeff","surface tension coefficient");
    parse_args.Add("-cut_cell",&fluids_parameters.second_order_cut_cell_method,"use second order cut cell");
    parse_args.Add("-oscillation_mode",&oscillation_mode,"modes","number of oscillation nodes");
    parse_args.Add("-solid_refinement",&solid_refinement,"level","solid refinement");
    parse_args.Add("-make_ellipse",&make_ellipse,"use ellipse for initial configuration");
    parse_args.Add("-epsilon",&circle_perturbation,"eps","circle perturbation");
    parse_args.Add("-m",&m,"unit","length unit");
    parse_args.Add("-s",&s,"unit","time unit");
    parse_args.Add("-kg",&kg,"unit","mass unit");
    parse_args.Add("-linear_force",&linear_force,"force","linear force");
    parse_args.Add("-rand",&rand,"rand","rand");
    parse_args.Parse();
    solids_tests.data_directory=data_directory;
    last_frame=100;
    frame_rate=24;

    frame_rate/=s;

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

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=200;
    solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    //solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
    fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;

    fluids_parameters.use_slip=true;
    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_preconditioner_for_slip_system=true;

    fluids_parameters.viscosity*=kg/s;
    if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=true;
    fluids_parameters.surface_tension*=kg*m/(s*s);

    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.store_particle_ids=true;
    // T default_removed_positive_particle_buoyancy_constant=fluids_parameters.removed_positive_particle_buoyancy_constant;
    fluids_parameters.removed_positive_particle_buoyancy_constant=0;
    //solid_body_collection.print_residuals=true;

    fluids_parameters.use_coupled_implicit_viscosity=use_viscous_forces;

    max_dt*=s;
    exact_dt*=s;
    solids_parameters.use_post_cg_constraints=true;

    fluids_parameters.gravity=TV();
    fluids_parameters.use_levelset_viscosity=false;
    fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.number_particles_per_cell=16;
    fluids_parameters.write_particles=true;
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;

    switch(test_number){
        case 1:
            fluids_parameters.grid->Initialize(TV_INT(resolution+1,resolution+1),RANGE<TV>(TV(0*m,0*m),TV((T).04*m,(T).04*m)));
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        case 2:
            fluids_parameters.grid->Initialize(TV_INT()+5,RANGE<TV>(TV()+10*m,TV()+11*m));
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
            fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            break;
        case 3:
            fluids_parameters.grid->Initialize(TV_INT(resolution+1,resolution+1),RANGE<TV>(TV(0*m,0*m),TV((T).04*m,(T).04*m)));
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        case 4:
        case 5:
            fluids_parameters.grid->Initialize(TV_INT(resolution+1,resolution+1),RANGE<TV>(TV(0*m,0*m),TV(1*m,1*m)));
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        case 6:
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            (*fluids_parameters.grid).Initialize(TV_INT(resolution+1,resolution+1),RANGE<TV>(TV(0,(T)0),TV(1,m)));
            if(!use_viscosity) fluids_parameters.viscosity=100;
            solid_density=150;
            solid_width=(T)1/3;
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}
        
    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);
    output_directory=LOG::sprintf("Standard_Tests/Test_%d_Resolution_%d",test_number,resolution);
}
#ifdef _MSC_VER
#pragma warning(pop)
#endif
template<class T> STANDARD_TESTS<T>::
~STANDARD_TESTS()
{
    delete deformable_collisions;
}
template<class T> void STANDARD_TESTS<T>::
After_Initialization() {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Initialize_Advection()
{
    fluids_parameters.particle_levelset_evolution->Levelset_Advection(1).
        Use_Semi_Lagrangian_Collidable_Advection(*fluids_parameters.collision_bodies_affecting_fluid,fluids_parameters.collidable_phi_replacement_value,
            fluids_parameters.incompressible->valid_mask);
    fluids_parameters.particle_levelset_evolution->Levelset(1).Set_Collision_Body_List(*fluids_parameters.collision_bodies_affecting_fluid);
    fluids_parameters.incompressible->Use_Semi_Lagrangian_Collidable_Advection(*fluids_parameters.collision_bodies_affecting_fluid);
    fluids_parameters.incompressible->collision_body_list=fluids_parameters.collision_bodies_affecting_fluid;

    fluids_parameters.particle_levelset_evolution->runge_kutta_order_particles=3;
    fluids_parameters.particle_levelset_evolution->runge_kutta_order_levelset=3;
    fluids_parameters.particle_levelset_evolution->Use_Hamilton_Jacobi_Weno_Advection();
    fluids_parameters.incompressible->advection=new ADVECTION_HAMILTON_JACOBI_ENO<TV,T>;
    fluids_parameters.particle_levelset_evolution->Use_Reinitialization();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Initialize_Phi()
{
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    if(test_number==1 || test_number==3){
        SPHERE<TV> object(TV((T).02*m,(T).02*m),(T).01*m);
        for(CELL_ITERATOR<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            TV X=it.Location();
            if(make_ellipse) X*=TV((T)1.1,(T).9);
            phi(it.index)=object.Signed_Distance(X);}}
    else if(test_number==2) phi.Fill(1);
    else if(test_number==4 || test_number==5){
        for(CELL_ITERATOR<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
            TV dX=it.Location()-TV((T)(.5*m),(T)(.5*m));
            T distance=dX.Magnitude();
            T angle=atan2(dX.y,dX.x);
            T radius=circle_radius+circle_perturbation*cos(oscillation_mode*angle);
            phi(it.index)=distance-radius;}}
    else if(test_number==6) phi.Fill(-1);
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Preprocess_Substep(const T dt,const T time)
{
    current_dt=dt;
    if(SURFACE_TENSION_FORCE<TV>* force=solid_body_collection.deformable_body_collection.template Find_Force<SURFACE_TENSION_FORCE<TV>*>()) force->dt=dt;
    if(LINEAR_POINT_ATTRACTION<TV>* force=solid_body_collection.deformable_body_collection.template Find_Force<LINEAR_POINT_ATTRACTION<TV>*>()) force->dt=dt;
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Preprocess_Frame(const int frame)
{
    if(fluids_parameters.use_slip){
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).run_self_tests=run_self_tests;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).output_iterators=output_iterators;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_matrix=print_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_each_matrix=print_each_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_poisson_matrix=print_poisson_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).use_viscous_forces=use_viscous_forces;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_index_map=print_index_map;}
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
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
template<class T> void STANDARD_TESTS<T>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    BASE::Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time)
{
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
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
        default: LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    if(test_number!=3) deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    SEGMENTED_CURVE_2D<T>* curve=rebuild_curve;
    if(!curve) curve=solid_body_collection.deformable_body_collection.template Find_Structure<SEGMENTED_CURVE_2D<T>*>();
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
template<class T> void STANDARD_TESTS<T>::
Kang_Circle(bool use_surface)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    fluids_parameters.gravity.y=-(T)0*m/(s*s);
    fluids_parameters.density=(T)1000*kg/(m*m);
    fluids_parameters.domain_walls[0][0]=true;
    fluids_parameters.domain_walls[0][1]=true;
    fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.cfl=(T).9;
//    solids_parameters.write_static_variables_every_frame=true;
    use_massless_structure=use_surface;

    solid_body_collection.Set_CFL_Number(10);

//    fluids_parameters.second_order_cut_cell_method=false;
    fluids_parameters.use_particle_levelset=true;

    if(!rebuild_surface && use_surface){
        SPHERE<TV> object(TV((T).02*m,(T).02*m),(T).01*m);
        front_tracked_structure=&solids_tests.Copy_And_Add_Structure(*TESSELLATION::Tessellate_Boundary(object,solid_refinement));
        solid_body_collection.deformable_body_collection.particles.Add_Elements(20*resolution);
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
template<class T> void STANDARD_TESTS<T>::
Oscillating_Circle(bool use_surface)
{
    (*fluids_parameters.grid).Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
    fluids_parameters.cfl=(T).9;
//    solids_parameters.write_static_variables_every_frame=true;
    use_massless_structure=use_surface;

    solid_body_collection.Set_CFL_Number(10);

    fluids_parameters.use_particle_levelset=true;

    circle_radius=(T)1/3*m;
    circle_perturbation=circle_radius/20;
    oscillation_mode*=circle_radius;
    fluids_parameters.surface_tension=(T)2/3*kg*m/(s*s);
    fluids_parameters.density=27*kg/(m*m);
    T omega=sqrt(oscillation_mode*(oscillation_mode*oscillation_mode-1)*fluids_parameters.surface_tension/(fluids_parameters.density*circle_radius*circle_radius*circle_radius));
    T period=(T)(2*pi/omega);
    LOG::cout<<"Analytic period of oscillation: "<<period<<std::endl;

    // Forces here
    if(linear_force){
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
        solid_body_collection.deformable_body_collection.particles.Add_Elements(20*resolution);
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
// Function Solid_Circle
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Solid_Circle()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SPHERE<TV> object(TV(),1);
    solids_tests.Copy_And_Add_Structure(*TESSELLATION::Tessellate_Boundary(object,solid_refinement));
    particles.mass+=(T)1/particles.mass.m;
    if(make_ellipse) for(int i=0;i<particles.X.m;i++) particles.X(i)*=TV((T).9,(T)1.1);
    fluids_parameters.viscosity=rand;
    if(rand){
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1234);
        ARRAY<TV> dX(particles.X.m);
        random.Fill_Uniform(dX,-rand*m,rand*m);
        particles.X+=dX;}
}
//#####################################################################
// Function Adjust_Phi_With_Objects
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
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
template<class T> void STANDARD_TESTS<T>::
Sync_Particle_To_Level_Set(int p)
{
    const LEVELSET<TV>& levelset=fluids_parameters.particle_levelset_evolution->Levelset(1);
    TV& X=front_tracked_structure->particles.X(p);
    X-=levelset.Extended_Phi(X)*levelset.Extended_Normal(X);
}
//#####################################################################
// Function Sync_Front_Tracked_Particles_To_Level_Set
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
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
template<class T> void STANDARD_TESTS<T>::
Divide_Segment(int e)
{
    VECTOR<int,2> seg=front_tracked_structure->mesh.elements(e);
    TV XL=front_tracked_structure->particles.X(seg.x),XR=front_tracked_structure->particles.X(seg.y);
    if((XL-XR).Magnitude_Squared()<sqr(fluids_parameters.grid->dX.Min())) return;

    LOG::cout<<"a e "<<front_tracked_structure->mesh.elements<<std::endl;
    LOG::cout<<"a p "<<particle_segments<<std::endl;
    int ne=front_tracked_structure->mesh.elements.m+1;
    int p=particle_segments.Append(VECTOR<int,2>(e,ne));
    LOG::cout<<p<<"  "<<front_tracked_structure->particles.Size()<<std::endl;
    PHYSBAM_ASSERT(p<=front_tracked_structure->particles.Size());
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
template<class T> void STANDARD_TESTS<T>::
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
template<class T> void STANDARD_TESTS<T>::
Swap_Segments(int e,int f)
{
    exchange(particle_segments(front_tracked_structure->mesh.elements(e).x).y,particle_segments(front_tracked_structure->mesh.elements(f).x).y);
    exchange(particle_segments(front_tracked_structure->mesh.elements(e).y).x,particle_segments(front_tracked_structure->mesh.elements(f).y).x);
    exchange(front_tracked_structure->mesh.elements(e),front_tracked_structure->mesh.elements(f));
}
//#####################################################################
// Function Remove_Particle
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
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
template<class T> T STANDARD_TESTS<T>::
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
template<class T> void STANDARD_TESTS<T>::
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
template<class T> void STANDARD_TESTS<T>::
Limit_Dt(T& dt,const T time)
{
    if(max_dt && dt>max_dt) dt=max_dt;
    if(exact_dt) dt=exact_dt;
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Limit_Solids_Dt(T& dt,const T time)
{
    if(max_dt && dt>max_dt) dt=max_dt;
    if(exact_dt) dt=exact_dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);
    if(debug_particles.Size() || frame==0){
        Create_Directory(LOG::sprintf("%s/%i",output_directory.c_str(),frame));
        Write_To_File(this->stream_type,LOG::sprintf("%s/%i/debug_particles",output_directory.c_str(),frame),debug_particles);}
}
//#####################################################################
// Function Initialize_Surface_Particles
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Initialize_Surface_Particles(int number)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    number_surface_particles=number;
    PHYSBAM_ASSERT(particles.Size()==0);
    particles.Add_Elements(number_surface_particles);
    rebuild_curve=SEGMENTED_CURVE_2D<T>::Create(particles);
    free_particles=FREE_PARTICLES<TV>::Create(particles);
    solid_body_collection.deformable_body_collection.Add_Structure(free_particles);
    free_particles->nodes=IDENTITY_ARRAY<>(number_surface_particles);
}
//#####################################################################
// Function Rebuild_Surface
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
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
    const ARRAY<T,TV_INT>& phi=fluids_parameters.particle_levelset_evolution->Levelset(1).phi;
    for(FACE_ITERATOR<TV> it(*fluids_parameters.grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face=it.Full_Index();
        TV_INT cell1=face.First_Cell_Index(),cell2=face.Second_Cell_Index();
        T phi1=phi(cell1),phi2=phi(cell2);
        if((phi1>0)==(phi2>0)) continue;
        T theta=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
        TV X0=fluids_parameters.grid->X(cell1),X1=fluids_parameters.grid->X(cell2);
        TV X=(1-theta)*X0+theta*X1;

        typename MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::ENTRY entry={face,phi1<=0?1:2,X};
        fluid_interpolation_entries.Append(entry);

        PHYSBAM_ASSERT(next<=number_surface_particles+1);
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));
        particles.X(next)=X;
        solid_interpolation_entries.Append(next);
        used_cells.Get_Or_Insert(cell2).Append(next);
        cell2(3-face.axis)++;
        used_cells.Get_Or_Insert(cell2).Append(next);
        next++;}

    for(typename HASHTABLE<TV_INT,ARRAY<int> >::ITERATOR it(used_cells);it.Valid();it.Next()){
        ARRAY<int>& array=it.Data();
        PHYSBAM_ASSERT(array.m==2);
        rebuild_curve->mesh.elements.Append(VECTOR<int,2>(array(1),array(2)));
        Add_Debug_Particle((rebuild_curve->particles.X(array(1))+rebuild_curve->particles.X(array(2)))/2,VECTOR<T,3>(1,0,0));}
    PHYSBAM_ASSERT(rebuild_curve->particles.X.Get_Array_Pointer()==particles.X.Get_Array_Pointer());
}
//#####################################################################
// Function Substitute_Coupling_Matrices
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Substitute_Coupling_Matrices(KRYLOV_SYSTEM_BASE<T>& coupled_system,T dt,T current_velocity_time,T current_position_time,bool velocity_update,bool leakproof_solve)
{
    if(use_massless_structure){
        SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>& system=dynamic_cast<SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>&>(coupled_system);
        fsi=new FLUID_TO_SOLID_INTERPOLATION_CUT<TV>(system.index_map,*front_tracked_structure,fluids_parameters.density);
        system.fluid_to_solid_interpolation=fsi;
        MATRIX_FLUID_GRADIENT_CUT<TV>* gradient=new MATRIX_FLUID_GRADIENT_CUT<TV>(system.index_map);
        system.fluid_gradient=gradient;
        fsi->fluid_mass=&system.fluid_mass;
        fsi->gradient=gradient;}

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
template<class T> void STANDARD_TESTS<T>::
Advance_One_Time_Step_Begin_Callback(const T dt,const T time)
{
    if(rebuild_curve) Rebuild_Surface();
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Update_Time_Varying_Material_Properties(const T time)
{
    if(!time) return;
    if(rebuild_curve) Rebuild_Surface();
}
//#####################################################################
// Function Store_Debug_Particles
//#####################################################################
template<class T> GEOMETRY_PARTICLES<VECTOR<T,2> >*  STANDARD_TESTS<T>::
Store_Debug_Particles(GEOMETRY_PARTICLES<TV>* particle)
{
    static GEOMETRY_PARTICLES<TV>* stored_particles=0;
    GEOMETRY_PARTICLES<TV>* tmp=stored_particles;
    if(particle) stored_particles=particle;
    return tmp;
}
//#####################################################################
// Function FSI_Analytic_Test
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
FSI_Analytic_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    fluids_parameters.collision_bodies_affecting_fluid->use_collision_face_neighbors=true;
    TV solid_gravity=TV(0,-(T)9.8*m/(s*s));
    fluids_parameters.surface_tension=0;
    debug_particles.template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);

    fluids_parameters.gravity.y=-(T)9.8*m;
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
    rigid_body.Frame()=FRAME<TV>(m*TV((T).5,.5));
    rigid_body.Set_Coefficient_Of_Restitution((T)0);
    rigid_body.Set_Mass(solid_density*solid_width);

    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity*m));
    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body);

    T solid_mass=solid_body_collection.rigid_body_collection.rigid_body_particles.mass(3);
    T rho=fluids_parameters.density;
    TV size=fluids_parameters.grid->domain.Edge_Lengths();
    size.x=(size.x-solid_width)/2;
    analytic_solution=-(solid_mass*-solid_gravity.y+rho*size.x*size.y*fluids_parameters.gravity)*size.x/(2*size.y*fluids_parameters.viscosity);
    LOG::cout<<"analytic_solution "<<analytic_solution<<std::endl;

    Create_Directory(LOG::sprintf("%s/%i",output_directory.c_str(),0));
    Write_To_File(this->stream_type,LOG::sprintf("%s/%i/debug_particles",output_directory.c_str(),0),debug_particles);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Postprocess_Frame(const int frame)
{
    if(test_number==6){
        T v=fluid_collection.incompressible_fluid_collection.face_velocities(FACE_INDEX<2>(2,fluids_parameters.grid->counts/2));
        LOG::cout<<"middle-velocity "<<v<<"   error from analytic solution "<<(v/analytic_solution-1)<<std::endl;}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T> void STANDARD_TESTS<T>::
Postprocess_Substep(const T dt,const T time)
{
}
template class STANDARD_TESTS<float>;
template class STANDARD_TESTS<double>;
