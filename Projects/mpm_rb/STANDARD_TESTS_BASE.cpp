//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Constitutive_Models/ST_VENANT_KIRCHHOFF_HENCKY_STRAIN.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
#include <Solids/Collisions/PENALTY_FORCE_COLLECTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_DRUCKER_PRAGER.h>
#include <Hybrid_Methods/Forces/MPM_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_GRAVITY.h>
#include <Hybrid_Methods/Forces/MPM_PLASTIC_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_CLAMP.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include "STANDARD_TESTS_BASE.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :MPM_EXAMPLE_RB<TV>(stream_type_input),test_number(0),resolution(32),user_resolution(false),stored_last_frame(0),
    user_last_frame(false),order(2),seed(1234),particles_per_cell(1<<TV::m),regular_seeding(false),
    no_regular_seeding(false),scale_mass(1),scale_E(1),scale_speed(1),
    penalty_collisions_stiffness((T)1e4),penalty_collisions_separation((T)1e-4),penalty_collisions_length(1),
    penalty_damping_stiffness(0),use_penalty_collisions(false),use_plasticity(true),
    use_theta_c(false),use_theta_s(false),use_hardening_factor(false),use_max_hardening(false),
    theta_c(0),theta_s(0),hardening_factor(0),max_hardening(0),use_implicit_plasticity(false),no_implicit_plasticity(false),
    hardening_mast_case(0),use_hardening_mast_case(false),override_output_directory(false),
    m(1),s(1),kg(1),forced_collision_type(-1),friction(0),friction_is_set(false),sigma_Y(0),use_cohesion(false),write_output_files(0),read_output_files(0),
    dump_collision_objects(false),tests(STREAM_TYPE(0.f),data_directory,solid_body_collection)
{
    T framerate=0;
    bool use_quasi_exp_F_update=false;
    bool no_affine=false;
    bool use_separate=false,use_slip=false,use_stick=false;
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Add("-restart",&restart,"frame","restart frame");
    parse_args.Add("-resolution",&resolution,&user_resolution,"resolution","grid resolution");
    parse_args.Add("-substeps",&write_substeps_level,"level","output-substep level");
    parse_args.Add("-substeps_delay",&substeps_delay_frame,"frame","delay substeps until after this frame");
    parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
    parse_args.Add("-test_diff",&test_diff,"test analytic derivatives");
    parse_args.Add("-test_system",&test_system,"test Krylov system properties");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Add("-o",&output_directory,&override_output_directory,"dir","Output directory");
    parse_args.Add("-mass_contour",&mass_contour,"contour","Draw mass contour as a scale to particle average mass");
    parse_args.Add("-framerate",&framerate,"rate","Number of frames per second");
    parse_args.Add("-frame_dt",&frame_dt,"rate","Number of frames per second");
    parse_args.Add("-min_dt",&min_dt,"dt","Minimum time step size");
    parse_args.Add("-max_dt",&max_dt,"dt","Maximum time step size");
    parse_args.Add("-order",&order,"order","Interpolation basis order");
    parse_args.Add_Not("-no_affine",&use_affine,"Use affine PIC");
    parse_args.Add("-affine",&use_affine,"Use affine PIC");
    parse_args.Add("-midpoint",&use_midpoint,"Use midpoint rule");
    parse_args.Add("-symplectic_euler",&use_symplectic_euler,"Use forward euler for grid update");
    parse_args.Add("-print_stats",&print_stats,"Print momentum/energy stats");
    parse_args.Add("-only_write_particles",&only_write_particles,"Only write particle data (ignore grid data, restart data etc)");
    parse_args.Add("-flip",&flip,&no_affine,"frac","Flip ratio");
    parse_args.Add("-cfl",&cfl,"cfl","CFL number");
    parse_args.Add("-cfl_F",&cfl_F,"cfl","CFL number for F update");
    parse_args.Add("-cfl_c",&cfl_sound,"cfl","CFL number for sound speed");
    parse_args.Add("-newton_tolerance",&newton_tolerance,"tol","Newton tolerance");
    parse_args.Add("-newton_iterations",&newton_iterations,"iter","Newton iterations");
    parse_args.Add("-angle_tol",&angle_tol,"tol","Newton iterations");
    parse_args.Add("-solver_tolerance",&solver_tolerance,"tol","Solver tolerance");
    parse_args.Add("-solver_iterations",&solver_iterations,"iter","Solver iterations");
    parse_args.Add("-seed",&seed,"seed","Random number seed");
    parse_args.Add("-particles_per_cell",&particles_per_cell,"num","Number of particles per cell");
    parse_args.Add("-scale_mass",&scale_mass,"scale","Scale mass of particles");
    parse_args.Add("-scale_E",&scale_E,"scale","Scale stiffness");
    parse_args.Add("-scale_speed",&scale_speed,"scale","Scale initial speed of simulated object");
    parse_args.Add("-penalty_stiffness",&penalty_collisions_stiffness,"tol","penalty collisions stiffness");
    parse_args.Add("-penalty_separation",&penalty_collisions_separation,"tol","penalty collisions separation");
    parse_args.Add("-penalty_length",&penalty_collisions_length,"tol","penalty collisions length scale");
    parse_args.Add("-penalty_damping",&penalty_damping_stiffness,"tol","penalty damping stiffness");
    parse_args.Add("-regular_seeding",&regular_seeding,"use regular particle seeding");
    parse_args.Add_Not("-no_regular_seeding",&no_regular_seeding,"use regular particle seeding");
    parse_args.Add("-lag_Dp",&lag_Dp,"use early gradient transfer for Cp");
    parse_args.Add("-use_exp_F",&use_quasi_exp_F_update,"Use an approximation of the F update that prevents inversion");
    parse_args.Add("-use_plasticity",&use_plasticity,"Use plasticity in the F update");
    parse_args.Add("-theta_c",&theta_c,&use_theta_c,"theta_c","Critical compression coefficient for plasticity");
    parse_args.Add("-theta_s",&theta_s,&use_theta_s,"theta_s","Critical stretch coefficient for plasticity");
    parse_args.Add("-hardening",&hardening_factor,&use_hardening_factor,"hardening factor","Hardening factor for plasticity");
    parse_args.Add("-max_hardening",&max_hardening,&use_max_hardening,"max hardening coefficient","Maximum hardening coefficient for plasticity");
    parse_args.Add("-use_penalty_collisions",&use_penalty_collisions,"Use penalty collisions objects");
    parse_args.Add("-use_implicit_plasticity",&use_implicit_plasticity,"Use implicit plasticity");
    parse_args.Add("-no_implicit_plasticity",&no_implicit_plasticity,"Disable implicit plasticity");
    parse_args.Add("-mast_case",&hardening_mast_case,&use_hardening_mast_case,"mast_case","The case number from the Mast thesis for hardening parameters");
    parse_args.Add("-m",&m,"scale","meter scale");
    parse_args.Add("-s",&s,"scale","second scale");
    parse_args.Add("-kg",&kg,"scale","kilogram scale");
    parse_args.Add("-slip",&use_slip,"force slip collisions");
    parse_args.Add("-stick",&use_stick,"force stick collisions");
    parse_args.Add("-separate",&use_separate,"force separating collisions");
    parse_args.Add("-friction",&friction,&friction_is_set,"friction","Coefficient of friction");
    parse_args.Add("-cohesion",&sigma_Y,&use_cohesion,"cohesion","sigma_Y in Drucker Prager");
    parse_args.Add("-T",&extra_T,"float","extra float argument");
    parse_args.Add("-I",&extra_int,"int","extra int argument");
    parse_args.Add("-dump_collisions",&dump_collision_objects,"dump out collision objects");
    parse_args.Add("-test_output_prefix",&test_output_prefix,&use_test_output,"","prefix to use for test output");
    parse_args.Add("-strong_cfl",&use_strong_cfl,"limit dt based on final change in X and F");
    parse_args.Add("-sound_cfl",&use_sound_speed_cfl,"limit dt based on final change in X and F");
    parse_args.Add("-reflection_bc",&reflection_bc_flags,"flags","Flags indicating which walls should be reflection BC");
    parse_args.Add("-coll_pair",&pairwise_collisions,"use pairwise collisions");
    parse_args.Add("-coll_proj",&projected_collisions,"use pairwise collisions");
    parse_args.Add("-rd_stiffness",&rd_penalty_stiffness,"stiffness","rigid-deformable penalty force stiffness");
    parse_args.Add("-rd_friction",&rd_penalty_friction,"friction","rigid-deformable penalty force friction");
    parse_args.Add("-rd",&use_rd,"enable rigid-deformable penalty force friction");
    parse_args.Add("-rr",&use_rr,"enable rigid-rigid penalty force friction");
    parse_args.Add("-di",&use_di,"enable deformable-object penalty force friction");
    parse_args.Add("-grad_ls",&use_gradient_magnitude_objective,"do line searches on norm of gradient");
    parse_args.Add("-debug_newton",&debug_newton,"Enable diagnostics in Newton's method");
    parse_args.Add("-rd_k",&rd_k,&use_rd_k,"stiffness","override stiffness for rigid-deformable penalty force friction");
    parse_args.Add("-rr_k",&rr_k,&use_rr_k,"stiffness","override stiffness for rigid-rigid penalty force friction");
    parse_args.Add("-di_k",&di_k,&use_di_k,"stiffness","override stiffness for deformable-object penalty force friction");
    parse_args.Add("-rd_mu",&rd_mu,&use_rd_mu,"friction","override friction for rigid-deformable penalty force friction");
    parse_args.Add("-rr_mu",&rr_mu,&use_rr_mu,"friction","override friction for rigid-rigid penalty force friction");
    parse_args.Add("-di_mu",&di_mu,&use_di_mu,"friction","override friction for deformable-object penalty force friction");
    parse_args.Add("-bisection",&this->use_bisection,"use bisection relaxation");
    
    parse_args.Parse(true);
    PHYSBAM_ASSERT((int)use_slip+(int)use_stick+(int)use_separate<=1);
    if(use_slip) forced_collision_type=COLLISION_TYPE::slip;
    if(use_stick) forced_collision_type=COLLISION_TYPE::stick;
    if(use_separate) forced_collision_type=COLLISION_TYPE::separate;

    unit_p=kg*pow<2-TV::m>(m)/(s*s);
    unit_rho=kg*pow<-TV::m>(m);
    unit_mu=kg*pow<2-TV::m>(m)/s;
    mass_contour*=m/kg;
    min_dt*=s;
    max_dt*=s;
    penalty_collisions_separation*=m;
    penalty_collisions_length*=m;
    newton_tolerance*=kg*m/(s*s); // TODO: check me.

    if(no_affine) use_affine=false;

    if(framerate) frame_dt=1/framerate;
    frame_dt*=s;

    quad_F_coeff=(T)0;
    if(use_quasi_exp_F_update) quad_F_coeff=(T).5;

#ifdef USE_OPENMP
    omp_set_num_threads(threads);
#pragma omp parallel
#pragma omp single
    {
        if(omp_get_num_threads()!=threads) PHYSBAM_FATAL_ERROR();
        LOG::cout<<"Running on "<<threads<<" threads"<<std::endl;
    }
#else
    PHYSBAM_ASSERT(threads==1);
#endif

    gather_scatter.threads=threads;
    stored_last_frame=last_frame;
    random.Set_Seed(seed);

    particles.Store_Fp(use_plasticity);
    particles.Store_B(use_affine);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
~STANDARD_TESTS_BASE()
{
    if(destroy) destroy();
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Seed_Particles_Poisson(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell,const char* name)
{
    POISSON_DISK<TV> poisson_disk(1);
    ARRAY<TV> X;
    poisson_disk.Set_Distance_By_Volume(grid.dX.Product()/particles_per_cell);
    poisson_disk.Sample(random,object,X);

    T volume=grid.dX.Product()/particles_per_cell;
    T mass=density*volume;
    LOG::printf("MPM OBJECT %s BEGIN %d\n",name,particles.number);
    for(int i=0;i<X.m;i++)
        Add_Particle(X(i),V,dV,mass,volume);
    LOG::printf("MPM OBJECT %s END %d\n",name,particles.number);
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Seed_Particles_Uniform(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,const GRID<TV>& seed_grid,const char* name)
{
    T volume=seed_grid.dX.Product();
    T mass=density*volume;
    LOG::printf("MPM OBJECT %s BEGIN %d\n",name,particles.number);
    for(CELL_ITERATOR<TV> it(seed_grid);it.Valid();it.Next()){
        TV X=it.Location();
        if(object.Lazy_Inside(X))
            Add_Particle(X,V,dV,mass,volume);}
    LOG::printf("MPM OBJECT %s END %d\n",name,particles.number);
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Seed_Particles(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell,const char* name)
{
    if(!regular_seeding) return Seed_Particles_Poisson(object,V,dV,density,particles_per_cell,name);

    object.Update_Box();
    RANGE<TV_INT> range=grid.Clamp_To_Cell(object.Box(),3).Intersect(grid.Cell_Indices());
    TV LB=grid.Node(range.min_corner);
    TV UB=grid.Node(range.max_corner);
    T scale=pow<1,TV::m>((T)particles_per_cell);
    GRID<TV> seed_grid(range.Edge_Lengths()*scale,RANGE<TV>(LB,UB),true);
    Seed_Particles_Uniform(object,V,dV,density,seed_grid,name);
}
//#####################################################################
// Function Perturb
//#####################################################################
template<class TV> auto STANDARD_TESTS_BASE<TV>::
Perturb(T a) -> T
{
    return random.Get_Uniform_Number(1-a,1+a);
}
//#####################################################################
// Function Uniform
//#####################################################################
template<class TV> auto STANDARD_TESTS_BASE<TV>::
Uniform(T a,T b) -> T
{
    return random.Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Add_Particle
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Particle(const TV& X,std::function<TV(const TV&)> V,std::function<MATRIX<T,TV::m>(const TV&)> dV,
    const T mass,const T volume)
{
    int p=particles.Add_Element();
    particles.valid(p)=true;
    particles.X(p)=X;
    if(V) particles.V(p)=V(X);
    particles.F(p)=MATRIX<T,TV::m>()+1;
    if(particles.store_Fp) particles.Fp(p).Set_Identity_Matrix();
    if(particles.store_B && dV){
        if(lag_Dp) particles.B(p)=dV(X);
        else
            particles.B(p)=dV(X)*weights->Dp_Inverse(X).Inverse();}
    if(particles.store_S) particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
    particles.mass(p)=mass;
    particles.volume(p)=volume;
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
    (*color_attribute)(p)=VECTOR<T,3>(1,1,1);
}
//#####################################################################
// Function Add_Lambda_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Lambda_Particles(ARRAY<int>* affected_particles,T E,T nu,T density,bool no_mu,T porosity,T saturation_level)
{
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
    ARRAY<int> lambda_particles(affected_particles->m);
    T volume_lambda=particles.volume(0)*porosity*saturation_level;
    T mass_lambda=density*volume_lambda;
    T lambda=E*nu/((1+nu)*(1-2*nu));
    for(int k=0;k<affected_particles->m;k++){
        int p=particles.Add_Element();
        lambda_particles(k)=p;
        int i=(*affected_particles)(k);
        particles.valid(p)=true;
        particles.X(p)=particles.X(i);
        particles.V(p)=particles.V(i);
        particles.F(p)=particles.F(i);
        if(particles.store_Fp) particles.Fp(p)=particles.Fp(i); 
        if(particles.store_B) particles.B(p)=particles.B(i);
        if(particles.store_S) particles.S(p)=particles.S(i);
        particles.mass(p)=mass_lambda;
        particles.volume(p)=volume_lambda;
        particles.mu(p)=(T)0;
        particles.mu0(p)=(T)0;
        particles.lambda(p)=lambda;
        particles.lambda0(p)=lambda;
        (*color_attribute)(p)=VECTOR<T,3>(1,0,0);}
    Add_Fixed_Corotated(E,nu,&lambda_particles,no_mu);
}
//#####################################################################
// Function Add_Gravity
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Gravity(TV g,ARRAY<int>* affected_particles)
{
    return Add_Force(*new MPM_GRAVITY<TV>(force_helper,g,gather_scatter,affected_particles));
}
//#####################################################################
// Function Add_Gravity2
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Gravity2(TV g,ARRAY<int>* affected_particles)
{
    if(affected_particles)
        return Add_Force(*new DEFORMABLE_GRAVITY<TV>(particles,affected_particles,g));
    return Add_Force(*new DEFORMABLE_GRAVITY<TV>(particles,true,g));
}
//#####################################################################
// Function Add_Fixed_Corotated
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Fixed_Corotated(T E,T nu,ARRAY<int>* affected_particles,bool no_mu)
{
    COROTATED_FIXED<T,TV::m>* coro=new COROTATED_FIXED<T,TV::m>(E,nu);
    if(no_mu){nu=0;coro->Zero_Out_Mu();}
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model=*coro;
    MPM_FINITE_ELEMENTS<TV>& fe=*new MPM_FINITE_ELEMENTS<TV>(force_helper,constitutive_model,gather_scatter,affected_particles);
    Set_Lame_On_Particles(E,nu,affected_particles);
    return Add_Force(fe);
}
//#####################################################################
// Function Add_Neo_Hookean
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Neo_Hookean(T E,T nu,ARRAY<int>* affected_particles)
{
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model=*new NEO_HOOKEAN<T,TV::m>(E,nu);
    MPM_FINITE_ELEMENTS<TV>& fe=*new MPM_FINITE_ELEMENTS<TV>(force_helper,constitutive_model,gather_scatter,affected_particles);
    Set_Lame_On_Particles(E,nu,affected_particles);
    return Add_Force(fe);
}
//#####################################################################
// Function Add_St_Venant_Kirchhoff_Hencky_Strain
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_St_Venant_Kirchhoff_Hencky_Strain(T E,T nu,ARRAY<int>* affected_particles,bool no_mu)
{
    ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,TV::m>* hencky=new ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,TV::m>(E,nu);
    if(no_mu){nu=0;hencky->Zero_Out_Mu();}
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model=*hencky;
    MPM_FINITE_ELEMENTS<TV>& fe=*new MPM_FINITE_ELEMENTS<TV>(force_helper,constitutive_model,gather_scatter,affected_particles);
    Set_Lame_On_Particles(E,nu,affected_particles);
    return Add_Force(fe);
}
//#####################################################################
// Function Add_Drucker_Prager
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Drucker_Prager(T E,T nu,T a0,T a1,T a3,T a4,ARRAY<int>* affected_particles,bool no_mu,T sigma_Y)
{
    ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,TV::m>* hencky=new ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,TV::m>(E,nu);
    if(no_mu){nu=0;hencky->Zero_Out_Mu();}
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model=*hencky;
    MPM_DRUCKER_PRAGER<TV>& plasticity=*new MPM_DRUCKER_PRAGER<TV>(particles,0,a0,a1,a3,a4);
    plasticity.use_implicit=use_implicit_plasticity;
    PARTICLE_GRID_FORCES<TV>* fe=0;
    if(use_implicit_plasticity){
        fe=new MPM_PLASTIC_FINITE_ELEMENTS<TV>(force_helper,constitutive_model,gather_scatter,affected_particles,plasticity);
        this->asymmetric_system=true;}
    else{
        MPM_FINITE_ELEMENTS<TV>* mfe=new MPM_FINITE_ELEMENTS<TV>(force_helper,constitutive_model,gather_scatter,affected_particles);
        plasticity.gather_scatter=&mfe->gather_scatter;
        fe=mfe;}
    plasticity.Initialize_Particles(affected_particles,sigma_Y);
    plasticity_models.Append(&plasticity);
    Set_Lame_On_Particles(E,nu,affected_particles);
    return Add_Force(*fe);
}
//#####################################################################
// Function Add_Drucker_Prager
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Drucker_Prager(T E,T nu,T phi_F,ARRAY<int>* affected_particles,bool no_mu,T sigma_Y)
{
    return Add_Drucker_Prager(E,nu,phi_F,0,0,0,affected_particles,no_mu,sigma_Y);
}
//#####################################################################
// Function Add_Drucker_Prager
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Drucker_Prager_Case(T E,T nu,int case_num,ARRAY<int>* affected_particles,bool no_mu)
{
    ARRAY<ARRAY<T>> mast_constants({
        {35,0,0.2,10},
        {35,4,0.29,10},
        {35,9,0.3,10},
        {35,13,0.27,10},
        {35,0,0.2,6.57},
        {35,0,0.2,3.33},
        {35,0,0.2,0},
        {38.33,0,0.2,13.33},
        {41.67,0,0.2,16.67},
        {45,0,0.2,20}});
    ARRAY<T>& a=mast_constants(case_num);
    return Add_Drucker_Prager(E,nu,a(0),a(1),a(2),a(3),affected_particles,no_mu);
}
//#####################################################################
// Function Add_Clamped_Plasticity
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Clamped_Plasticity(ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& icm,T theta_c,T theta_s,T max_hardening,
    T hardening_factor,ARRAY<int>* affected_particles)
{
    MPM_PLASTICITY_CLAMP<TV>& plasticity=*new MPM_PLASTICITY_CLAMP<TV>(particles,0,theta_c,theta_s,
        max_hardening,hardening_factor);
    plasticity.use_implicit=use_implicit_plasticity;
    PARTICLE_GRID_FORCES<TV>* fe=0;
    if(use_implicit_plasticity){
        fe=new MPM_PLASTIC_FINITE_ELEMENTS<TV>(force_helper,icm,gather_scatter,affected_particles,plasticity);
        this->asymmetric_system=true;}
    else{
        MPM_FINITE_ELEMENTS<TV>* mfe=new MPM_FINITE_ELEMENTS<TV>(force_helper,icm,gather_scatter,affected_particles);
        plasticity.gather_scatter=&mfe->gather_scatter;
        fe=mfe;}
    plasticity.Initialize_Particles(affected_particles);
    plasticity_models.Append(&plasticity);
    return Add_Force(*fe);
}
//#####################################################################
// Function Add_Walls
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Walls(int flags,COLLISION_TYPE type,T friction,T inset,bool penalty) // -x +x -y +y [ -z +z ], as bit flags
{
    RANGE<TV> range=grid.domain.Thickened(grid.dX*(ghost*2+1));
    for(int a=0;a<TV::m;a++)
        for(int s=0;s<2;s++)
            if(flags&(1<<(a*2+s))){
                RANGE<TV> wall=range;
                if(s) wall.max_corner(a)=grid.domain.min_corner(a)+inset;
                else wall.min_corner(a)=grid.domain.max_corner(a)-inset;
                if(penalty) Add_Penalty_Collision_Object(wall);
                else Add_Collision_Object(wall,type,friction);}
}
//#####################################################################
// Function Seed_Lagrangian_Particles
//#####################################################################
template<class TV> template<class T_STRUCTURE> T_STRUCTURE& STANDARD_TESTS_BASE<TV>::
Seed_Lagrangian_Particles(T_STRUCTURE& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,bool use_constant_mass,bool destroy_after)
{
    int old_particles_number=particles.number;
    T_STRUCTURE& new_object=tests.Copy_And_Add_Structure(object,0,destroy_after);
    tests.Set_Mass_Of_Particles(new_object,density,use_constant_mass);
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
    for(int p=old_particles_number;p<particles.number;p++){
        TV X=particles.X(p);
        particles.valid(p)=true;
        if(V) particles.V(p)=V(X);
        if(particles.store_B && dV){
            if(lag_Dp) particles.B(p)=dV(X);
            else
                particles.B(p)=dV(X)*weights->Dp_Inverse(X).Inverse();}
        particles.F(p)=MATRIX<T,TV::m>()+1;
        if(particles.store_Fp) particles.Fp(p).Set_Identity_Matrix();
        if(particles.store_S) particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
        particles.volume(p)=particles.mass(p)/density;
        (*color_attribute)(p)=VECTOR<T,3>(1,1,1);}
    return new_object;
}
//#####################################################################
// Function Add_Fixed_Corotated
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Fixed_Corotated(T_VOLUME& object,T E,T nu)
{
    return Add_Force(*Create_Finite_Volume(object,new COROTATED_FIXED<T,TV::m>(E,nu,0)));
}
//#####################################################################
// Function Add_Neo_Hookean
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Neo_Hookean(T_VOLUME& object,T E,T nu)
{
    return Add_Force(*Create_Finite_Volume(object,new NEO_HOOKEAN<T,TV::m>(E,nu,0,(T).25)));
}
//#####################################################################
// Function Add_Penalty_Collision_Object
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Penalty_Collision_Object(IMPLICIT_OBJECT<TV>* io,const T coefficient_of_friction)
{
    IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>* pf=new IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>(particles,io,penalty_collisions_stiffness*kg/(m*s*s),penalty_collisions_separation,penalty_collisions_length);
    pf->coefficient_of_friction=coefficient_of_friction;
    this->Add_Force(*pf);
}
//#####################################################################
// Function Set_Lame_On_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_Lame_On_Particles(T E,T nu,ARRAY<int>* affected_particles)
{
    particles.Store_Lame(true);
    T mu=E/(2*(1+nu));
    T lambda=E*nu/((1+nu)*(1-2*nu));
    if(affected_particles)
#pragma omp parallel for
        for(int k=0;k<affected_particles->m;k++){
            int p=(*affected_particles)(k);
            particles.mu(p)=mu;
            particles.mu0(p)=mu;
            particles.lambda(p)=lambda;
            particles.lambda0(p)=lambda;}
    else{
        particles.mu.Fill(mu);
        particles.mu0.Fill(mu);
        particles.lambda.Fill(lambda);
        particles.lambda0.Fill(lambda);}
}
//#####################################################################
// Function Set_Grid
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_Grid(const RANGE<TV>& domain,TV_INT resolution_scale,int default_resolution)
{
    if(!user_resolution) resolution=default_resolution;
    grid.Initialize(resolution_scale*resolution,domain,true);
    Set_Weights(order);
}
//#####################################################################
// Function Set_Grid
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_Grid(const RANGE<TV>& domain,TV_INT resolution_scale,TV_INT resolution_padding,
    int resolution_multiple,int default_resolution)
{
    if(!user_resolution) resolution=default_resolution;
    int scaled_resolution=(resolution+resolution_multiple-1)/resolution_multiple;
    grid.Initialize(resolution_scale*scaled_resolution+resolution_padding,domain,true);
    Set_Weights(order);
}
//#####################################################################
// Function Add_Collision_Object
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Collision_Object(IMPLICIT_OBJECT<TV>* io)
{
    if(!use_di) return;
    if(use_di && !pfd){
        pfd=new PENALTY_FORCE_COLLECTION<TV>(solid_body_collection,simulated_particles,this->move_rb_diff);
        pfd->Init(use_di,false,use_rd,use_rr);}
    pfd->di_penalty->ios.Append(io);
}
//#####################################################################
// Function Add_Collision_Object
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,std::function<FRAME<TV>(T)> func_frame,std::function<TWIST<TV>(T)> func_twist)
{
    collision_objects.Append(new MPM_COLLISION_IMPLICIT_OBJECT<TV>(io,type,friction,func_frame,func_twist));
}
//#####################################################################
// Function Seed_Particles_Surface
//#####################################################################
template<class TV> template<class T_STRUCTURE> void STANDARD_TESTS_BASE<TV>::
Seed_Particles_Volume(T_STRUCTURE& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density)
{
    TETRAHEDRALIZED_VOLUME<T>& tv=object;
    ARRAY<int> indices;
    tests.Copy_And_Add_Structure(tv.Get_Boundary_Object(),&indices,false);
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
    T total_volume=tv.Total_Volume();
    T total_mass=total_volume*density;
    for(int i=0;i<indices.m;i++){
        int p=indices(i);
        TV X=particles.X(p);
        particles.valid(p)=true;
        particles.F(p)=MATRIX<T,TV::m>()+1;
        if(particles.store_Fp) particles.Fp(p)=MATRIX<T,TV::m>()+1;
        if(particles.store_S) particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
        if(V) particles.V(p)=V(X);
        (*color_attribute)(p)=TV(1,1,1);
        if(particles.store_B && dV){
            if(lag_Dp) particles.B(p)=dV(X);
            else
                particles.B(p)=dV(X)*weights->Dp_Inverse(X).Inverse();}
        particles.mass(p)=total_mass/indices.m;
        particles.volume(p)=total_volume/indices.m;}
}
//#####################################################################
// Function Seed_Particles_Surface
//#####################################################################
template<class TV> template<class T_STRUCTURE> void STANDARD_TESTS_BASE<TV>::
Seed_Particles_Surface(const T_STRUCTURE& object,IMPLICIT_OBJECT<TV>& io,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
{
    const TRIANGULATED_SURFACE<T>& ts=object;
    ARRAY<int> indices;
    tests.Copy_And_Add_Structure(ts,&indices,false);

    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
    T volume=grid.dX.Product()/particles_per_cell;
    T mass=density*volume;
    for(int i=0;i<indices.m;i++){
        int p=indices(i);
        TV X=particles.X(p);
        particles.valid(p)=true;
        particles.F(p)=MATRIX<T,TV::m>()+1;
        if(particles.store_Fp) particles.Fp(p)=MATRIX<T,TV::m>()+1;
        if(particles.store_S) particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
        if(V) particles.V(p)=V(X);
        (*color_attribute)(p)=TV(1,1,1);
        if(particles.store_B && dV){
            if(lag_Dp) particles.B(p)=dV(X);
            else
                particles.B(p)=dV(X)*weights->Dp_Inverse(X).Inverse();}
        particles.mass(p)=mass;
        particles.volume(p)=volume;}

    ARRAY<TV> X(ts.particles.X);
    POISSON_DISK<TV> poisson_disk(1);
    poisson_disk.Set_Distance_By_Volume(volume);
    poisson_disk.Sample(random,io,X);
    for(int p=ts.particles.X.m;p<X.m;p++)
        Add_Particle(X(p),V,dV,mass,volume);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Write_Output_Files(const int frame)
{
    Write_To_File(stream_type,LOG::sprintf("%s/%d/random_number",output_directory.c_str(),frame),random);
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Read_Output_Files(const int frame)
{
    Read_From_File(LOG::sprintf("%s/%d/random_number",output_directory.c_str(),frame),random);
    BASE::Read_Output_Files(frame);
}
template class STANDARD_TESTS_BASE<VECTOR<float,2> >;
template class STANDARD_TESTS_BASE<VECTOR<float,3> >;
template class STANDARD_TESTS_BASE<VECTOR<double,2> >;
template class STANDARD_TESTS_BASE<VECTOR<double,3> >;
template TRIANGULATED_AREA<double>& STANDARD_TESTS_BASE<VECTOR<double,2> >::Seed_Lagrangian_Particles<TRIANGULATED_AREA<double> >(TRIANGULATED_AREA<double>&,std::function<VECTOR<double,2> (VECTOR<double,2> const&)>,std::function<MATRIX<double,2,2> (VECTOR<double,2> const&)>,double,bool,bool);
template TRIANGULATED_AREA<float>& STANDARD_TESTS_BASE<VECTOR<float,2> >::Seed_Lagrangian_Particles<TRIANGULATED_AREA<float> >(TRIANGULATED_AREA<float>&,std::function<VECTOR<float,2> (VECTOR<float,2> const&)>,std::function<MATRIX<float,2,2> (VECTOR<float,2> const&)>,float,bool,bool);
template OPENSUBDIV_SURFACE<VECTOR<double,3>,3>& STANDARD_TESTS_BASE<VECTOR<double,3> >::Seed_Lagrangian_Particles<OPENSUBDIV_SURFACE<VECTOR<double,3>,3> >(OPENSUBDIV_SURFACE<VECTOR<double,3>,3>&,std::function<VECTOR<double,3> (VECTOR<double,3> const&)>,std::function<MATRIX<double,3,3> (VECTOR<double,3> const&)>,double,bool,bool);
template OPENSUBDIV_SURFACE<VECTOR<float,3>,3>& STANDARD_TESTS_BASE<VECTOR<float,3> >::Seed_Lagrangian_Particles<OPENSUBDIV_SURFACE<VECTOR<float,3>,3> >(OPENSUBDIV_SURFACE<VECTOR<float,3>,3>&,std::function<VECTOR<float,3> (VECTOR<float,3> const&)>,std::function<MATRIX<float,3,3> (VECTOR<float,3> const&)>,float,bool,bool);
template SEGMENTED_CURVE_2D<double>& STANDARD_TESTS_BASE<VECTOR<double,2> >::Seed_Lagrangian_Particles<SEGMENTED_CURVE_2D<double> >(SEGMENTED_CURVE_2D<double>&,std::function<VECTOR<double,2> (VECTOR<double,2> const&)>,std::function<MATRIX<double,2,2> (VECTOR<double,2> const&)>,double,bool,bool);
template SEGMENTED_CURVE_2D<float>& STANDARD_TESTS_BASE<VECTOR<float,2> >::Seed_Lagrangian_Particles<SEGMENTED_CURVE_2D<float> >(SEGMENTED_CURVE_2D<float>&,std::function<VECTOR<float,2> (VECTOR<float,2> const&)>,std::function<MATRIX<float,2,2> (VECTOR<float,2> const&)>,float,bool,bool);
template TRIANGULATED_SURFACE<double>& STANDARD_TESTS_BASE<VECTOR<double,3> >::Seed_Lagrangian_Particles<TRIANGULATED_SURFACE<double> >(TRIANGULATED_SURFACE<double>&,std::function<VECTOR<double,3> (VECTOR<double,3> const&)>,std::function<MATRIX<double,3,3> (VECTOR<double,3> const&)>,double,bool,bool);
template TRIANGULATED_SURFACE<float>& STANDARD_TESTS_BASE<VECTOR<float,3> >::Seed_Lagrangian_Particles<TRIANGULATED_SURFACE<float> >(TRIANGULATED_SURFACE<float>&,std::function<VECTOR<float,3> (VECTOR<float,3> const&)>,std::function<MATRIX<float,3,3> (VECTOR<float,3> const&)>,float,bool,bool);

template void STANDARD_TESTS_BASE<VECTOR<double,3> >::Seed_Particles_Surface<TRIANGULATED_SURFACE<double> >(TRIANGULATED_SURFACE<double> const&,IMPLICIT_OBJECT<VECTOR<double,3> >&,std::function<VECTOR<double,3> (VECTOR<double,3> const&)>,std::function<MATRIX<double,3,3> (VECTOR<double,3> const&)>,double,double);
template void STANDARD_TESTS_BASE<VECTOR<double,3> >::Seed_Particles_Volume<TETRAHEDRALIZED_VOLUME<double> >(TETRAHEDRALIZED_VOLUME<double> &,std::function<VECTOR<double,3> (VECTOR<double,3> const&)>,std::function<MATRIX<double,3,3> (VECTOR<double,3> const&)>,double);
template void STANDARD_TESTS_BASE<VECTOR<float,3> >::Seed_Particles_Surface<TRIANGULATED_SURFACE<float> >(TRIANGULATED_SURFACE<float> const&,IMPLICIT_OBJECT<VECTOR<float,3> >&,std::function<VECTOR<float,3> (VECTOR<float,3> const&)>,std::function<MATRIX<float,3,3> (VECTOR<float,3> const&)>,float,float);
template void STANDARD_TESTS_BASE<VECTOR<float,3> >::Seed_Particles_Volume<TETRAHEDRALIZED_VOLUME<float> >(TETRAHEDRALIZED_VOLUME<float> &,std::function<VECTOR<float,3> (VECTOR<float,3> const&)>,std::function<MATRIX<float,3,3> (VECTOR<float,3> const&)>,float);
}
