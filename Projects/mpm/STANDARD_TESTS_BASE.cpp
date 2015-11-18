//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Constitutive_Models/ST_VENANT_KIRCHHOFF_HENCKY_STRAIN.h>
#include <Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FINITE_ELEMENTS.h>
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
    :MPM_EXAMPLE<TV>(stream_type_input),test_number(0),resolution(32),user_resolution(false),stored_last_frame(0),
    user_last_frame(false),order(2),seed(1234),particles_per_cell(1<<TV::m),regular_seeding(false),
    scale_mass(1),scale_E(1),scale_speed(1),
    penalty_collisions_stiffness((T)1e4),penalty_collisions_separation((T)1e-4),penalty_collisions_length(1),
    penalty_damping_stiffness(0),tests(stream_type_input,deformable_body_collection)
{
    T framerate=24;
    bool use_quasi_exp_F_update=false;
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Add("-restart",&restart,"frame","restart frame");
    parse_args.Add("-resolution",&resolution,&user_resolution,"resolution","grid resolution");
    parse_args.Add("-substeps",&write_substeps_level,"level","output-substep level");
    parse_args.Add("-substeps_delay",&substeps_delay_frame,"frame","delay substeps until after this frame");
    parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
    parse_args.Add("-test_diff",&test_diff,"test analytic derivatives");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Add("-o",&output_directory,"dir","Output directory");
    parse_args.Add("-mass_contour",&mass_contour,"contour","Draw mass contour as a scale to particle average mass");
    parse_args.Add("-framerate",&framerate,"rate","Number of frames per second");
    parse_args.Add("-min_dt",&min_dt,"dt","Minimum time step size");
    parse_args.Add("-max_dt",&max_dt,"dt","Maximum time step size");
    parse_args.Add("-order",&order,"order","Interpolation basis order");
    parse_args.Add("-affine",&use_affine,"Use affine PIC");
    parse_args.Add("-use_f2p",&use_f2p,"Use face to particle (direct) transfer");
    parse_args.Add("-midpoint",&use_midpoint,"Use midpoint rule");
    parse_args.Add("-symplectic_euler",&use_symplectic_euler,"Use forward euler for grid update");
    parse_args.Add("-particle_collision",&use_particle_collision,"Use particle collision");
    parse_args.Add("-print_stats",&print_stats,"Print momentum/energy stats");
    parse_args.Add("-flip",&flip,"frac","Flip ratio");
    parse_args.Add("-cfl",&cfl,"cfl","CFL number");
    parse_args.Add("-newton_tolerance",&newton_tolerance,"tol","Newton tolerance");
    parse_args.Add("-newton_iterations",&newton_iterations,"iter","Newton iterations");
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
    parse_args.Add("-use_early_gradient_transfer",&use_early_gradient_transfer,"use early gradient transfer for Cp");
    parse_args.Add("-incompressible",&incompressible,"Make simulated media incompressible");
    parse_args.Add("-kkt",&kkt,"Use KKT solver");
    parse_args.Add("-use_exp_F",&use_quasi_exp_F_update,"Use an approximation of the F update that prevents inversion");
    parse_args.Add("-use_plasticity",&use_plasticity,"Use plasticity in the F update");
    parse_args.Add("-theta_c",&theta_c,"theta_c","Critical compression coefficient for plasticity");
    parse_args.Add("-theta_s",&theta_s,"theta_s","Critical stretch coefficient for plasticity");
    parse_args.Add("-hardening",&hardening_factor,"hardening factor","Hardening factor for plasticity");
    parse_args.Add("-max_hardening",&max_hardening,"max hardening coefficient","Maximum hardening coefficient for plasticity");

    parse_args.Parse(true);

    frame_dt=1/framerate;

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
    particles.Store_B(use_affine && !incompressible);
    particles.Store_C(use_affine && (incompressible || kkt));
    particles.Store_One_Over_Lambda(kkt);

    if(order==1) Set_Weights(new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(grid,threads));
    else if(order==2) Set_Weights(new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(grid,threads));
    else if(order==3) Set_Weights(new PARTICLE_GRID_WEIGHTS_SPLINE<TV,3>(grid,threads));
    else PHYSBAM_FATAL_ERROR("Unrecognized interpolation order");
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
~STANDARD_TESTS_BASE()
{
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Seed_Particles(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell)
{
    POISSON_DISK<TV> poisson_disk(1);
    ARRAY<TV> X;
    poisson_disk.Set_Distance_By_Volume(grid.dX.Product()/particles_per_cell);
    poisson_disk.Sample(random,object,X);

    object.Update_Box();
    RANGE<TV_INT> range=grid.Cell_Indices(0).Intersect(grid.Clamp_To_Cell(object.Box(),3));

    T volume=grid.dX.Product()/particles_per_cell;
    T mass=density*volume;
    for(int i=0;i<X.m;i++)
        Add_Particle(X(i),V,dV,mass,volume);
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Seed_Particles(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,const GRID<TV>& seed_grid)
{
    T volume=seed_grid.dX.Product();
    T mass=density*volume;
    for(CELL_ITERATOR<TV> it(seed_grid);it.Valid();it.Next()){
        TV X=it.Location();
        if(object.Lazy_Inside(X))
            Add_Particle(X,V,dV,mass,volume);}
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Seed_Particles_Helper(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell)
{
    if(!regular_seeding) return Seed_Particles(object,V,dV,density,particles_per_cell);

    object.Update_Box();
    RANGE<TV_INT> range=grid.Clamp_To_Cell(object.Box(),3).Intersect(grid.Cell_Indices());
    TV LB=grid.Node(range.min_corner);
    TV UB=grid.Node(range.max_corner);
    T scale=pow<1,TV::m>((T)particles_per_cell);
    GRID<TV> seed_grid(range.Edge_Lengths()*scale,RANGE<TV>(LB,UB),true);
    Seed_Particles(object,V,dV,density,seed_grid);
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
    if(particles.store_B && dV) particles.B(p)=dV(X)*weights->Dp(X);
    if(particles.store_C && dV) particles.C(p)=dV(X);
    if(particles.store_S) particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
    particles.mass(p)=mass;
    particles.volume(p)=volume;
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    (*color_attribute)(p)=VECTOR<T,3>(1,1,1);
}
//#####################################################################
// Function Add_Gravity
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Gravity(TV g)
{
    return Add_Force(*new DEFORMABLE_GRAVITY<TV>(particles,true,g));
}
//#####################################################################
// Function Add_Fixed_Corotated
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Fixed_Corotated(T E,T nu,ARRAY<int>* affected_particles,bool no_mu)
{
    COROTATED_FIXED<T,TV::m>* coro=new COROTATED_FIXED<T,TV::m>(E,nu);
    if(no_mu) coro->Zero_Out_Mu();
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model=*coro;
    MPM_FINITE_ELEMENTS<TV>& fe=*new MPM_FINITE_ELEMENTS<TV>(force_helper,constitutive_model,gather_scatter,affected_particles,use_variable_coefficients);
    return Add_Force(fe);
}
//#####################################################################
// Function Add_Neo_Hookean
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_Neo_Hookean(T E,T nu,ARRAY<int>* affected_particles)
{
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model=*new NEO_HOOKEAN<T,TV::m>(E,nu);
    MPM_FINITE_ELEMENTS<TV>& fe=*new MPM_FINITE_ELEMENTS<TV>(force_helper,constitutive_model,gather_scatter,affected_particles,use_variable_coefficients);
    return Add_Force(fe);
}
//#####################################################################
// Function Add_St_Venant_Kirchhoff_Hencky_Strain
//#####################################################################
template<class TV> int STANDARD_TESTS_BASE<TV>::
Add_St_Venant_Kirchhoff_Hencky_Strain(T E,T nu,ARRAY<int>* affected_particles,bool no_mu)
{
    ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,TV::m>* hencky=new ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,TV::m>(E,nu);
    if(no_mu) hencky->Zero_Out_Mu();
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model=*hencky;
    MPM_FINITE_ELEMENTS<TV>& fe=*new MPM_FINITE_ELEMENTS<TV>(force_helper,constitutive_model,gather_scatter,affected_particles,use_variable_coefficients);
    return Add_Force(fe);
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
    for(int p=old_particles_number;p<particles.number;p++){
        TV X=particles.X(p);
        particles.valid(p)=true;
        if(V) particles.V(p)=V(X);
        if(particles.store_B && dV) particles.B(p)=dV(X)*weights->Dp(X);
        if(particles.store_C && dV) particles.C(p)=dV(X);
        particles.F(p)=MATRIX<T,TV::m>()+1;
        if(particles.store_S) particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
        particles.volume(p)=particles.mass(p)/density;}
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
Add_Penalty_Collision_Object(IMPLICIT_OBJECT<TV>* io)
{
    this->Add_Force(*new IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>(particles,io,penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length));
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
}
