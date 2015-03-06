//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include "STANDARD_TESTS_BASE.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
STANDARD_TESTS_BASE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
    :MPM_EXAMPLE<TV>(stream_type),test_number(0),resolution(32),stored_last_frame(0),user_last_frame(false),
    order(2),seed(1234),particles_per_cell(1<<TV::m),scale_mass(1),scale_E(1)
{
    T framerate=24;
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Add("-restart",&restart,"frame","restart frame");
    parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
    parse_args.Add("-substeps",&write_substeps_level,"level","output-substep level");
    parse_args.Add("-substeps_delay",&substeps_delay_frame,"frame","delay substeps until after this frame");
    parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
    parse_args.Add("-test_diff",&test_diff,"test analytic derivatives");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Add("-o",&output_directory,"dir","Output directory");
    parse_args.Add("-framerate",&framerate,"rate","Number of frames per second");
    parse_args.Add("-min_dt",&min_dt,"dt","Minimum time step size");
    parse_args.Add("-max_dt",&max_dt,"dt","Maximum time step size");
    parse_args.Add("-order",&order,"order","Interpolation basis order");
    parse_args.Add("-affine",&use_affine,"Use affine PIC");
    parse_args.Add("-midpoint",&use_midpoint,"Use midpoint rule");
    parse_args.Add("-forward_euler",&use_forward_euler,"Use forward euler for grid update");
    parse_args.Add("-particle_collision",&use_particle_collision,"Use particle collision");
    parse_args.Add("-print_stats",&print_stats,"Print momentum/energy stats");
    parse_args.Add("-flip",&flip,"frac","Flip ratio");
    parse_args.Add("-cfl",&cfl,"cfl","CFL number");
    parse_args.Add("-newton_tolerance",&newton_tolerance,"tol","Newton tolerance");
    parse_args.Add("-newton_iterations",&newton_iterations,"iter","Newton iterations");
    parse_args.Add("-solver_tolerance",&solver_tolerance,"tol","Solver tolerance");
    parse_args.Add("-solver_iterations",&solver_iterations,"iter","Solver iterations");
    parse_args.Add("-test_diff",&test_diff,"Test derivatives");
    parse_args.Add("-threads",&threads,"num","Number of threads");
    parse_args.Add("-seed",&seed,"seed","Random number seed");
    parse_args.Add("-particles_per_cell",&particles_per_cell,"num","Number of particles per cell");
    parse_args.Add("-scale_mass",&scale_mass,"scale","Scale mass of particles");
    parse_args.Add("-scale_E",&scale_E,"scale","Scale stiffness");
    parse_args.Parse(true);

    frame_dt=1/framerate;

#ifdef USE_OPENMP
    omp_set_num_threads(number_of_threads);
#pragma omp parallel
#pragma omp single
    {
        if(omp_get_num_threads()!=number_of_threads) PHYSBAM_FATAL_ERROR();
        LOG::cout<<"Running on "<<number_of_threads<<" threads"<<std::endl;
    }
#endif

    stored_last_frame=last_frame;
    random.Set_Seed(seed);

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
Seed_Particles(IMPLICIT_OBJECT<TV>& object,boost::function<TV(const TV&)> V,
    boost::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell)
{
    POISSON_DISK<TV> poisson_disk(1);
    ARRAY<TV> X;
    poisson_disk.Set_Distance_By_Volume(grid.dX.Product()/particles_per_cell);
    poisson_disk.Sample(random,object,X);

    object.Update_Box();
    RANGE<TV_INT> range=grid.Cell_Indices(0).Intersect(grid.Clamp_To_Cell(object.Box(),3));

    T volume=grid.dX.Product()/particles_per_cell;
    T mass=density*volume;
    for(int i=0;i<X.m;i++){
        int p=particles.Add_Element();
        particles.valid(p)=true;
        particles.X(p)=X(i);
        particles.V(p)=V(X(i));
        if(use_affine) particles.B(p)=dV(X(i))*weights->Dp(X(i));
        particles.F(p)=MATRIX<T,TV::m>()+1;
        particles.mass(p)=mass;
        particles.volume(p)=volume;}
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
Add_Fixed_Corotated(T E,T nu,ARRAY<int>* affected_particles)
{
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model=*new COROTATED_FIXED<T,TV::m>(E,nu);
    MPM_FINITE_ELEMENTS<TV>& fe=*new MPM_FINITE_ELEMENTS<TV>(particles,constitutive_model,gather_scatter,affected_particles);
    return Add_Force(fe);
}
template class STANDARD_TESTS_BASE<VECTOR<float,2> >;
template class STANDARD_TESTS_BASE<VECTOR<float,3> >;
template class STANDARD_TESTS_BASE<VECTOR<double,2> >;
template class STANDARD_TESTS_BASE<VECTOR<double,3> >;
}
