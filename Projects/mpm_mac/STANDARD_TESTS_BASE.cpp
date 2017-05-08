//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_SYSTEM.h>
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
    :MPM_MAC_EXAMPLE<TV>(stream_type_input),test_number(0),resolution(32),user_resolution(false),stored_last_frame(0),
    user_last_frame(false),order(2),seed(1234),particles_per_cell(1<<TV::m),regular_seeding(false),
    no_regular_seeding(false),scale_mass(1),override_output_directory(false),
    m(1),s(1),kg(1),forced_collision_type(-1),write_output_files(0),read_output_files(0),dump_collision_objects(false),
    test_diff(false),bc_periodic(false),use_periodic_test_shift(false)
{
    T framerate=24;
    bool use_quasi_exp_F_update=false;
    bool no_affine=false;
    bool use_separate=false,use_slip=false,use_stick=false;
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Add("-restart",&restart,"frame","restart frame");
    parse_args.Add("-resolution",&resolution,&user_resolution,"resolution","grid resolution");
    parse_args.Add("-substeps",&write_substeps_level,"level","output-substep level");
    parse_args.Add("-substeps_delay",&substeps_delay_frame,"frame","delay substeps until after this frame");
    parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Add("-o",&output_directory,&override_output_directory,"dir","Output directory");
    parse_args.Add("-framerate",&framerate,"rate","Number of frames per second");
    parse_args.Add("-min_dt",&min_dt,"dt","Minimum time step size");
    parse_args.Add("-max_dt",&max_dt,"dt","Maximum time step size");
    parse_args.Add("-order",&order,"order","Interpolation basis order");
    parse_args.Add_Not("-no_affine",&use_affine,"Use affine PIC");
    parse_args.Add("-affine",&use_affine,"Use affine PIC");
    parse_args.Add("-print_stats",&print_stats,"Print momentum/energy stats");
    parse_args.Add("-only_write_particles",&only_write_particles,"Only write particle data (ignore grid data, restart data etc)");
    parse_args.Add("-cfl",&cfl,"cfl","CFL number");
    parse_args.Add("-solver_tolerance",&solver_tolerance,"tol","Solver tolerance");
    parse_args.Add("-solver_iterations",&solver_iterations,"iter","Solver iterations");
    parse_args.Add("-seed",&seed,"seed","Random number seed");
    parse_args.Add("-particles_per_cell",&particles_per_cell,"num","Number of particles per cell");
    parse_args.Add("-scale_mass",&scale_mass,"scale","Scale mass of particles");
    parse_args.Add("-regular_seeding",&regular_seeding,"use regular particle seeding");
    parse_args.Add_Not("-no_regular_seeding",&no_regular_seeding,"use regular particle seeding");
    parse_args.Add("-use_early_gradient_transfer",&use_early_gradient_transfer,"use early gradient transfer for Cp");
    parse_args.Add("-use_exp_F",&use_quasi_exp_F_update,"Use an approximation of the F update that prevents inversion");
    parse_args.Add("-m",&m,"scale","meter scale");
    parse_args.Add("-s",&s,"scale","second scale");
    parse_args.Add("-kg",&kg,"scale","kilogram scale");
    parse_args.Add("-slip",&use_slip,"force slip collisions");
    parse_args.Add("-stick",&use_stick,"force stick collisions");
    parse_args.Add("-separate",&use_separate,"force separating collisions");
    parse_args.Add("-T",&extra_T,"float","extra float argument");
    parse_args.Add("-I",&extra_int,"int","extra int argument");
    parse_args.Add("-dump_collisions",&dump_collision_objects,"dump out collision objects");
    parse_args.Add("-use_volume",&use_particle_volumes,"use particle volumes to avoid boiling");
    parse_args.Add("-move_mass_inside",&move_mass_inside,"move mass and momentum transferred outside back in");
    parse_args.Add("-move_mass_inside_nearest",&move_mass_inside_nearest,"move mass and momentum transferred outside back in");
    parse_args.Add("-test_system",&test_system,"run sanity checks on Krylov systems");
    parse_args.Add("-print_matrix",&print_matrix,"print Krylov matrix");
    parse_args.Add("-test_output_prefix",&test_output_prefix,&use_test_output,"","prefix to use for test output");
    parse_args.Add_Not("-no_preconditioner",&this->projection_system.use_preconditioner,"disable preconditioner");
    parse_args.Add("-flip",&flip,"frac","Use this fraction of flip in transfers");
    parse_args.Add("-test_diff",&test_diff,"test analytic derivatives");
    parse_args.Add("-bc_periodic",&bc_periodic,"set boundary condition periodic");
    parse_args.Add("-test_periodic",&use_periodic_test_shift,"test periodic bc");

    parse_args.Parse(true);
    PHYSBAM_ASSERT((int)use_slip+(int)use_stick+(int)use_separate<=1);
    if(use_slip) forced_collision_type=COLLISION_TYPE::slip;
    if(use_stick) forced_collision_type=COLLISION_TYPE::stick;
    if(use_separate) forced_collision_type=COLLISION_TYPE::separate;

    if(use_periodic_test_shift){
        auto old_bf=begin_frame,old_ef=end_frame;
        begin_frame=[this,old_bf](int frame)
            {
                if(frame==0)
                    for(int i=0;i<TV::m;i++)
                        periodic_test_shift(i)=random.Get_Uniform_Integer(0,grid.numbers_of_cells(i))*grid.dX(i);
                for(int i=0;i<particles.X.m;i++)
                    particles.X(i)=wrap(particles.X(i)+periodic_test_shift,grid.domain.min_corner,grid.domain.max_corner);
                if(old_bf) old_bf(frame);
            };
        end_frame=[this,old_ef](int frame)
            {
                if(old_ef) old_ef(frame);
                for(int i=0;i<particles.X.m;i++)
                    particles.X(i)=wrap(particles.X(i)-periodic_test_shift,grid.domain.min_corner,grid.domain.max_corner);
            };}

    unit_p=kg*pow<2-TV::m>(m)/(s*s);
    unit_rho=kg*pow<-TV::m>(m);
    unit_mu=kg*pow<2-TV::m>(m)/s;
    min_dt*=s;
    max_dt*=s;

    if(no_affine) use_affine=false;

    framerate/=s;
    frame_dt=1/framerate;

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

    stored_last_frame=last_frame;
    random.Set_Seed(seed);

    particles.Store_B(use_affine);
    particles.Store_C(false);
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
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
{
    POISSON_DISK<TV> poisson_disk(1);
    ARRAY<TV> X;
    poisson_disk.Set_Distance_By_Volume(grid.dX.Product()/particles_per_cell);
    poisson_disk.Sample(random,object,X);

    T volume=grid.dX.Product()/particles_per_cell;
    T mass=density*volume;
    for(int i=0;i<X.m;i++)
        Add_Particle(X(i),V,dV,mass,volume);
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Seed_Particles_Uniform(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
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
Seed_Particles(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
{
    if(test_diff) Test_dV(V,dV);

    if(!regular_seeding) return Seed_Particles_Poisson(object,V,dV,density,particles_per_cell);

    object.Update_Box();
    RANGE<TV_INT> range=grid.Clamp_To_Cell(object.Box(),3).Intersect(grid.Cell_Indices());
    TV LB=grid.Node(range.min_corner);
    TV UB=grid.Node(range.max_corner);
    T scale=pow<1,TV::m>((T)particles_per_cell);
    GRID<TV> seed_grid(range.Edge_Lengths()*scale,RANGE<TV>(LB,UB),true);
    Seed_Particles_Uniform(object,V,dV,density,seed_grid);
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
    if(particles.store_B && dV) particles.B(p)=dV(X)*weights(0)->Dp(X);
    if(particles.store_C && dV) particles.C(p)=dV(X);
    if(particles.store_S) particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
    particles.mass(p)=mass;
    particles.volume(p)=volume;
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    (*color_attribute)(p)=VECTOR<T,3>(1,1,1);
}
//#####################################################################
// Function Add_Walls
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Walls(int flags,COLLISION_TYPE type,T inset) // -x +x -y +y [ -z +z ], as bit flags
{
    RANGE<TV> range=grid.domain.Thickened(grid.dX*(ghost*2+1));
    for(int a=0;a<TV::m;a++)
        for(int s=0;s<2;s++)
            if(flags&(1<<(a*2+s))){
                RANGE<TV> wall=range;
                if(s) wall.max_corner(a)=grid.domain.min_corner(a)+inset;
                else wall.min_corner(a)=grid.domain.max_corner(a)-inset;
                Add_Collision_Object(wall,type,0);}
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
// Function Test_dV
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Test_dV(std::function<TV(const TV&)> V,std::function<MATRIX<T,TV::m>(const TV&)> dV) const
{
    PHYSBAM_ASSERT(V || !dV);
    if(!dV) return;
    RANDOM_NUMBERS<T> rand;
    T eps=(T)1e-6;
    TV x0,dx;
    rand.Fill_Uniform(x0,-(T)1,(T)1);
    rand.Fill_Uniform(dx,-eps,eps);
    TV v0=V(x0),v1=V(x0+dx);
    MATRIX<T,TV::m> dv0=dV(x0),dv1=dV(x0+dx);
    TV a=(v1-v0)/eps,b=(T).5/eps*(dv0+dv1)*dx,c=a-b;
    T ma=a.Magnitude(),mb=b.Magnitude(),mc=c.Magnitude();
    LOG::printf("dV %g %g %g  rel %g\n",ma,mb,mc,mc/max(ma,mb,(T)1e-30));
}
//#####################################################################
// Function Set_Phases
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_Phases(const ARRAY<T,PHASE_ID>& phase_densities)
{
    phases.Resize(phase_densities.m);
    for(PHASE_ID i(0);i<phase_densities.m;i++)
        phases(i).density=phase_densities(i);
}
template class STANDARD_TESTS_BASE<VECTOR<float,2> >;
template class STANDARD_TESTS_BASE<VECTOR<float,3> >;
template class STANDARD_TESTS_BASE<VECTOR<double,2> >;
template class STANDARD_TESTS_BASE<VECTOR<double,3> >;
}
