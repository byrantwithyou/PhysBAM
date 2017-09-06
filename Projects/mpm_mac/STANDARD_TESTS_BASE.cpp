//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Fourier_Transforms/FFT.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
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
#include <Hybrid_Methods/Seeding/MPM_PARTICLE_SOURCE.h>
#include <fstream>
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
    m(1),s(1),kg(1),forced_collision_type(-1),dump_collision_objects(false),
    test_diff(false),bc_periodic(false),use_analytic_field(false),mu(0),poisson_disk(*new POISSON_DISK<TV>(1))

{
    T framerate=0;
    bool use_quasi_exp_F_update=false;
    bool use_separate=false,use_slip=false,use_stick=false;
    bool print_stats=false;
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Add("-restart",&restart,"frame","restart frame");
    parse_args.Add("-resolution",&resolution,&user_resolution,"resolution","grid resolution");
    parse_args.Add("-substeps",&write_substeps_level,"level","output-substep level");
    parse_args.Add("-substeps_delay",&substeps_delay_frame,"frame","delay substeps until after this frame");
    parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Add("-o",&output_directory,&override_output_directory,"dir","Output directory");
    parse_args.Add("-framerate",&framerate,"rate","Number of frames per second");
    parse_args.Add("-frame_dt",&frame_dt,"rate","Number of frames per second");
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
    parse_args.Add("-lag_Dp",&lag_Dp,"use early gradient transfer for Cp");
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
    parse_args.Add_Not("-use_mass",&use_constant_density,"use particle volumes to avoid boiling");
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
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-analyze_u_modes",&analyze_u_modes,"Perform FFT analysis on velocity");
    parse_args.Add("-dump_modes_freq",&dump_modes_freq,"num","Dump FFT modes every num time steps");
    parse_args.Add("-max_ke",&max_ke,"value","Normalization for FFT images");
    parse_args.Add("-extrap",&extrap_type,&use_extrap,"type","Velocity extrapolation");
    parse_args.Add("-clamp",&clamp_particles,"clamp particles on wall");

    parse_args.Parse(true);
    PHYSBAM_ASSERT((int)use_slip+(int)use_stick+(int)use_separate<=1);
    if(use_slip) forced_collision_type=COLLISION_TYPE::slip;
    if(use_stick) forced_collision_type=COLLISION_TYPE::stick;
    if(use_separate) forced_collision_type=COLLISION_TYPE::separate;

    unit_p=kg*pow<2-TV::m>(m)/(s*s);
    unit_rho=kg*pow<-TV::m>(m);
    unit_mu=kg*pow<2-TV::m>(m)/s;
    min_dt*=s;
    max_dt*=s;
    mu*=unit_mu;

    if(framerate) frame_dt=1/framerate;
    frame_dt*=s;

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

    if(analyze_u_modes){
        Add_Callbacks(false,"p2g",[this](){
                static bool done=false;
                if(!done) Velocity_Fourier_Analysis();
                done=true;});
        Add_Callbacks(false,"time-step",[this](){
                Velocity_Fourier_Analysis();});}

    if(use_periodic_test_shift){
        auto shift_func=[this](int sign){
                TV shift=TV(sign*periodic_test_shift)*grid.dX;
                for(int i=0;i<particles.X.m;i++)
                    particles.X(i)=wrap(particles.X(i)+shift,grid.domain.min_corner,grid.domain.max_corner);};
        Add_Callbacks(true,"time-step",[=](){shift_func(1);});
        Add_Callbacks(false,"time-step",[=](){shift_func(-1);});
        Add_Callbacks(false,"initialize",[=](){
                RANDOM_NUMBERS<T> local_random;
                for(int i=0;i<TV::m;i++)
                    periodic_test_shift(i)=local_random.Get_Uniform_Integer(0,grid.numbers_of_cells(i)-1);
                LOG::printf("shift %P\n",periodic_test_shift);});}

    if(print_stats){
        auto stats_p=[=](const char* name){Add_Callbacks(false,name,[=](){Print_Particle_Stats(name);});};
        auto stats_g=[=](const char* name){Add_Callbacks(false,name,[=](){Print_Grid_Stats(name);});};
        stats_p("simulated-particles");
        stats_g("p2g");
        stats_g("forces");
        stats_g("projection");
        stats_g("viscosity");
        stats_p("g2p");}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
~STANDARD_TESTS_BASE()
{
    analytic_velocity.Delete_Pointers_And_Clean_Memory();
    analytic_pressure.Delete_Pointers_And_Clean_Memory();
    for(int i=0;i<destroy.m;i++) destroy(i)();
    delete &poisson_disk;
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Seed_Particles_Poisson(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
{
    ARRAY<TV> X;
    for(int i=0;i<TV::m;i++) poisson_disk.is_periodic(i)=(bc_type(i)==BC_PERIODIC);
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
    Add_Particle(X,V?V(X):TV(),dV?dV(X):MATRIX<T,TV::m>(),mass,volume);
}
//#####################################################################
// Function Add_Particle
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Particle(const TV& X,const TV& V,const MATRIX<T,TV::m>& dV,const T mass,const T volume)
{
    int p=particles.Add_Element_From_Deletion_List();
    particles.valid(p)=true;
    particles.X(p)=X;
    particles.V(p)=V;
    particles.F(p)=MATRIX<T,TV::m>()+1;
    if(particles.store_Fp) particles.Fp(p).Set_Identity_Matrix();
    if(particles.store_B){
        if(lag_Dp) particles.B(p)=dV;
        else
            for(int a=0;a<TV::m;a++)
                particles.B(p).Set_Row(a,weights(a)->Dp_Inverse(X).Inverse_Times(dV.Row(a)));}
    if(particles.store_S) particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
    particles.mass(p)=mass;
    particles.volume(p)=volume;
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    (*color_attribute)(p)=VECTOR<T,3>(1,1,1);
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
//#####################################################################
// Function Check_Analytic_Velocity
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Check_Analytic_Velocity() const
{
    if(!use_analytic_field) return;
    T max_error=0,l2_error=0;
    int num_l2_samples=0;
    for(PHASE_ID i(0);i<phases.m;i++){
        const PHASE& ph=phases(i);
        for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
            if(ph.mass(it.Full_Index())){
                T u=ph.velocity(it.Full_Index());
                TV v=analytic_velocity(i)->v(it.Location(),time);
                T e=abs(u-v(it.face.axis));
                max_error=std::max(max_error,e);
                l2_error+=sqr(e);
                num_l2_samples++;}}}
    if(num_l2_samples) l2_error/=num_l2_samples;
    l2_error=sqrt(l2_error);
    LOG::printf("grid velocity error: L-inf: %g L-2: %g N: %i\n",max_error,l2_error,num_l2_samples);
    max_error=0;
    l2_error=0;
    num_l2_samples=0;
    for(int p=0;p<particles.X.m;p++){
        if(!particles.valid(p)) continue;
        PHASE_ID pid(0);
        if(particles.store_phase) pid=particles.phase(p);
        TV e=particles.V(p)-analytic_velocity(pid)->v(particles.X(p),time);
        max_error=std::max(max_error,e.Max_Abs());
        l2_error+=e.Magnitude_Squared();
        num_l2_samples+=TV::m;
        PHYSBAM_ASSERT(l2_error<=num_l2_samples*max_error*max_error+1e-10);
    }
    if(num_l2_samples) l2_error/=num_l2_samples;
    l2_error=sqrt(l2_error);
    LOG::printf("particle velocity error: L-inf: %g L-2: %g N: %i\n",max_error,l2_error,num_l2_samples);
}
//#####################################################################
// Function Compute_Analytic_Force
//#####################################################################
template<class TV> TV STANDARD_TESTS_BASE<TV>::
Compute_Analytic_Force(PHASE_ID p,const TV& X,T time) const
{
    if(!use_analytic_field) return TV();
    const ANALYTIC_VECTOR<TV>* av=analytic_velocity(p);
    const ANALYTIC_SCALAR<TV>* ap=analytic_pressure(p);
    T mu=phases(p).viscosity;
    T density=phases(p).density;
    return av->dt(X,time)+av->dX(X,time)*av->v(X,time)-mu/density*av->L(X,time)+1/density*ap->dX(X,time);
}
//#####################################################################
// Function Dump_Image
//#####################################################################
template<class T> static void
Dump_Image(const std::string& file,const ARRAY<T,VECTOR<int,2> >& ke)
{
    INTERPOLATED_COLOR_MAP<T> icm;
    icm.Initialize_Colors((T)1e-8,1,true,true,false);
    
    ARRAY<VECTOR<T,3>,VECTOR<int,2> > image(ke.domain);
    for(RANGE_ITERATOR<2> it(image.domain);it.Valid();it.Next())
        image(it.index)=icm(ke(it.index));

    ARRAY<VECTOR<T,3>,VECTOR<int,2> > bar(VECTOR<int,2>(1000,1));
    for(int i=0;i<1000;i++)
        bar(VECTOR<int,2>(i,0))=icm.colors.Value((i/(T)999)*(icm.mx-icm.mn)+icm.mn);
    PNG_FILE<T>::Write("bar.png",bar);

    PNG_FILE<T>::Write(file,image);
}
//#####################################################################
// Function Dump_Image
//#####################################################################
template<class T,int d> static void
Dump_Image(const std::string& file,const ARRAY<T,VECTOR<int,d> >& ke)
{
}
//#####################################################################
// Function Velocity_Fourier_Analysis
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Velocity_Fourier_Analysis() const
{
    static int id=-1;
    id++;
    LOG::printf("analysis id %i\n",id);
    FFT<TV> fft;
    ARRAY<T,TV_INT> ua(grid.Domain_Indices()),ke(ua.domain);
    ARRAY<std::complex<T>,TV_INT> out(ua.domain);
    const PHASE& ph=phases(PHASE_ID());
    
    TV coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    TV_INT counts=grid.numbers_of_cells,hi=counts/2,lo=hi-counts;
    ARRAY<T> bins(rint((coefficients*TV((hi-1).Componentwise_Max(-lo))).Magnitude())+1);
    int taylor_modes=extra_int.m>=1?extra_int(0):1;
    T total_taylor=0;
    for(int a=0;a<TV::m;a++){
        ua.Put(ph.velocity.Component(a),ua);
        fft.Transform(ua,out);
        out/=sqrt((T)ke.domain.Size());
        for(RANGE_ITERATOR<TV::m> it(ke.domain);it.Valid();it.Next()){
            TV k=coefficients*TV(wrap(it.index,lo,hi));
            ke(wrap(it.index+counts/2,TV_INT(),counts))=sqr(abs(out(it.index))/max_ke);
            bins(rint(k.Magnitude()))+=ke(it.index);}
        if(id%dump_modes_freq==0)
            Dump_Image(LOG::sprintf("%s/modes-%i-%c.png",output_directory,id,"xyz"[a]),ke);
        if(taylor_modes>0)
            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>::Unit_Box()*2);it.Valid();it.Next())
                total_taylor+=sqr(abs(out(wrap((it.index*2-1)*taylor_modes,TV_INT(),grid.numbers_of_cells))));
        else total_taylor+=sqr(abs(out(TV_INT())));}
    std::ofstream fout(LOG::sprintf("%s/bins-%i.txt",output_directory,id).c_str());
    for(int i=0;i<bins.m;i++)
        fout<<i<<" "<<bins(i)<<std::endl;

    T l2_u=0;
    int num_l2_u=0;
    auto valid=[&](FACE_INDEX<TV::m> face){return ph.mass(face) && (this->psi_N.domain_indices.Empty() || !this->psi_N(face));};
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        if(bc_type(2*it.face.axis)==BC_PERIODIC)
            if(it.face.index(it.face.axis)==grid.numbers_of_cells(it.face.axis))
                continue;
        if(!valid(it.face)) continue;
        l2_u+=sqr(ph.velocity(it.face));
        num_l2_u++;}
    if(num_l2_u) l2_u/=num_l2_u;
    if(num_l2_u) total_taylor/=num_l2_u;

    T l2_omega=0;
    int num_l2_omega=0;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        bool ok=true;
        MATRIX<T,TV::m> dV;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++){
                if(i==j) continue;
                FACE_INDEX<TV::m> A(i,it.index),B(A),C(A),D(A);
                C.index(i)++;
                D.index(i)++;
                A.index(j)++;
                B.index(j)--;
                C.index(j)++;
                D.index(j)--;
                if(valid(A) && valid(B) && valid(C) && valid(D))
                    dV(i,j)=(T).25*grid.one_over_dX(j)*(
                        ph.velocity(A)-ph.velocity(B)+
                        ph.velocity(C)-ph.velocity(D));
                else ok=false;}
        if(!ok) continue;
        l2_omega+=dV.Contract_Permutation_Tensor().Magnitude_Squared();
        num_l2_omega++;}
    if(num_l2_omega) l2_omega/=num_l2_omega;

    LOG::printf("l2 velocity %P  l2 vorticity %P\n",l2_u,l2_omega);
    LOG::printf("taylor total %P\n",total_taylor);
}
//#####################################################################
// Function Add_Source
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Source(const TV& X0,const TV& n,IMPLICIT_OBJECT<TV>* io,
    std::function<void(TV X,T ts,T t,SOURCE_PATH<TV>& p)> path,T density,
    T particles_per_cell,bool owns_io)
{
    MPM_PARTICLE_SOURCE<TV>* source=new MPM_PARTICLE_SOURCE<TV>(
        poisson_disk,random,X0,n,io,path);
    T volume=grid.dX.Product()/particles_per_cell;
    T mass=density*volume;
    // NOTE: assumes initial particles are already added.
    source->Seed_Points(particles.X);
    Add_Callbacks(false,"time-step",[=](){
            ARRAY<TV> X,V;
            ARRAY<MATRIX<T,TV::m> > dV;
            source->Seed(time,dt,X,V,&dV);
            for(int i=0;i<X.m;i++)
                Add_Particle(X(i),V(i),dV(i),mass,volume);});
    if(owns_io) destroy.Append([io](){delete io;});
}
//#####################################################################
// Function Setup_Analytic_Boundary_Conditions
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Setup_Analytic_Boundary_Conditions()
{
    auto bc_v=[this](const TV& X,int axis,PHASE_ID p,T t)
        {return analytic_velocity(p)->v(X,t)(axis);};
    for(int i=0;i<bc_velocity.m;i++)
        if(bc_type(i)==BC_WALL)
            bc_velocity(i)=bc_v;
    bc_pressure=[this](TV_INT c,PHASE_ID p,T t)
        {return analytic_pressure(p)->f(grid.Center(c),t);};
}
template class STANDARD_TESTS_BASE<VECTOR<float,2> >;
template class STANDARD_TESTS_BASE<VECTOR<float,3> >;
template class STANDARD_TESTS_BASE<VECTOR<double,2> >;
template class STANDARD_TESTS_BASE<VECTOR<double,3> >;
}
