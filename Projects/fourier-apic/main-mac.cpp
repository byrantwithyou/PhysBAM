//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Fourier_Transforms/FFT.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <fftw3.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

template<class TV>
struct FOURIER_EXAMPLE: public MPM_MAC_EXAMPLE<TV>
{
    FOURIER_EXAMPLE()
        :MPM_MAC_EXAMPLE<TV>(STREAM_TYPE((RW)0))
    {
    }

    void Initialize() override {}
};

// X is particle seeding in unit box; replicate it in every grid cell
void Replicate_Particles(MPM_PARTICLES<TV>& particles,const GRID<TV>& grid,const ARRAY<TV>& X)
{
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        TV offset=grid.Node(it.index);
        for(int i=0;i<X.m;i++){
            int p=particles.Add_Element();
            particles.X(p)=offset+grid.dX*X(i);
            particles.mass(p)=1;
            particles.valid(p)=true;}}
}

void Sample_Box_Regularly(ARRAY<TV>& X,int particles_per_dim)
{
    GRID<TV> init_grid(TV_INT()+particles_per_dim,RANGE<TV>::Unit_Box(),true);
    for(CELL_ITERATOR<TV> it(init_grid);it.Valid();it.Next())
        X.Append(it.Location());
}

void Sample_Box_Random(RANDOM_NUMBERS<T>& rand,ARRAY<TV>& X,int number_of_particles)
{
    for(int i=0;i<number_of_particles;i++)
        X.Append(rand.Get_Uniform_Vector(RANGE<TV>::Unit_Box()));
}

int main(int argc, char* argv[])
{
#ifdef USE_OPENMP
    omp_set_num_threads(1);
#endif

    bool use_affine=false;
    T flip=0;
    int order=3;
    int resolution=16;
    int size=256;
    int particles_per_dim=2;
    int irregular_seeding=0;
    int seed=-1;
    std::string output_filename="eigen.png";
    std::string viewer_directory="output";
    bool dump_particles=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-resolution",&resolution,"num","transfer resolution");
    parse_args.Add("-size",&size,"num","analyze transfer as though this resolution");
    parse_args.Add("-affine",&use_affine,"use affine transfers");
    parse_args.Add("-flip",&flip,"num","flip ratio");
    parse_args.Add("-order",&order,"order","interpolation order");
    parse_args.Add("-ppd",&particles_per_dim,"num","particles per cell per dimension");
    parse_args.Add("-o",&output_filename,"file.png","filename for output image");
    parse_args.Add("-v",&viewer_directory,"dir","viewer directory");
    parse_args.Add("-irreg",&irregular_seeding,"num","each cell is seeded identicially with num particles");
    parse_args.Add("-seed",&seed,"seed","random number generator seed (-1 = timer)");
    parse_args.Add("-dump_particles",&dump_particles,"Output particle distribution");
    parse_args.Parse();

    PHYSBAM_ASSERT(resolution<=size);
    
    FOURIER_EXAMPLE<TV> example;
    MPM_MAC_DRIVER<TV> driver(example);

    TV_INT center=TV_INT()+resolution/2;
    RANDOM_NUMBERS<T> rand;
    if(seed!=-1) rand.Set_Seed(seed);
    
    example.grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
    ARRAY<TV> unit_X;
    if(irregular_seeding) Sample_Box_Random(rand,unit_X,irregular_seeding);
    else Sample_Box_Regularly(unit_X,particles_per_dim);
    Replicate_Particles(example.particles,example.grid,unit_X);

    example.dt=0;
    example.flip=flip;
    example.use_affine=use_affine;
    example.Set_Weights(order);
    example.phases.Resize(PHASE_ID(1));
    auto& ph=example.phases(PHASE_ID());
    ph.density=1;
    ph.viscosity=0;
    ph.Initialize(example.grid,example.weights,example.ghost,example.threads);
    ph.velocity.Resize(example.grid.Cell_Indices(3));
    ph.mass.Resize(example.grid.Cell_Indices(3));
    ph.gather_scatter->Prepare_Scatter(example.particles);
    example.periodic_boundary.is_periodic.Fill(true);
    example.bc_type.Fill(example.BC_PERIODIC);
    example.particles.Store_B(example.use_affine);
    example.particles.Store_C(false);
    ph.velocity(FACE_INDEX<TV::m>(0,center))=1;
    driver.Update_Simulated_Particles();
    driver.Update_Particle_Weights();
    driver.Grid_To_Particle();
    driver.Particle_To_Grid();

    ARRAY<T,TV_INT> row(TV_INT()+size);
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>::Unit_Box()*resolution);it.Valid();it.Next()){
        TV_INT index=it.index-center;
        for(int i=0;i<TV::m;i++) if(index(i)<0) index(i)+=size;
        row(index)=ph.velocity(FACE_INDEX<TV::m>(0,it.index));
        PHYSBAM_ASSERT(!ph.velocity(FACE_INDEX<TV::m>(1,it.index)));}

    ARRAY<std::complex<T>,TV_INT> out(row.domain);

    GRID<TV> fft_grid(out.domain.Edge_Lengths(),RANGE<TV>::Unit_Box(),true);
    FFT<TV> fft;
    fft.Transform(row,out);
    
    INTERPOLATED_COLOR_MAP<T> icm;
    icm.colors.Add_Control_Point(1.00001,VECTOR<T,3>(1,1,1));
    icm.colors.Add_Control_Point(1,VECTOR<T,3>(.5,0,0));
    icm.colors.Add_Control_Point(1-.01,VECTOR<T,3>(1,0,0));
    icm.colors.Add_Control_Point(1-.02,VECTOR<T,3>(1,.5,0));
    icm.colors.Add_Control_Point(1-.04,VECTOR<T,3>(1,1,0));
    icm.colors.Add_Control_Point(1-.08,VECTOR<T,3>(0,1,0));
    icm.colors.Add_Control_Point(1-.16,VECTOR<T,3>(0,1,1));
    icm.colors.Add_Control_Point(1-.32,VECTOR<T,3>(0,0,1));
    icm.colors.Add_Control_Point(1-.64,VECTOR<T,3>(.5,0,1));
    icm.colors.Add_Control_Point(0,VECTOR<T,3>(0,0,0));

    ARRAY<VECTOR<T,3>,TV_INT> image(out.domain);
    for(RANGE_ITERATOR<TV::m> it(image.domain);it.Valid();it.Next())
        image(it.index)=icm(abs(out(it.index)));

    PNG_FILE<T>::Write(output_filename,image);

    if(dump_particles){
        VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),example.grid,viewer_directory);
        for(int i=0;i<example.particles.X.m;i++)
            Add_Debug_Particle(example.particles.X(i),VECTOR<T,3>(1,0,0));
        Flush_Frame<TV>("particles");
        Flush_Frame<TV>("end");
    }
    
    return 0;
}

