//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Vectors/VECTOR.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Fourier_Transforms/FFT_2D.h>
#include <Grid_Tools/Fourier_Transforms/FFTW.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <fftw3.h>

using namespace PhysBAM;

typedef double RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

template<class TV>
struct FOURIER_EXAMPLE: public MPM_EXAMPLE<TV>
{
    FOURIER_EXAMPLE()
        :MPM_EXAMPLE<TV>(STREAM_TYPE((RW)0))
    {
    }

    void Initialize() override {}
};

int main(int argc, char* argv[])
{
#ifdef USE_OPENMP
    omp_set_num_threads(1);
#endif

    bool use_affine=false;
    T flip=0;
    int order=3;
    int resolution=16;
    int particles_per_dim=2;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-resolution",&resolution,"num","resolution");
    parse_args.Add("-affine",&use_affine,"use affine transfers");
    parse_args.Add("-flip",&flip,"num","flip ratio");
    parse_args.Add("-order",&order,"order","interpolation order");
    parse_args.Add("-ppd",&particles_per_dim,"num","particles per cell per dimension");
    parse_args.Parse();

    FOURIER_EXAMPLE<TV> example;
    MPM_DRIVER<TV> driver(example);

    TV_INT center=TV_INT()+resolution/2;

    example.grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
    GRID<TV> init_grid(TV_INT()+particles_per_dim*resolution,RANGE<TV>::Unit_Box(),true);
    for(CELL_ITERATOR<TV> it(init_grid);it.Valid();it.Next()){
        int p=example.particles.Add_Element();
        example.particles.X(p)=it.Location();
        example.particles.mass(p)=1;
        example.particles.valid(p)=true;}
    example.velocity_new.Resize(example.grid.Cell_Indices(3));
    example.velocity.Resize(example.grid.Cell_Indices(3));
    example.gather_scatter.Prepare_Scatter(example.particles);
    example.use_affine=use_affine;
    example.flip=flip;
    example.velocity_new(center).x=1;
    example.velocity_friction=example.velocity_new;
    example.dt=0;
    example.Set_Weights(order);
    driver.Update_Particle_Weights();
    driver.Grid_To_Particle();
    driver.Particle_To_Grid();

    ARRAY<std::complex<T>,TV_INT> row(example.grid.Cell_Indices(0));
    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next())
    {
        TV_INT new_index=it.index+center;
        for(int i=0;i<TV::m;i++) new_index(i)%=resolution;
        row(it.index)=example.velocity(new_index).x;
    }
    
    ARRAY<std::complex<T>,TV_INT> out(example.grid.Cell_Indices(0));

    fftw_plan plan = fftw_plan_dft_2d(example.grid.counts.x, example.grid.counts.y,
        (fftw_complex*)row.array.base_pointer, (fftw_complex*)out.array.base_pointer, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

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
    
    ARRAY<VECTOR<T,3>,TV_INT> image(example.grid.Cell_Indices(0));
    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next())
        image(it.index)=icm(abs(out(it.index)));

    PNG_FILE<T>::Write("eigen.png",image);

    return 0;
}

