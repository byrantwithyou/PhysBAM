//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Math_Tools/pow.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <fstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

template<class TV>
struct FOURIER_EXAMPLE:public MPM_MAC_EXAMPLE<TV>
{
    FOURIER_EXAMPLE()
        :MPM_MAC_EXAMPLE<TV>(STREAM_TYPE((RW)0))
    {
    }

    void Initialize() override {}
};

RANDOM_NUMBERS<T> rand_;

template<class T,class TV> void
Add_Particle(const MPM_MAC_EXAMPLE<TV>& example,const TV& X,const TV& V,const MATRIX<T,TV::m>& dV,T mass,T volume)
{
    int p=example.particles.Add_Element_From_Deletion_List();
    example.particles.valid(p)=true;
    example.particles.X(p)=X;
    example.particles.V(p)=V;
    example.particles.F(p)=MATRIX<T,TV::m>()+1;
    if(example.particles.store_Fp) example.particles.Fp(p).Set_Identity_Matrix();
    if(example.particles.store_B)
        for(int a=0;a<TV::m;a++)
            example.particles.B(p).Set_Row(a,example.weights(a)->Dp_Inverse(X).Inverse_Times(dV.Row(a)));
    if(example.particles.store_S) example.particles.S(p)=SYMMETRIC_MATRIX<T,TV::m>()+1;
    example.particles.mass(p)=mass;
    example.particles.volume(p)=volume;
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=example.particles.template Get_Array<VECTOR<T,3> >("color");
    (*color_attribute)(p)=VECTOR<T,3>(1,1,1);
}

template<class T,class TV> void
Seed_Particles(const MPM_MAC_EXAMPLE<TV>& example,std::function<TV(const TV&)> V,
    std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
{
    ARRAY<TV> X;
    POISSON_DISK<TV> poisson_disk(1);
    for(int i=0;i<TV::m;i++) poisson_disk.is_periodic(i)=true;
    poisson_disk.Set_Distance_By_Volume(example.grid.dX.Product()/particles_per_cell);
    poisson_disk.Sample(rand_,*Make_IO(example.grid.domain),X);

    T volume=example.grid.dX.Product()/particles_per_cell;
    T mass=density*volume;
    for(int i=0;i<X.m;i++)
        Add_Particle(example,X(i),V(X(i)),dV(X(i)),mass,volume);
}

int main(int argc, char* argv[])
{
#ifdef USE_OPENMP
    omp_set_num_threads(1);
#endif

    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    int order=3;
    int resolution=16;
    int particles_per_dim=2;
    int seed=-1;
    std::string output_filename="eigen.png";
    std::string viewer_directory="output";
    T inv_rate=0;

    FOURIER_EXAMPLE<TV> example;
    MPM_MAC_DRIVER<TV> driver(example);

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-resolution",&resolution,"num","transfer resolution");
    parse_args.Add_Not("-no_affine",&example.use_affine,"use affine transfers");
    parse_args.Add("-flip",&example.flip,"num","flip ratio");
    parse_args.Add("-order",&order,"order","interpolation order");
    parse_args.Add("-ppd",&particles_per_dim,"num","particles per cell per dimension");
    parse_args.Add("-o",&viewer_directory,"dir","viewer directory");
    parse_args.Add("-seed",&seed,"seed","random number generator seed (-1 = timer)");
    parse_args.Add("-mls",&example.use_mls_xfers,"use affine transfers");
    parse_args.Add("-inv_rate",&inv_rate,"rate","invalidate a sample with probability");
    parse_args.Parse();
    Create_Directory(viewer_directory);

    if(seed!=-1) rand_.Set_Seed(seed);
    
    example.grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),example.grid,viewer_directory);

    Seed_Particles<T,TV>(example,
        [](const TV& X){return TV();},
        [](const TV& X){return MATRIX<T,TV::m>();},
        (T)100,pow<TV::m>((T)particles_per_dim));
//    Flush_Frame<TV>("initial");

    // TODO: sample particles

    LOG::printf("example.use_affine %i\n",example.use_affine);
    
    example.dt=0;
    if(example.use_affine)
        for(int i=0;i<TV::m;i++)
            example.Dp_inv(i).Resize(example.particles.X.m);
    example.Set_Weights(order);
    example.density=1;
    example.viscosity=0;
    example.velocity.Resize(example.grid.Cell_Indices(3));
    example.volume.Resize(example.grid.Cell_Indices(3));
    example.mass.Resize(example.grid.Cell_Indices(3));
    example.valid_xfer_data.Resize(example.grid.Cell_Indices(3));
    example.gather_scatter=new GATHER_SCATTER<TV>(example.grid,example.simulated_particles);
    example.gather_scatter->threads=example.threads;
    example.gather_scatter->face_weights=example.weights;
    example.gather_scatter->weights=example.weights(0);
    example.gather_scatter->Prepare_Scatter(example.particles);
    example.periodic_boundary.is_periodic.Fill(true);
    example.side_bc_type.Fill(example.BC_PERIODIC);
    example.particles.Store_B(example.use_affine);
    driver.Update_Simulated_Particles();
    driver.Update_Particle_Weights();

    driver.Write_Output_Files(0);

    TV A(.2,-.1);
    MATRIX<T,TV::m> M(.2,.3,.1,-.3);

    for(FACE_ITERATOR<TV> it(example.grid,example.ghost);it.Valid();it.Next())
    {
        TV X=it.Location();
        TV V=M*X+A;
        if(rand_.Get_Uniform_Number(0,1)<inv_rate)
        {
            Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
            example.valid_xfer_data(it.face)=false;
        }
        else
        {
            example.velocity(it.face)=V(it.face.axis);
            example.valid_xfer_data(it.face)=true;
        }
    }

    driver.Write_Output_Files(1);
    driver.Grid_To_Particle();
    driver.Write_Output_Files(2);
    driver.Particle_To_Grid();
    driver.Write_Output_Files(3);

    for(int i=0;i<example.particles.number;i++)
    {
        TV X=example.particles.X(i);
        TV V=M*X+A;
//        LOG::printf("%P %P\n",example.particles.V(i)-V,example.particles.B(i)-M/(sqr(resolution)*3));
    }

    ARRAY<T> err;
    for(FACE_ITERATOR<TV> it(example.grid,-3);it.Valid();it.Next())
    {
        TV X=it.Location();
        TV V=M*X+A;
        example.velocity(it.face)-=V(it.face.axis);
        err.Append(example.velocity(it.face));
    }
    driver.Write_Output_Files(4);
    LOG::printf("%P\n",err);
    
    return 0;
}

