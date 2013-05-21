//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __MPLE_DRIVER__
#define __MPLE_DRIVER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Fourier_Transforms/FFT_2D.h>
#include <PhysBAM_Tools/Fourier_Transforms/FFT_3D.h>
#include <PhysBAM_Tools/Fourier_Transforms/FFT_POLICY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include "MPLE_DOUBLE_WELL.h"
#include "MPLE_ITERATOR.h"
#include "MPLE_POINT.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

int omp_get_thread_num_helper()
{
    int tid;
#ifdef USE_OPENMP
    tid=omp_get_thread_num();
#else
    tid=0;
#endif
    return tid;
}

int omp_get_num_threads_helper()
{
    int threads;
#ifdef USE_OPENMP
#pragma omp parallel
#pragma omp master
    threads=omp_get_num_threads();
#else
    threads=1;
#endif
    return threads;
}

namespace PhysBAM{

template<class TV,int w>
class MPLE_DRIVER: public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MARCHING_CUBES<TV>::T_SURFACE T_SURFACE;
    typedef typename FFT_POLICY<TV>::FFT T_FFT;    

public:
    
    ARRAY<MPLE_POINT<TV,w> > points;     // data
    GRID<TV> grid;                       // grid

    ARRAY<T,TV_INT> u;                   // segmentation function
    ARRAY<T,TV_INT> force;               // force
    ARRAY<COMPLEX<T>,TV_INT> fft;        // fft representation
    ARRAY<T,TV_INT> source;              // source
    ARRAY<TV,TV_INT> location;           // flat index to location
    ARRAY<TV_INT,TV_INT> index;          // flat to vector index
    ARRAY<TV_INT,TV_INT> fft_index;      // flat to vector index

    T_FFT* transform;
    
    int array_m;
    int fft_array_m;
    
    T cfl;
    T spread;
    T rescale;
    T contour_value;
    T frame_dt;
    int frames;
    int timesteps;
    T mu,nu,epsilon;
    T dt,one_over_dx_squared;
    T cell_volume;
    T one_over_cell_volume;
    int threads;
    
    MPLE_DRIVER():transform(0),cfl((T)1),spread((T)1),rescale((T)1),contour_value((T).5),frame_dt((T)1/24),frames(100),mu(5e-4),nu(.05){}

    ~MPLE_DRIVER()
    {delete transform;}

    void Initialize()
    {
        RANGE<TV_INT> fft_range(RANGE<TV_INT>(TV_INT(),grid.Node_Indices().Edge_Lengths()));
        fft_range.max_corner(TV::m-1)/=2;
        fft_range.max_corner(TV::m-1)++;

        u.Resize(grid.Node_Indices(),false);
        fft.Resize(fft_range,false);
        fft_index.Resize(fft_range,false);
        force.Resize(grid.Node_Indices(),false);
        source.Resize(grid.Node_Indices(),false);
        location.Resize(grid.Node_Indices(),false);
        index.Resize(grid.Node_Indices(),false);
        array_m=location.array.m;

        transform=new T_FFT(grid);

        int k=0;
        for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next(),k++){
            location.array(k)=it.Location();
            index.array(k)=it.Node_Index();}

        k=0;
        for(RANGE_ITERATOR<TV_INT::m> it(fft_range);it.Valid();it.Next(),k++)
        {fft_index.array(k)=it.index;}
        fft_array_m=fft_index.array.m;

        for(int i=0;i<array_m;i++){
            u.array(i)=0;
            source.array(i)=0;}

#pragma omp parallel for        
        for(int i=0;i<points.m;i++)
            points(i).Update_Base_And_Weights(grid);

        for(int i=0;i<points.m;i++)
            for(MPLE_ITERATOR<TV,w> it(points(i));it.Valid();it.Next()){
                source(it.Node())+=it.Weight();
                u(it.Node())=1;}

        T max_value=0;
        for(int i=0;i<array_m;i++)
            max_value=max(max_value,source.array(i));

#pragma omp parallel for        
        for(int i=0;i<array_m;i++){
            source.array(i)*=rescale;
            source.array(i)=min(source.array(i),max_value);}

        dt=cfl*sqr(grid.dX(0))/(2*(TV::m+1));
        timesteps=(int)(frame_dt/dt);
        epsilon=spread*grid.dX(0);
        one_over_dx_squared=(T)1/sqr(grid.dX(0));
        cell_volume=grid.dX.Product();
        one_over_cell_volume=(T)1/cell_volume;
        PHYSBAM_ASSERT(grid.dX.Min()==grid.dX.Max());

        threads=omp_get_num_threads_helper();
        LOG::cout<<"Running on "<<threads<<" threads."<<std::endl;
    }        

    void Dump_Surface(SEGMENTED_CURVE_2D<T>& surface)
    {
        for(int i=0;i<surface.mesh.elements.m;i++){
            VECTOR<int,2>& s=surface.mesh.elements(i);
            Add_Debug_Object(VECTOR<VECTOR<T,2>,2>(surface.particles.X(s.x),surface.particles.X(s.y)),VECTOR<T,3>(0,.25,1));}
    }

    void Dump_Surface(TRIANGULATED_SURFACE<T>& surface)
    {
        for(int i=0;i<surface.mesh.elements.m;i++)
            Add_Debug_Object(VECTOR<TV,TV::m>(surface.particles.X.Subset(surface.mesh.elements(i))),TV(0,.25,1),TV(1,.25,0));
    }

    void Write(const char* title)
    {
        for(int i=0;i<points.m;i++)
            Add_Debug_Particle<TV>(points(i).X,VECTOR<T,3>(1,0,0));
        
        for(int i=0;i<array_m;i++){
            T value=u.array(i);
            if(value){
                Add_Debug_Particle(location.array(i),value>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(0,.5,1));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,value);}}
        
        T_SURFACE surface;
        MARCHING_CUBES<TV>::Create_Surface(surface,grid,u,contour_value);
        Dump_Surface(surface);

        Flush_Frame<TV>(title);
    }

    void Compute_Force()
    {
        // compute integral of u
        ARRAY<T> int_u_per_thread(threads);
#pragma omp parallel for
        for(int i=0;i<array_m;i++){
            int tid=omp_get_thread_num_helper();
            int_u_per_thread(tid)+=u.array(i);}
        T int_u=0;
        for(int i=0;i<threads;i++)
            int_u+=int_u_per_thread(i);
        int_u*=cell_volume;

        // compute integral of uw
        ARRAY<T> int_uw_per_thread(threads);
#pragma omp parallel for
        for(int i=0;i<points.m;i++){
            int tid=omp_get_thread_num_helper();
            for(MPLE_ITERATOR<TV,w> it(points(i));it.Valid();it.Next())
                int_uw_per_thread(tid)+=u(it.Node())*it.Weight();}
        T int_uw=0;
        for(int i=0;i<threads;i++)
            int_uw+=int_uw_per_thread(i);

        T c1=int_uw/int_u;

#pragma omp parallel for schedule(guided)
        for(int i=0;i<array_m;i++){
            force.array(i)=
                u.array(i)-
                (dt/epsilon)*MPLE_DOUBLE_WELL<T>::Gradient(u.array(i))+
                dt*mu*(u.array(i)*source.array(i)*one_over_cell_volume+(source.array(i)*one_over_cell_volume-c1));}
    }

    void Transform_Force()
    {
        transform->Transform(force,fft);
    }

    TV Index(const TV_INT& input)
    {
        TV result;
        for(int i=0;i<TV::m;i++){
            if(input(i)<grid.Counts()(i)/2+1) result(i)=input(i);
            else result(i)=input(i)-grid.Counts()(i);}
        return result;
    }

    void Compute_Transformed_U()
    {
#pragma omp parallel for
        for(int i=0;i<fft_array_m;i++){
            fft.array(i)*=(T)1/(1+dt*2*epsilon*(Index(fft_index.array(i))/grid.domain.Edge_Lengths()*(2*M_PI)).Magnitude_Squared());}
    }

    void Update_U()
    {transform->Inverse_Transform(fft,u);}

    void Clamp_U()
    {
#pragma omp parallel for
        for(int i=0;i<array_m;i++)
            u.array(i)=max(u.array(i),(T)0);
    }

    void Advance_Frame(int frame)
    {
        LOG::cout<<"Frame "<<frame<<std::endl;
        for(int i=0;i<timesteps;i++){
            Compute_Force();
            Transform_Force();
            Compute_Transformed_U();
            Update_U();
            Clamp_U();
            Write("substep");}
    }

    void Run()
    {
        Initialize();
        Write("Frame 0");
        for(int k=1;k<=frames;k++){
            Advance_Frame(k);
            char buff[100];
            sprintf(buff,"Frame %d",k);
            Write(buff);}
    }
};
}

#endif
