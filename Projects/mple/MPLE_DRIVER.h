//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __MPLE_DRIVER__
#define __MPLE_DRIVER__

#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Fourier_Transforms/FFT_2D.h>
#include <Grid_Tools/Fourier_Transforms/FFT_3D.h>
#include <Grid_Tools/Fourier_Transforms/FFT_POLICY.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
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
class MPLE_DRIVER:public NONCOPYABLE
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
    ARRAY<std::complex<T>,TV_INT> fft;        // fft representation
    ARRAY<T,TV_INT> source;              // source
    ARRAY<TV,TV_INT> location;           // flat index to location
    ARRAY<TV_INT,TV_INT> index;          // flat to vector index
    ARRAY<TV_INT,TV_INT> fft_index;      // flat to vector index

    T_FFT* transform;
    
    int array_m;
    int fft_array_m;
    
    std::string output_directory;
    T cfl;
    T spread;
    T rescale;
    T identity;
    T contour_value;
    int frames;
    T mu,epsilon;
    T dt,one_over_dx_squared;
    T cell_volume;
    T one_over_cell_volume;
    int threads;
    
    MPLE_DRIVER():transform(0),output_directory("output"),cfl((T)1),spread((T)1),rescale((T)1),identity((T)0),contour_value((T).5),frames(100),mu(5e-4){}

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
        epsilon=spread*grid.dX(0);
        one_over_dx_squared=(T)1/sqr(grid.dX(0));
        cell_volume=grid.dX.Product();
        one_over_cell_volume=(T)1/cell_volume;
        PHYSBAM_ASSERT(grid.dX.Min()==grid.dX.Max());

        threads=omp_get_num_threads_helper();
        LOG::cout<<"Running on "<<threads<<" threads."<<std::endl;
    }        

    void Dump_Surface(SEGMENTED_CURVE_2D<T>& surface,std::ofstream& out)
    {
        for(int i=0;i<surface.mesh.elements.m;i++){
            VECTOR<int,2>& s=surface.mesh.elements(i);
            Add_Debug_Object(VECTOR<VECTOR<T,2>,2>(surface.particles.X(s.x),surface.particles.X(s.y)),VECTOR<T,3>(0,.25,1));}
    }

    void Emit_Vector(std::ofstream& out,const VECTOR<T,3>& v, const char* str="")
    {
        out<<"<"<<v(0);
        for(int i=1;i<3;i++) out<<","<<v(i);
        out<<">"<<str;
    }

    void Dump_Surface(TRIANGULATED_SURFACE<T>& surface,std::ofstream& out)
    {
        surface.Update_Vertex_Normals();
        for(int i=0;i<surface.mesh.elements.m;i++){
            out<<"smooth_triangle{";
            VECTOR<TV,3> X(surface.particles.X.Subset(surface.mesh.elements(i)));
            VECTOR<TV,3> N;
            Add_Debug_Object(X,TV(0,.25,1),TV(1,.25,0));
            if(surface.face_vertex_normals) N=(*surface.face_vertex_normals)(i);
            else N=VECTOR<TV,3>(surface.vertex_normals->Subset(surface.mesh.elements(i)));
            Emit_Vector(out,X(0),",");
            Emit_Vector(out,N(0),",");
            Emit_Vector(out,X(1),",");
            Emit_Vector(out,N(1),",");
            Emit_Vector(out,X(2),",");
            Emit_Vector(out,N(2),"");
            out<<"}"<<std::endl;;}
    }

    void Write(const int frame)
    {
        char buff[100];
        sprintf(buff,"Frame %d",frame);

        FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(LOG::sprintf("%s/%d",output_directory,frame));

        for(int i=0;i<points.m;i++)
            Add_Debug_Particle<TV>(points(i).X,VECTOR<T,3>(1,0,0));
        
        for(int i=0;i<array_m;i++){
            T value=u.array(i);
            if(value){
                Add_Debug_Particle(location.array(i),value>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(0,.5,1));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,value);}}
        
        T_SURFACE surface;
        MARCHING_CUBES<TV>::Create_Surface(surface,grid,u,contour_value);
        std::ofstream out(LOG::sprintf("%s/%d/surface",output_directory,frame).c_str());
        Dump_Surface(surface,out);
        out.close();

        Flush_Frame<TV>(buff);
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
                (dt*identity+1)*u.array(i)-
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
            if(input(i)<grid.counts(i)/2+1) result(i)=input(i);
            else result(i)=input(i)-grid.counts(i);}
        return result;
    }

    void Compute_Transformed_U()
    {
#pragma omp parallel for
        for(int i=0;i<fft_array_m;i++){
            fft.array(i)*=(T)1/(dt*identity+1+dt*2*epsilon*(Index(fft_index.array(i))/grid.domain.Edge_Lengths()*(2*M_PI)).Magnitude_Squared());}
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
        Compute_Force();
        Transform_Force();
        Compute_Transformed_U();
        Update_U();
        Clamp_U();
    }

    void Dump_Points()
    {
        FILE_UTILITIES::Create_Directory(output_directory);
        std::ofstream out(LOG::sprintf("%s/particles",output_directory).c_str());
        for(int i=0;i<points.m;i++){
            const MPLE_POINT<TV,w>& point=points(i);
            out<<"sphere{<";
            for(int j=0;j<TV::m;j++){
                out<<point.X(j);
                if(j<TV::m-1) out<<",";}
            out<<">sphere_radius}"<<std::endl;
        }
        out.close();
    }
    
    void Run()
    {
        Initialize();
        Dump_Points();
        Write(0);
        for(int k=1;k<=frames;k++){
            Advance_Frame(k);
            Write(k);}
    }
};
}

#endif
