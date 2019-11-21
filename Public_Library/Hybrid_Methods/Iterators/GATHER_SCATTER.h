//#####################################################################
// Copyright 2015, Craig Schroeder, Andre Pradhana, and Gergely Klar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GATHER_SCATTER__
#define __GATHER_SCATTER__
#include <Core/Log/SCOPE.h>
#include <Core/Utilities/TIMER.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_FACE_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{

template<class TV> class MPM_EXAMPLE;
template<class TV> class MPM_OBJECTIVE;
template<class TV> class PARTICLE_GRID_WEIGHTS;

// TODO: probably need to leave an extra cell width to avoid collisions for MAC transfers.

template<class TV>
class GATHER_SCATTER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH T_SCRATCH;
    typedef typename PARTICLE_GRID_FACE_ITERATOR<TV>::SCRATCH T_FACE_SCRATCH;
    const ARRAY<int>& simulated_particles;
    const PARTICLE_GRID_WEIGHTS<TV>* weights;
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> face_weights;
    const GRID<TV>& grid;
    ARRAY<ARRAY<int> > bins;
    int threads;
    int partitions;

    GATHER_SCATTER(const GRID<TV>& grid,const ARRAY<int>& simulated_particles);
    ~GATHER_SCATTER();

    void Prepare_Scatter(const MPM_PARTICLES<TV>& particles);
    int Compute_Optimal_Bins(ARRAY<int>& bin_ends,ARRAY<int>& counts,ARRAY<int>& sum_counts,int min_width,int max_bins);

private:
    template<class DATA,class REDUCE> static enable_if_t<sizeof(((*(REDUCE*)0)(*(DATA*)0),1)),int> Test_Reduce(int);
    template<class DATA,class REDUCE> static char Test_Reduce(...);
public:
    
    // Safe to write to grid data
    template<class DATA,class REDUCE,class... Args>
    enable_if_t<sizeof(Test_Reduce<DATA,REDUCE>(0))==1>
    Scatter(bool want_gradient,REDUCE&& reduce,Args&&... args)
    {
        Scatter_Reduce<DATA>(want_gradient,[](const DATA&){},reduce,args...);
    }

    // Safe to write to grid data
    template<class DATA,class REDUCE,class... Args>
    enable_if_t<sizeof(Test_Reduce<DATA,REDUCE>(0))!=1>
    Scatter(bool want_gradient,REDUCE&& reduce,Args&&... args)
    {
        Scatter_Reduce<DATA>(want_gradient,reduce,args...);
    }

    // Safe to write to grid data
    template<class DATA,class... Args>
    void Scatter_Reduce(bool want_gradient,Args&&... args)
    {
        if(threads>=2){
#pragma omp parallel
            Scatter_Parallel<DATA>(want_gradient,args...);}
        else Scatter_Serial<DATA>(want_gradient,args...);
    }

    template<class DATA,class REDUCE,class... Args>
    void Scatter_Serial(bool want_gradient,REDUCE&& reduce,Args&&... args)
    {
        T_FACE_SCRATCH face_scratch;
        DATA data((DATA()));
        for(int k=0;k<simulated_particles.m;k++){
            int p=simulated_particles(k);
            Helper(face_scratch,want_gradient,p,data,args...);}
        if(reduce) reduce(data);
    }

    template<class DATA,class REDUCE,class... Args>
    void Scatter_Parallel(bool want_gradient,REDUCE&& reduce,Args&&... args)
    {
        for(int pass=0;pass<partitions;++pass){
#pragma omp for
            for(int i=0;i<threads;++i){
                DATA data((DATA()));
                int bin_id=i*partitions+pass;
                T_FACE_SCRATCH face_scratch;
                for(int k=0;k<bins(bin_id).m;k++){
                    int p=bins(bin_id)(k);
                    Helper(face_scratch,want_gradient,p,data,args...);}
#pragma omp critical
                reduce(data);
            }
#pragma omp barrier
        }
    }
    
    // Safe to write to grid data
    template<class DATA,class REDUCE,class... Args>
    enable_if_t<sizeof(Test_Reduce<DATA,REDUCE>(0))==1>
    Gather(bool want_gradient,REDUCE&& reduce,Args&&... args)
    {
        Gather_Reduce<DATA>(want_gradient,[](const DATA&){},reduce,args...);
    }

    // Safe to write to grid data
    template<class DATA,class REDUCE,class... Args>
    enable_if_t<sizeof(Test_Reduce<DATA,REDUCE>(0))!=1>
    Gather(bool want_gradient,REDUCE&& reduce,Args&&... args)
    {
        Gather_Reduce<DATA>(want_gradient,reduce,args...);
    }

    // Safe to write to particle data
    template<class DATA,class... Args>
    void Gather_Reduce(bool want_gradient,Args&&... args)
    {
#pragma omp parallel
        Gather_Parallel<DATA>(want_gradient,args...);
    }

    template<class DATA,class REDUCE,class... Args>
    void Gather_Parallel(bool want_gradient,REDUCE&& reduce,Args&&... args)
    {
#pragma omp for
        for(int tid=0;tid<threads;tid++){
            int a=tid*simulated_particles.m/threads;
            int b=(tid+1)*simulated_particles.m/threads;
            DATA data((DATA()));
            T_FACE_SCRATCH face_scratch;
            for(int k=a;k<b;k++){
                int p=simulated_particles(k);
                Helper(face_scratch,want_gradient,p,data,args...);}
#pragma omp critical
            reduce(data);
        }
    }

    template<class DATA,class F,class... Args>
    enable_if_t<sizeof(((*(F*)0)(0,*(DATA*)0),true))>
    Helper(T_FACE_SCRATCH& face_scratch,bool want_gradient,int p,DATA& data,F f,Args&&... args)
    {
        f(p,data);
        Helper(face_scratch,want_gradient,p,data,args...);
    }

    template<class DATA,class F,class... Args>
    enable_if_t<sizeof(((*(F*)0)(0,*(PARTICLE_GRID_ITERATOR<TV>*)0,*(DATA*)0),true))>
    Helper(T_FACE_SCRATCH& face_scratch,bool want_gradient,int p,DATA& data,F f,Args&&... args)
    {
        for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,face_scratch(0));it.Valid();it.Next())
            f(p,it,data);
        Helper(face_scratch,want_gradient,p,data,args...);
    }

    template<class DATA,class F,class... Args>
    enable_if_t<sizeof(((*(F*)0)(0,*(PARTICLE_GRID_FACE_ITERATOR<TV>*)0,*(DATA*)0),true))>
    Helper(T_FACE_SCRATCH& face_scratch,bool want_gradient,int p,DATA& data,F f,Args&&... args)
    {
        for(PARTICLE_GRID_FACE_ITERATOR<TV> it(face_weights,p,want_gradient,face_scratch);it.Valid();it.Next())
            f(p,it,data);
        Helper(face_scratch,want_gradient,p,data,args...);
    }

    template<class DATA>
    void Helper(T_FACE_SCRATCH& face_scratch,bool want_gradient,int p,DATA& data)
    {
    }
    
//#####################################################################
};
}
#endif
