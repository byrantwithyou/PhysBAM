//#####################################################################
// Copyright 2015, Craig Schroeder, Andre Pradhana, and Gergely Klar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GATHER_SCATTER__
#define __GATHER_SCATTER__
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Utilities/TIMER.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{

template<class TV> class MPM_EXAMPLE;
template<class TV> class MPM_OBJECTIVE;
template<class TV> class PARTICLE_GRID_WEIGHTS;

template<class TV>
class GATHER_SCATTER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH T_SCRATCH;
    const ARRAY<int>& simulated_particles;
    const PARTICLE_GRID_WEIGHTS<TV>* weights;
    const GRID<TV>& grid;
    ARRAY<ARRAY<int> > bins;
    int threads;
    int partitions;

    GATHER_SCATTER(const GRID<TV>& grid,const ARRAY<int>& simulated_particles);
    ~GATHER_SCATTER();

    void Prepare_Scatter(const MPM_PARTICLES<TV>& particles);
    int Compute_Optimal_Bins(ARRAY<int>& bin_ends,ARRAY<int>& counts,ARRAY<int>& sum_counts,int min_width,int max_bins);

    // Safe to write to grid data
    template<class DATA,class FUNC>
    void Scatter(FUNC func,bool want_gradient)
    {
        Scatter<DATA>([](int p,DATA data){},func,[](int p,DATA data){},want_gradient);
    }

    template<class DATA,class PARTICLE_FUNC,class FUNC,class PARTICLE_FUNC2>
    void Scatter(PARTICLE_FUNC particle_func,FUNC func,PARTICLE_FUNC2 particle_func2,bool want_gradient)
    {
        if(threads>=2){
            for(int pass=0;pass<partitions;++pass){
#pragma omp parallel for
                for(int i=0;i<threads;++i){
                    DATA data((DATA()));
                    int bin_id=i*partitions+pass;
                    T_SCRATCH scratch;
                    for(int k=0;k<bins(bin_id).m;k++){
                        int p=bins(bin_id)(k);
                        particle_func(p,data);
                        for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,scratch);it.Valid();it.Next())
                            func(p,it,data);
                        particle_func2(p,data);}}}}
        else{
            T_SCRATCH scratch;
            DATA data((DATA()));
            for(int k=0;k<simulated_particles.m;k++){
                int p=simulated_particles(k);
                particle_func(p,data);
                for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,scratch);it.Valid();it.Next())
                    func(p,it,data);
                particle_func2(p,data);}}
    }

    // Safe to write to particle data
    template<class DATA,class FUNC>
    void Gather(FUNC func,bool want_gradient)
    {
        Gather<DATA>([](int p,int thread_id){},func,[](int p,int thread_id){},want_gradient);
    }

    template<class DATA,class PARTICLE_FUNC,class FUNC,class PARTICLE_FUNC2>
    void Gather(PARTICLE_FUNC particle_func,FUNC func,PARTICLE_FUNC2 particle_func2,bool want_gradient)
    {
#pragma omp parallel for
        for(int tid=0;tid<threads;tid++){
            int a=tid*simulated_particles.m/threads;
            int b=(tid+1)*simulated_particles.m/threads;
            DATA data((DATA()));
            T_SCRATCH scratch;
            for(int k=a;k<b;k++){
                int p=simulated_particles(k);
                particle_func(p,data);
                for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,scratch);it.Valid();it.Next())
                    func(p,it,data);
                particle_func2(p,data);}}
    }
    
//#####################################################################
};
}
#endif
