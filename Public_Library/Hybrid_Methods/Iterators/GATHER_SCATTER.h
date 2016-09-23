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
    template<class DATA,class... Args>
    void Scatter(bool want_gradient,Args&&... args)
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
                        Helper_P(scratch,want_gradient,p,data,args...);}}}}
        else{
            T_SCRATCH scratch;
            DATA data((DATA()));
            for(int k=0;k<simulated_particles.m;k++){
                int p=simulated_particles(k);
                Helper_P(scratch,want_gradient,p,data,args...);}}
    }

    // Safe to write to particle data
    template<class DATA,class... Args>
    void Gather(bool want_gradient,Args&&... args)
    {
#pragma omp parallel for
        for(int tid=0;tid<threads;tid++){
            int a=tid*simulated_particles.m/threads;
            int b=(tid+1)*simulated_particles.m/threads;
            DATA data((DATA()));
            T_SCRATCH scratch;
            for(int k=a;k<b;k++){
                int p=simulated_particles(k);
                Helper_P(scratch,want_gradient,p,data,args...);}}
    }

    template<class DATA,class F,class... Args>
    void Helper_P(T_SCRATCH& scratch,bool want_gradient,int p,DATA& data,F f,Args&&... args)
    {
        f(p,data);
        Helper_G(scratch,want_gradient,p,data,args...);
    }

    template<class DATA,class... Args>
    void Helper_P(T_SCRATCH& scratch,bool want_gradient,int p,DATA& data,int f,Args&&... args)
    {
        Helper_G(scratch,want_gradient,p,data,args...);
    }

    template<class DATA>
    void Helper_P(T_SCRATCH& scratch,bool want_gradient,int p,DATA& data)
    {
    }

    template<class DATA,class F,class... Args>
    void Helper_G(T_SCRATCH& scratch,bool want_gradient,int p,DATA& data,F f,Args&&... args)
    {
        for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,scratch);it.Valid();it.Next())
            f(p,it,data);
        Helper_P(scratch,want_gradient,p,data,args...);
    }

    template<class DATA,class... Args>
    void Helper_G(T_SCRATCH& scratch,bool want_gradient,int p,DATA& data,int f,Args&&... args)
    {
        Helper_P(scratch,want_gradient,p,data,args...);
    }

    template<class DATA>
    void Helper_G(T_SCRATCH& scratch,bool want_gradient,int p,DATA& data)
    {
    }
    
//#####################################################################
};
}
#endif
