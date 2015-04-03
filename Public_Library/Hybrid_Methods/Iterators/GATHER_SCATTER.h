//#####################################################################
// Copyright 2015, Craig Schroeder, Andre Pradhana, and Gergely Klar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GATHER_SCATTER__
#define __GATHER_SCATTER__
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Utilities/TIMER.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <omp.h>
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
    const ARRAY<int>& simulated_particles;
    const PARTICLE_GRID_WEIGHTS<TV>* weights;
    ARRAY<ARRAY<int> > bins;
    int threads;
    int scatter_threads;
    int min_width;
    int partitions;
    int splitting_direction;

    GATHER_SCATTER(const ARRAY<int>& simulated_particles);
    ~GATHER_SCATTER();

    void Initialize(const GRID<TV>& grid)
    {
        splitting_direction=0;
        partitions=weights->Order()<3?2:3;
        min_width=std::max(2,(weights->Order()+1)/2);
        scatter_threads=std::min(threads,(grid.numbers_of_cells(splitting_direction)+min_width*partitions-1)/(min_width*partitions));
    }

    void Set_Weights(const PARTICLE_GRID_WEIGHTS<TV>* weights_input);

    void Prepare_Scatter(const MPM_PARTICLES<TV>& particles,const GRID<TV>& grid)
    {
        if(scatter_threads<2) return;
        ARRAY<ARRAY<int> > precount(grid.numbers_of_cells(splitting_direction));
        ARRAY<omp_lock_t> locks(grid.numbers_of_cells(splitting_direction));
        for(int i=0;i<grid.numbers_of_cells(splitting_direction);++i)
            omp_init_lock(&locks(i));
#pragma omp parallel for
        for(int k=0;k<simulated_particles.m;++k){
            int p=simulated_particles(k);
            TV_INT cell=grid.Cell(particles.X(p),0);
            omp_set_lock(&locks(cell(splitting_direction)));
            precount(cell(splitting_direction)).Append(p);
            omp_unset_lock(&locks(cell(splitting_direction)));}
        int ideal_num_part_per_bin=simulated_particles.m/(scatter_threads*partitions);
        int count=0,width=0;
        bins.Remove_All();
        int bin=0;
        bins.Append(ARRAY<int>());
        for (int i=0;i<precount.m;++i){
            count+=precount(i).m;
            for (int j=0;j<precount(i).m;j++)
                bins(bin).Append(precount(i)(j)); 
            if(++width>=min_width && count>=ideal_num_part_per_bin && i!=precount.m-1){
                bins.Append(ARRAY<int>());
                ++bin;
                width=0;
                count=0;}}
        if((bins.end()-1)->m==0) bins.Remove_End();
        if(bins.m/partitions<scatter_threads) scatter_threads=(bins.m+partitions-1)/partitions;
    }

    // Safe to write to grid data
    template<class FUNC,class DATA=int>
    void Scatter(FUNC func,bool want_gradient)
    {
        Scatter([](int p,DATA data){},func,[](int p,DATA data){},want_gradient);
    }

    template<class PARTICLE_FUNC,class FUNC,class PARTICLE_FUNC2,class DATA=int>
    void Scatter(PARTICLE_FUNC particle_func,FUNC func,PARTICLE_FUNC2 particle_func2,bool want_gradient,DATA data_=DATA())
    {
#ifdef USE_OPENMP
        if(scatter_threads<2){
            DATA data(data_);
            for(int k=0;k<simulated_particles.m;k++){
                int p=simulated_particles(k);
                particle_func(p,data);
                for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,0);it.Valid();it.Next())
                    func(p,it,data);
                particle_func2(p,data);}
        }else{
            omp_set_num_threads(scatter_threads);
            for(int pass=0;pass<partitions;++pass){
#               pragma omp parallel for schedule(static)
                for (int i=0;i<scatter_threads;++i){
                    int tid=omp_get_thread_num();
                    DATA data(data_);
                    int bin_id=tid*partitions+pass;
                    typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH scratch;
                    if (bin_id<bins.m){
                        for(int k=0;k<bins(bin_id).m;k++){
                            int p=bins(bin_id)(k);
                            particle_func(p,data);
                            for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,scratch);it.Valid();it.Next())
                                func(p,it,data);
                            particle_func2(p,data);}}}}}
#else
        DATA data(data_);
        for(int k=0;k<simulated_particles.m;k++){
            int p=simulated_particles(k);
            particle_func(p,data);
            for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,0);it.Valid();it.Next())
                func(p,it,data);
        particle_func2(p,data);}
#endif
    }

    // Safe to write to particle data
    template<class FUNC>
    void Gather(FUNC func,bool want_gradient)
    {
        Gather([](int p,int thread_id){},func,[](int p,int thread_id){},want_gradient);
    }

    template<class PARTICLE_FUNC,class FUNC,class PARTICLE_FUNC2>
    void Gather(PARTICLE_FUNC particle_func,FUNC func,PARTICLE_FUNC2 particle_func2,bool want_gradient)
    {
#ifdef USE_OPENMP
        omp_set_num_threads(threads);
#pragma omp parallel for
        for(int k=0;k<simulated_particles.m;k++){
            int tid=omp_get_thread_num();
            int p=simulated_particles(k);
            particle_func(p,tid);
            for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,tid);it.Valid();it.Next())
                func(p,it,tid);
            particle_func2(p,tid);}
#else
        Scatter(particle_func,func,particle_func2,want_gradient); // Works because no threading
#endif
    }
    
//#####################################################################
};
}
#endif
