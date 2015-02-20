//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GATHER_SCATTER__
#define __GATHER_SCATTER__
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
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
    int threads;

    GATHER_SCATTER(const ARRAY<int>& simulated_particles);
    ~GATHER_SCATTER();

    void Set_Weights(const PARTICLE_GRID_WEIGHTS<TV>* weights_input);

    // Safe to write to grid data
    template<class FUNC>
    void Scatter(FUNC func,bool want_gradient)
    {
        Scatter([](int p,int thread_id){},func,[](int p,int thread_id){},want_gradient);
    }

    template<class PARTICLE_FUNC,class FUNC,class PARTICLE_FUNC2>
    void Scatter(PARTICLE_FUNC particle_func,FUNC func,PARTICLE_FUNC2 particle_func2,bool want_gradient)
    {
        for(int k=0;k<simulated_particles.m;k++){
            int p=simulated_particles(k);
            particle_func(p,0); // 0 = thread
            for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,0);it.Valid();it.Next())
                func(p,it,0);
            particle_func2(p,0);} // 0 = thread
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
        Scatter(particle_func,func,particle_func2,want_gradient); // Works because no threading
    }
    
//#####################################################################
};
}
#endif
