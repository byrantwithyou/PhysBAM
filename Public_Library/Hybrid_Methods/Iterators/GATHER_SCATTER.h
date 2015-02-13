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

    GATHER_SCATTER(const ARRAY<int>& simulated_particles,const PARTICLE_GRID_WEIGHTS<TV>* weights)
        :simulated_particles(simulated_particles),weights(weights)
    {}

    // Safe to write to grid data
    template<class FUNC>
    void Scatter(FUNC func,bool want_gradient)
    {
        Scatter([](int p,int thread_id){},func,want_gradient);
    }

    template<class PARTICLE_FUNC,class FUNC>
    void Scatter(PARTICLE_FUNC particle_func,FUNC func,bool want_gradient)
    {
        for(int k=0;k<simulated_particles.m;k++){
            int p=simulated_particles(k);
            particle_func(p,0); // 0 = thread
            for(PARTICLE_GRID_ITERATOR<TV> it(weights,p,want_gradient,0);it.Valid();it.Next())
                func(p,it,0);}
    }

    // Safe to write to particle data
    template<class FUNC>
    void Gather(FUNC func)
    {
        Gather([](int p,int thread_id){},func);
    }

    template<class PARTICLE_FUNC,class FUNC>
    void Gather(PARTICLE_FUNC particle_func,FUNC func)
    {
        Scatter(particle_func,func); // Works because no threading
    }
    
//#####################################################################
};
}
#endif
