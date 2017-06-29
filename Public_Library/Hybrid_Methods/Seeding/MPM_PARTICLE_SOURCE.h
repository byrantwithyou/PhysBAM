//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace MPM_PARTICLE_SOURCE
//#####################################################################
#ifndef __MPM_PARTICLE_SOURCE__
#define __MPM_PARTICLE_SOURCE__

#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/MATRIX.h>
#include <functional>

namespace PhysBAM{
template<class TV> class IMPLICIT_OBJECT;
template<class TV> class POISSON_DISK;
template<class T> class RANDOM_NUMBERS;
template<class T,int m,int n> class MATRIX;

// State of a particle at time t that was seeded at location X at time ts.
// xp = position
// vp = dxp/dt
// xp_ts,vp_ts = diff on ts
// xp_x,vp_x = diff on X
template<class TV>
struct SOURCE_PATH
{
    TV xp,vp,xp_ts,vp_ts;
    MATRIX<typename TV::SCALAR,TV::m> xp_x,vp_x;
};

template<class TV>
class MPM_PARTICLE_SOURCE
{
public:
    typedef typename TV::SCALAR T;
    
    TV X0,n; // Seed on plane containing point X0, normal direction n.
    IMPLICIT_OBJECT<TV>* io; // Only this part of the plane
    POISSON_DISK<TV>& poisson_disk; // Use this for sampling
    RANDOM_NUMBERS<T>& random;

    std::function<void(TV X,T ts,T t,SOURCE_PATH<TV>& p)> path;

private:
    ARRAY<TV> seed_X;
    FRAME<TV> frame; // Sample space to real space
    RANGE<TV> sample_box;
public:

    MPM_PARTICLE_SOURCE(POISSON_DISK<TV>& poisson_disk,RANDOM_NUMBERS<T>& random,
        const TV& X0,const TV& n,IMPLICIT_OBJECT<TV>* io,
        std::function<void(TV X,T ts,T t,SOURCE_PATH<TV>& p)> path);
    ~MPM_PARTICLE_SOURCE()=default;

    // Seed particles for the next dt time.  Return particle positions,
    // velocities, and (if desired) dv/dx.
    void Seed(T t0,T dt,ARRAY<TV>& X,ARRAY<TV>& V,ARRAY<MATRIX<T,TV::m> >* dV);

    // World space locations of seed points; try to fit new points to this.
    void Seed_Points(ARRAY_VIEW<const TV> X);
};
//#####################################################################
}
#endif
