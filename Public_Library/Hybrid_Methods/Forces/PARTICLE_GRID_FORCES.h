//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_GRID_FORCES
//#####################################################################
#ifndef __PARTICLE_GRID_FORCES__
#define __PARTICLE_GRID_FORCES__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class MPM_FORCE_HELPER;

template<class TV>
class PARTICLE_GRID_FORCES
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    const MPM_FORCE_HELPER<TV>& force_helper;
    const MPM_PARTICLES<TV>& particles;

    PARTICLE_GRID_FORCES(const MPM_FORCE_HELPER<TV>& force_helper);
    PARTICLE_GRID_FORCES(const PARTICLE_GRID_FORCES&) = delete;
    void operator=(const PARTICLE_GRID_FORCES&) = delete;
    virtual ~PARTICLE_GRID_FORCES();

//#####################################################################
    virtual void Precompute(const T time,const T dt,bool want_dE,bool want_ddE)=0;
    virtual T Potential_Energy(const T time) const=0;
    virtual void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const=0;
    virtual void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const=0;
//#####################################################################
};
}
#endif
