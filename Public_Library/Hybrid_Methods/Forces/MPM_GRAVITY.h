//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_GRAVITY
//#####################################################################
#ifndef __MPM_GRAVITY__
#define __MPM_GRAVITY__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class GATHER_SCATTER;
template<class T,int d> class ISOTROPIC_CONSTITUTIVE_MODEL;

template<class TV>
class MPM_GRAVITY:public PARTICLE_GRID_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef PARTICLE_GRID_FORCES<TV> BASE;
public:
    using BASE::particles;using BASE::force_helper;
    TV gravity;
    bool affect_all;
    GATHER_SCATTER<TV>& gather_scatter;

    MPM_GRAVITY(const MPM_FORCE_HELPER<TV>& force_helper,const TV& gravity,GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles);
    virtual ~MPM_GRAVITY();

//#####################################################################
    void Precompute(const T time,const T dt,bool want_dE,bool want_ddE) override;
    T Potential_Energy(const T time) const override;
    void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const override;
    void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const override;
//#####################################################################
};
}
#endif
