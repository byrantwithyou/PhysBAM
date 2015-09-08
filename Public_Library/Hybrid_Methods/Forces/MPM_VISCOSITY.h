//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_VISCOSITY
//#####################################################################
#ifndef __MPM_VISCOSITY__
#define __MPM_VISCOSITY__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Hybrid_Methods/Forces/OLDROYD_CONSTITUTIVE_MODEL.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class GATHER_SCATTER;
template<class TV> class MPM_FORCE_HELPER;

template<class TV>
class MPM_VISCOSITY:public PARTICLE_GRID_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef PARTICLE_GRID_FORCES<TV> BASE;
public:
    using BASE::particles;using BASE::force_helper;
    bool affect_all;
    GATHER_SCATTER<TV>& gather_scatter;
    mutable ARRAY<MATRIX<T,TV::m> > tmp;
    T stored_dt;
    T viscosity;

    MPM_VISCOSITY(MPM_FORCE_HELPER<TV>& force_helper,GATHER_SCATTER<TV>& gather_scatter_input,
        ARRAY<int>* affected_particles,const T viscosity);
    virtual ~MPM_VISCOSITY();

//#####################################################################
    void Precompute(const T time,const T dt) override;
    T Potential_Energy(const T time) const override;
    void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const override;
    void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const override;
//#####################################################################
};
}
#endif