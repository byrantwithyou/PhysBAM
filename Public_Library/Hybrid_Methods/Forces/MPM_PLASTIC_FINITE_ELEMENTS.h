//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_PLASTIC_FINITE_ELEMENTS
//#####################################################################
#ifndef __MPM_PLASTIC_FINITE_ELEMENTS__
#define __MPM_PLASTIC_FINITE_ELEMENTS__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class GATHER_SCATTER;
template<class T,int d> class ISOTROPIC_CONSTITUTIVE_MODEL;
template<class TV> class MPM_PLASTICITY_MODEL;

template<class TV>
class MPM_PLASTIC_FINITE_ELEMENTS:public PARTICLE_GRID_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef PARTICLE_GRID_FORCES<TV> BASE;
public:
    using BASE::particles;using BASE::force_helper;
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model;
    MPM_PLASTICITY_MODEL<TV>& plasticity;
    ARRAY<MATRIX<T,TV::m> > U,FV,PFT;
    bool affect_all;
    GATHER_SCATTER<TV>& gather_scatter;
    mutable ARRAY<MATRIX<T,TV::m> > tmp;

    struct HESSIAN_HELPER
    {
        MATRIX<T,TV::m> diag;
        typename TV::SPIN off_diag_sum,off_diag_diff;

        MATRIX<T,TV::m> Times(const MATRIX<T,TV::m>& M) const;
        void Set(const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dPi_dF);
    };

    ARRAY<HESSIAN_HELPER> hessian_helper;

    MPM_PLASTIC_FINITE_ELEMENTS(const MPM_FORCE_HELPER<TV>& force_helper,
        ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model,
        GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles,
        MPM_PLASTICITY_MODEL<TV>& plasticity);
    virtual ~MPM_PLASTIC_FINITE_ELEMENTS();

//#####################################################################
    void Precompute(const T time,const T dt,bool want_dE,bool want_ddE) override;
    T Potential_Energy(const T time) const override;
    void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const override;
    void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const override;
//#####################################################################
};
}
#endif
