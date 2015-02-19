//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_FINITE_ELEMENTS
//#####################################################################
#ifndef __MPM_FINITE_ELEMENTS__
#define __MPM_FINITE_ELEMENTS__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class GATHER_SCATTER;
template<class T,int d> class ISOTROPIC_CONSTITUTIVE_MODEL;

template<class TV>
class MPM_FINITE_ELEMENTS:public PARTICLE_GRID_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef PARTICLE_GRID_FORCES<TV> BASE;
public:
    using BASE::particles;
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> > sigma;
    ARRAY<MATRIX<T,TV::m> > U;
    bool affect_all;
    GATHER_SCATTER<TV>& gather_scatter;
    mutable ARRAY<MATRIX<T,TV::m> > tmp;
    ARRAY<DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,TV::m> > dPi_dF;

    MPM_FINITE_ELEMENTS(MPM_PARTICLES<TV>& particles,ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model,
        GATHER_SCATTER<TV>& gather_scatter_input,bool affect_all);
    virtual ~MPM_FINITE_ELEMENTS();

//#####################################################################
    virtual void Precompute(const T time)=0;
    virtual T Potential_Energy(const T time) const=0;
    virtual void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const=0;
    virtual void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const=0;
//#####################################################################
};
}
#endif
