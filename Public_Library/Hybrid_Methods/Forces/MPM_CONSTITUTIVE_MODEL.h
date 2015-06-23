//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_CONSTITUTIVE_MODEL
//#####################################################################
#ifndef __MPM_CONSTITUTIVE_MODEL__
#define __MPM_CONSTITUTIVE_MODEL__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;

template<class TV>
class MPM_CONSTITUTIVE_MODEL:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_CONSTITUTIVE_MODEL(MPM_PARTICLES<TV>& particles);
    virtual ~MPM_CONSTITUTIVE_MODEL();

//#####################################################################
    virtual T Psi(const T time)=0;
    virtual T Potential_Energy(const T time) const=0;
    virtual void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const=0;
    virtual void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const=0;
//#####################################################################
};
}
#endif
