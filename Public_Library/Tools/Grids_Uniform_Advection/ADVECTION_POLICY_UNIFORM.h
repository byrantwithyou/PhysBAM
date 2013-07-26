//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADVECTION_POLICY_UNIFORM__
#define __ADVECTION_POLICY_UNIFORM__

#include <Tools/Advection/ADVECTION_FORWARD.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_UNIFORM_FORWARD.h>

namespace PhysBAM{

template<class TV>
struct ADVECTION_POLICY
{
private:
    typedef typename TV::SCALAR T;
public:
    // normal
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T> ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ::PhysBAM::ADVECTION_CONSERVATIVE_ENO<TV,T> ADVECTION_CONSERVATIVE_ENO;
    typedef ADVECTION_HAMILTON_JACOBI_WENO<TV,T> ADVECTION_HAMILTON_JACOBI_WENO_SCALAR;
};

}
#endif
