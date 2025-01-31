//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_BOUNDARY_CONDITION_SOURCES
//#####################################################################
#ifndef __IMPLICIT_BOUNDARY_CONDITION_SOURCES__
#define __IMPLICIT_BOUNDARY_CONDITION_SOURCES__
#include <Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{
template<class TV> class FLUIDS_PARAMETERS_CALLBACKS;

template<class TV>
class IMPLICIT_BOUNDARY_CONDITION_SOURCES:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
public:
    FLUIDS_PARAMETERS_CALLBACKS<TV>& callbacks;

    IMPLICIT_BOUNDARY_CONDITION_SOURCES(FLUIDS_PARAMETERS_CALLBACKS<TV>& callbacks_input);
    virtual ~IMPLICIT_BOUNDARY_CONDITION_SOURCES();

    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time) override;
};
}
#endif
