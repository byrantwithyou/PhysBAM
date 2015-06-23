//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_TENSION_BOUNDARY_CONDITION
//#####################################################################
#ifndef __SURFACE_TENSION_BOUNDARY_CONDITION__
#define __SURFACE_TENSION_BOUNDARY_CONDITION__
#include <Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{
template<class TV>
class SURFACE_TENSION_BOUNDARY_CONDITION:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
    const LEVELSET<TV>& levelset;
public:
    T surface_tension_coefficient;

    SURFACE_TENSION_BOUNDARY_CONDITION(const LEVELSET<TV>& levelset,T surface_tension_boundary_condition_input);
    virtual ~SURFACE_TENSION_BOUNDARY_CONDITION();

    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time) override;
};
}
#endif
