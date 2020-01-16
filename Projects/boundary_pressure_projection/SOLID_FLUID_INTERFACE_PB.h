//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_FLUID_INTERFACE_PB__
#define __SOLID_FLUID_INTERFACE_PB__
#include <Core/Vectors/VECTOR.h>
#include "SOLID_FLUID_INTERFACE.h"
namespace PhysBAM{

template<class TV> class FLUID_SOLVER;
template<class TV> class SOLID_SOLVER;
template<class TV> class FLUID_BC;
template<class TV> class SOLID_BC;

template<class TV>
class SOLID_FLUID_INTERFACE_PB:public SOLID_FLUID_INTERFACE<TV>
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    ARRAY<bool> structure_is_thin;
    
    struct DEFORMABLE_ENTRY
    {
        int p;
        FACE_INDEX<TV::m> face;
        T w;
    };
    ARRAY<DEFORMABLE_ENTRY> deformable_weights;
    
    struct RIGID_ENTRY
    {
        int p;
        FACE_INDEX<TV::m> face;
        T w;
        typename TV::SPIN aw;
    };
    ARRAY<RIGID_ENTRY> rigid_weights;
    
    SOLID_FLUID_INTERFACE_PB()=default;
    virtual ~SOLID_FLUID_INTERFACE_PB()=default;

    virtual void Compute_BC(const FLUID_SOLVER<TV>* fluid_solver,SOLID_BC<TV>* solid_bc,T time,T dt) const override;
    virtual void Compute_BC(const SOLID_SOLVER<TV>* solid_solver,FLUID_BC<TV>* fluid_bc,T time,T dt) const override;

    virtual void Interpolate_Velocity(FLUID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v) override;
    virtual void Distribute_Force(SOLID_BOUNDARY_VECTOR<TV>* v, const FLUID_BOUNDARY_VECTOR<TV>* u) override;
    virtual void Get_Boundary(SOLID_BOUNDARY_VECTOR<TV>* v) override;
    virtual void Get_Boundary(FLUID_BOUNDARY_VECTOR<TV>* v) override;
    virtual void Compute_Coupling_Weights(const SOLID_SOLVER<TV>* solid_solver,const FLUID_SOLVER<TV>* fluid_solver,T time,T dt) override;
};
}
#endif
