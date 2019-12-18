//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_SOLVER_PB__
#define __FLUID_SOLVER_PB__

#include "FLUID_SOLVER.h"

namespace PhysBAM{

template<class TV> class SOLIDS_FLUIDS_DRIVER_UNIFORM;

template<class TV>
class FLUID_SOLVER_PB:public FLUID_SOLVER<TV>
{
public:
    typedef typename TV::SCALAR T;

    SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>* driver = 0;
    
    FLUID_SOLVER_PB();
    virtual ~FLUID_SOLVER_PB();

    virtual void Initialize() override;
    virtual void Write(int frame) const override;
    virtual void Read(int frame) override;
    virtual T Compute_Dt(T time) const override;
    virtual void Simulate_Time_Step(FLUID_BC<TV>* bc,T time,T dt) override;
    virtual void Predict_Time_Step(T time,T dt) override;
    virtual void Before_Time_Step(T time) override;
    virtual void After_Time_Step(T time,T dt) override;
    virtual void Before_Frame(int frame) override;
    virtual void After_Frame(int frame) override;

    virtual FLUID_STATE<TV>* Make_State() const override;
    virtual FLUID_BC<TV>* Make_BC() const override;
    virtual FLUID_BOUNDARY_VECTOR<TV>* Make_Boundary_Vector() const override;

    virtual void Save(FLUID_STATE<TV>* fluid_state) const override;
    virtual void Restore(const FLUID_STATE<TV>* fluid_state) override;
    virtual T Diff_u(const FLUID_STATE<TV>* fluid_state) const override;
    virtual T Diff_p(const FLUID_STATE<TV>* fluid_state) const override;

    virtual void Get_Constraints(ARRAY<FLUID_BOUNDARY_VECTOR<TV>*>& array) const override;
};
}
#endif
