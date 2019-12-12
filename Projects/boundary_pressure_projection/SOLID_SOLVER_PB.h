//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_SOLVER_PB__
#define __SOLID_SOLVER_PB__

#include "SOLID_SOLVER.h"

namespace PhysBAM{

template<class TV> class SOLIDS_FLUIDS_DRIVER_UNIFORM;

template<class TV>
class SOLID_SOLVER_PB:public SOLID_SOLVER<TV>
{
public:
    typedef typename TV::SCALAR T;

    SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>* driver = 0;
    
    SOLID_SOLVER_PB();
    virtual ~SOLID_SOLVER_PB();

    virtual void Initialize() override;
    virtual void Write(int frame) const override;
    virtual void Read(int frame) override;
    virtual T Compute_Dt(T time) const override;
    virtual void Simulate_Time_Step(SOLID_BC<TV>* bc,T time,T dt) override;
    virtual void Before_Time_Step(T time) override;
    virtual void After_Time_Step(T time,T dt) override;
    virtual void Before_Frame(int frame) override;
    virtual void After_Frame(int frame) override;
    virtual SOLID_STATE<TV>* Make_State() const override;
    virtual SOLID_BC<TV>* Make_BC() const override;

    virtual void Save(SOLID_STATE<TV>* solid_state) const override;
    virtual void Restore(const SOLID_STATE<TV>* solid_state) override;
    virtual T Diff_x(const SOLID_STATE<TV>* solid_state) const override;
    virtual T Diff_v(const SOLID_STATE<TV>* solid_state) const override;
};
}
#endif
