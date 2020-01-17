//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_SOLVER_PB__
#define __FLUID_SOLVER_PB__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include "FLUID_SOLVER.h"

namespace PhysBAM{

template<class TV> class SOLIDS_FLUIDS_DRIVER_UNIFORM;

template<class TV>
class FLUID_SOLVER_PB:public FLUID_SOLVER<TV>
{
public:
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::SCALAR T;

    SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>* driver = 0;

    mutable ARRAY<int,TV_INT> last_colors;

    FLUID_SOLVER_PB();
    virtual ~FLUID_SOLVER_PB();

    virtual void Initialize() override;
    virtual void Write(int frame) const override;
    virtual void Read(int frame) override;
    virtual T Compute_Dt(T time) const override;
    virtual void Simulate_Time_Step(FLUID_BOUNDARY_VECTOR<TV>* velocity,T time,T dt) override;
    virtual void Predict_Time_Step(T time,T dt) override;
    virtual void Before_Time_Step(T time) override;
    virtual void After_Time_Step(T time,T dt) override;
    virtual void Before_Frame(int frame) override;
    virtual void After_Frame(int frame) override;

    virtual FLUID_STATE<TV>* Make_State() const override;
    virtual FLUID_BC<TV>* Make_BC() const override;
    virtual FLUID_BOUNDARY_VECTOR<TV>* Make_Boundary_Vector() const override;
    virtual FLUID_REGIONS<TV>* Make_Regions() const override;

    virtual void Save(FLUID_STATE<TV>* fluid_state) const override;
    virtual void Restore(const FLUID_STATE<TV>* fluid_state) override;
    virtual T Diff_u(const FLUID_STATE<TV>* fluid_state) const override;
    virtual T Diff_p(const FLUID_STATE<TV>* fluid_state) const override;

    virtual void Get_Constraints(const SOLID_FLUID_INTERFACE<TV>* interface,
        ARRAY<FLUID_BOUNDARY_VECTOR<TV>*>& array,ARRAY<T>& rhs,
        FLUID_REGIONS<TV>* regions) const override;
    virtual void Compute_Region_Mapping(const FLUID_REGIONS<TV>* prev,
        const FLUID_REGIONS<TV>* next,ARRAY<int>& next_to_prev) const override;
    virtual void Get_Force(FLUID_BOUNDARY_VECTOR<TV>* force) const override;
};
}
#endif
