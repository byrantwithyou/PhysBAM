//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_SOLVER__
#define __FLUID_SOLVER__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class FLUID_STATE;
template<class TV> class FLUID_BC;

template<class TV>
class FLUID_SOLVER
{
public:
    typedef typename TV::SCALAR T;

    FLUID_SOLVER()=default;
    virtual ~FLUID_SOLVER()=default;

    virtual void Initialize()=0;
    virtual void Write(int frame) const=0;
    virtual void Read(int frame)=0;
    virtual T Compute_Dt(T time) const=0;
    virtual void Simulate_Time_Step(FLUID_BC<TV>* bc,T time,T dt)=0;
    virtual void Predict_Time_Step(T time,T dt)=0;
    virtual void Before_Time_Step(T time)=0;
    virtual void After_Time_Step(T time,T dt)=0;
    virtual void Before_Frame(int frame)=0;
    virtual void After_Frame(int frame)=0;
    virtual FLUID_STATE<TV>* Make_State() const=0;
    virtual FLUID_BC<TV>* Make_BC() const=0;

    virtual void Save(FLUID_STATE<TV>* fluid_state) const=0;
    virtual void Restore(const FLUID_STATE<TV>* fluid_state)=0;
    virtual T Diff_u(const FLUID_STATE<TV>* fluid_state) const=0;
    virtual T Diff_p(const FLUID_STATE<TV>* fluid_state) const=0;
};
}
#endif
