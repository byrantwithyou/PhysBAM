//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_SOLVER__
#define __SOLID_SOLVER__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class SOLID_STATE;
template<class TV> class SOLID_BC;
template<class TV> class SOLID_BOUNDARY_VECTOR;

template<class TV>
class SOLID_SOLVER
{
public:
    typedef typename TV::SCALAR T;

    SOLID_SOLVER()=default;
    virtual ~SOLID_SOLVER()=default;

    virtual void Initialize()=0;
    virtual void Write(int frame) const=0;
    virtual void Read(int frame)=0;
    virtual T Compute_Dt(T time) const=0;
    virtual void Simulate_Time_Step(SOLID_BC<TV>* bc,T time,T dt)=0;
    virtual void Before_Time_Step(T time)=0;
    virtual void After_Time_Step(T time,T dt)=0;
    virtual void Before_Frame(int frame)=0;
    virtual void After_Frame(int frame)=0;
    virtual SOLID_STATE<TV>* Make_State() const=0;
    virtual SOLID_BC<TV>* Make_BC() const=0;
    virtual SOLID_BOUNDARY_VECTOR<TV>* Make_Boundary_Vector() const=0;

    virtual void Save(SOLID_STATE<TV>* solid_state) const=0;
    virtual void Restore(const SOLID_STATE<TV>* solid_state)=0;
    virtual T Diff_x(const SOLID_STATE<TV>* solid_state) const=0;
    virtual T Diff_v(const SOLID_STATE<TV>* solid_state) const=0;

    virtual void Fill_Boundary_Vector(SOLID_BOUNDARY_VECTOR<TV>* v) const=0;
    virtual void Apply_Velocity_Change(T c,const SOLID_BOUNDARY_VECTOR<TV>* v) const=0;
    virtual T Inner_Product(const SOLID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v)=0;
    virtual void Mass_Inverse(SOLID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v)=0;
};
}
#endif
