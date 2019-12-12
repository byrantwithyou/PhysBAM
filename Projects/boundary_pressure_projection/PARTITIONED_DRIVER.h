//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARTITIONED_DRIVER__
#define __PARTITIONED_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <climits>
namespace PhysBAM{

template<class TV> class FLUID_SOLVER;
template<class TV> class SOLID_SOLVER;
template<class TV> class FLUID_STATE;
template<class TV> class SOLID_STATE;
template<class TV> class SOLID_FLUID_INTERFACE;
template<class TV> class FLUID_BC;
template<class TV> class SOLID_BC;

template<class TV>
class PARTITIONED_DRIVER
{
public:
    typedef typename TV::SCALAR T;

    // User sets this
    FLUID_SOLVER<TV> * fluid_solver = 0;
    SOLID_SOLVER<TV> * solid_solver = 0;
    SOLID_FLUID_INTERFACE<TV> * interface = 0;

    int last_frame = 0;
    T frame_dt = 1;
    T max_dt = FLT_MAX;
    T min_dt = 0;
    T fixed_dt = 0;
    int max_subiterations = INT_MAX;

    T utol=0;
    T ptol=0;
    T xtol=0;
    T vtol=0;
    T p0tol=0;

private:
    T time;
    FLUID_STATE<TV> * fluid_new_state = 0;
    FLUID_STATE<TV> * fluid_old_state = 0;
    FLUID_STATE<TV> * fluid_prev_state = 0;
    SOLID_STATE<TV> * solid_new_state = 0;
    SOLID_STATE<TV> * solid_old_state = 0;
    SOLID_STATE<TV> * solid_prev_state = 0;
    FLUID_BC<TV> * fluid_bc = 0;
    SOLID_BC<TV> * solid_bc = 0;

public:
    PARTITIONED_DRIVER();
    ~PARTITIONED_DRIVER();
    void Run();

private:
    void Initialize();
    void Simulate_Frame(int frame);
    void Simulate_Time_Step(T dt);
    void Write(int frame) const;
    void Read(int frame);
    T Compute_Dt(T target_time,bool& done) const;
    bool Is_Subiteration_Converged() const;
};
}
#endif
