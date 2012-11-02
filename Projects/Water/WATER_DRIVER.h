//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WATER_DRIVER__
#define __WATER_DRIVER__
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/KINEMATIC_EVOLUTION.h>
namespace PhysBAM{


template<class TV> class WATER_EXAMPLE;

template<class TV>
class WATER_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename ADVECTION_POLICY<GRID<TV> >::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef ARRAY<bool,FACE_INDEX<TV::m> > T_FACE_ARRAYS_BOOL;

protected:
    int current_frame;
    T time;
    int output_number;

    WATER_EXAMPLE<TV>& example;
    KINEMATIC_EVOLUTION<TV> kinematic_evolution;
public:
    THREAD_QUEUE* thread_queue;

    WATER_DRIVER(WATER_EXAMPLE<TV>& example);
    virtual ~WATER_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();
    void Advance_To_Target_Time(const T target_time);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Run(RANGE<TV_INT>& domain,const T dt,const T time);

//#####################################################################
};
}
#endif
