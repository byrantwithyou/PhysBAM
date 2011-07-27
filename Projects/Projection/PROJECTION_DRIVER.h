//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PROJECTION_DRIVER__
#define __PROJECTION_DRIVER__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{


template<class TV> class PROJECTION_EXAMPLE;

template<class TV>
class PROJECTION_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

protected:
    int current_frame;
    T time;
    int output_number;

    PROJECTION_EXAMPLE<TV>& example;
public:
    PROJECTION_DRIVER(PROJECTION_EXAMPLE<TV>& example);
    virtual ~PROJECTION_DRIVER();

    void Scalar_Advance(const T dt,const T time);
    void Project(const T dt,const T time);
    void Execute_Main_Program();
    void Initialize();
    void Advance_To_Target_Time(const T target_time);
    void Simulate_To_Frame(const int frame_input);
    void Advance_One_Time_Step_Implicit_Part(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T dt,const T time);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);

//#####################################################################
};
}
#endif
