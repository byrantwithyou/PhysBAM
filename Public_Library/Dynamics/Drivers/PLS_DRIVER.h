//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_DRIVER__
#define __PLS_DRIVER__
#include <Tools/Grids_Uniform_Advection/ADVECTION_POLICY_UNIFORM.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{


template<class TV> class PLS_EXAMPLE;

template<class TV>
class PLS_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename ADVECTION_POLICY<TV>::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;

protected:
    int current_frame;
    T time;
    int output_number;

    PLS_EXAMPLE<TV>& example;
public:

    PLS_DRIVER(PLS_EXAMPLE<TV>& example);
    virtual ~PLS_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();
    void Advance_To_Target_Time(const T target_time);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);

//#####################################################################
};
}
#endif
