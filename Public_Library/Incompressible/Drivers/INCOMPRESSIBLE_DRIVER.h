//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INCOMPRESSIBLE_DRIVER__
#define __INCOMPRESSIBLE_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_PDE/Advection/ADVECTION_UNIFORM_FORWARD.h>
#include <Rigids/Rigids_Evolution/KINEMATIC_EVOLUTION.h>
namespace PhysBAM{


template<class TV> class INCOMPRESSIBLE_EXAMPLE;

template<class TV>
class INCOMPRESSIBLE_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

protected:
    int current_frame;
    T time;
    int output_number;

    INCOMPRESSIBLE_EXAMPLE<TV>& example;
    KINEMATIC_EVOLUTION<TV> kinematic_evolution;
    ARRAY<T,FACE_INDEX<TV::m> > sum_jc;
public:

    INCOMPRESSIBLE_DRIVER(INCOMPRESSIBLE_EXAMPLE<TV>& example);
    virtual ~INCOMPRESSIBLE_DRIVER();

    void Scalar_Advance(const T dt,const T time);
    void Convect(const T dt,const T time);
    void Add_Forces(const T dt,const T time);
    void Project(const T dt,const T time);
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
