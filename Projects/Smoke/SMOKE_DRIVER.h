//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_DRIVER__
#define __SMOKE_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{


template<class TV> class SMOKE_EXAMPLE;

template<class TV>
class SMOKE_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
 
protected:
    int current_frame;
    T time;
    int output_number;
    int ghost;

    SMOKE_EXAMPLE<TV>& example;
public:
    SMOKE_DRIVER(SMOKE_EXAMPLE<TV>& example);
    virtual ~SMOKE_DRIVER();

    // EAPIC
    void Update_Particle_Weights();
    void Particle_To_Grid();
    void Grid_To_Particle();
    void Move_Particles(const T dt);

    void Add_Buoyancy_Force(const T dt,const T time);
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
    void Print_Max_Divergence(const char* str);
    T Max_Particle_Speed() const;
    
    T Max_C() const;
    TV Total_Particle_Linear_Momentum() const;
    T Total_Particle_Kinetic_Energy() const;

    T Max_Grid_Speed() const;
    TV Total_Grid_Momentum() const;
    T Total_Grid_Kinetic_Energy() const;
    T Total_Grid_Mass() const;
    T Total_Particle_Mass() const;
    T Min_Face_Mass()const;
    TV Total_Grid_Linear_Momentum() const;
//#####################################################################
};
}
#endif
