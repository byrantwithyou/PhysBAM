//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MAC_DRIVER__
#define __MPM_MAC_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{

template<class TV> class MPM_MAC_EXAMPLE;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class LAPLACE_UNIFORM;

template<class TV>
class MPM_MAC_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
public:

    int current_frame;
    int output_number;

    MPM_MAC_EXAMPLE<TV>& example;
    
    MPM_MAC_DRIVER(MPM_MAC_EXAMPLE<TV>& example);
    virtual ~MPM_MAC_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Update_Particle_Weights();
    void Particle_To_Grid();
    void Grid_To_Particle();
    void Pressure_Projection();
    void Apply_Forces();
    T Compute_Dt() const;
    T Max_Particle_Speed() const;
    T Grid_V_Upper_Bound() const;
    void Update_Simulated_Particles();
    void Print_Grid_Stats(const char* str,T dt);
    void Print_Particle_Stats(const char* str,T dt);
    void Print_Energy_Stats(const char* str);
//#####################################################################
};
}
#endif