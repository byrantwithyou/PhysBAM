//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_DRIVER__
#define __MPM_DRIVER__
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class MPM_EXAMPLE;
template<class TV> class MPM_OBJECTIVE;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class MPM_KRYLOV_VECTOR;

template<class TV>
class MPM_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:

    int current_frame;
    int output_number;

    MPM_EXAMPLE<TV>& example;
    MPM_OBJECTIVE<TV>& objective;
    MPM_KRYLOV_VECTOR<TV>& dv,&rhs;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    MPM_DRIVER(MPM_EXAMPLE<TV>& example);
    virtual ~MPM_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Update_Particle_Weights();
    void Particle_To_Grid();
    void Grid_To_Particle();
    void Apply_Forces();
    void Perform_Particle_Collision(int p);
    void Apply_Friction();
//#####################################################################
};
}
#endif
