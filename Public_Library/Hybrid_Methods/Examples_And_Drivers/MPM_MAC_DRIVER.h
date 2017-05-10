//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MAC_DRIVER__
#define __MPM_MAC_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Hybrid_Methods/Examples_And_Drivers/PHASE_ID.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/PHASE_ID.h>
#include <climits>
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
    typedef typename MPM_MAC_EXAMPLE<TV>::PHASE PHASE;
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
    void Prepare_Scatter();
    void Particle_To_Grid();
    void Particle_To_Grid(PHASE& ph) const;
    void Grid_To_Particle();
    void Grid_To_Particle(const PHASE& ph);
    void Build_Level_Sets();
    void Build_Level_Sets(PHASE& ph);
    void Pressure_Projection();
    void Apply_Forces();
    T Compute_Dt() const;
    T Max_Particle_Speed() const;
    T Max_Particle_Speed(const PHASE& ph) const;
    T Grid_V_Upper_Bound() const;
    T Grid_V_Upper_Bound(const PHASE& ph) const;
    void Update_Simulated_Particles();
    void Print_Grid_Stats(const char* str,T dt);
    void Print_Particle_Stats(const char* str,T dt);
    void Print_Energy_Stats(const char* str);
    void Compute_Poisson_Matrix();
    void Move_Mass_Momentum_Inside(PHASE& ph) const;
    void Move_Mass_Momentum_Inside_Nearest(PHASE& ph) const;
    template<class T2> void Fix_Periodic(ARRAY<T2,TV_INT>& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic(ARRAY<T2,FACE_INDEX<TV::m> >& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic_Accum(ARRAY<T2,TV_INT>& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic_Accum(ARRAY<T2,FACE_INDEX<TV::m> >& u,int ghost=INT_MAX) const;
private:
    T Face_Fraction(const FACE_INDEX<TV::m>& face_index,const ARRAY<T,TV_INT>& phi) const;
    void Face_Fraction(const FACE_INDEX<TV::m>& face_index,ARRAY<PAIR<PHASE_ID,T> >& face_fractions) const;
    PAIR<PHASE_ID,typename TV::SCALAR> Phase_Of_Cell(const TV_INT& cell_index) const;
    T Density(const FACE_INDEX<TV::m>& face_index) const;

    void Bump_Particles();
    TV Nearest_Point_On_Surface(const TV& p,const PHASE& ph,
        const ARRAY<TV,TV_INT>& gradient,
        const ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& Hessian) const;

    void Apply_BC(ARRAY<bool,FACE_INDEX<TV::m> >& psi_N);
    int Allocate_Projection_System_Variable(ARRAY<int,TV_INT>& cell_index,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N);
    void Compute_Laplacian(const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const ARRAY<int,TV_INT>& cell_index,int var);
    void Compute_Gradient(const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const ARRAY<int,TV_INT>& cell_index,int nvar);
    void Shift_Particle_Position_Periodic(TV shift);
    void Dump_Grid_ShiftTest(const std::string& var_name,const ARRAY<T,FACE_INDEX<TV::m> >& arr);
//#####################################################################
};
}
#endif
