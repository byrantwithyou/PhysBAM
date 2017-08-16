//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MAC_DRIVER__
#define __MPM_MAC_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/PHASE_ID.h>
#include <climits>
#include <functional>
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
    void Apply_Viscosity();
    void Move_Particles();
    T Compute_Dt() const;
    T Max_Particle_Speed() const;
    T Grid_V_Upper_Bound() const;
    void Update_Simulated_Particles();
    void Compute_Poisson_Matrix();
    void Move_Mass_Momentum_Inside(PHASE& ph) const;
    void Move_Mass_Momentum_Inside_Nearest(PHASE& ph) const;
    template<class T2> void Fix_Periodic(ARRAY<T2,TV_INT>& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic(ARRAY<T2,FACE_INDEX<TV::m> >& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic_Accum(ARRAY<T2,TV_INT>& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic_Accum(ARRAY<T2,FACE_INDEX<TV::m> >& u,int ghost=INT_MAX) const;
    T Face_Fraction(const FACE_INDEX<TV::m>& face_index,const ARRAY<T,TV_INT>& phi) const;
    void Face_Fraction(const FACE_INDEX<TV::m>& face_index,ARRAY<PAIR<PHASE_ID,T> >& face_fractions) const;
    PAIR<PHASE_ID,typename TV::SCALAR> Phase_Of_Cell(const TV_INT& cell_index) const;
    T Density(const FACE_INDEX<TV::m>& face_index) const;

    void Bump_Particles();
    TV Nearest_Point_On_Surface(const TV& p,const PHASE& ph,
        const ARRAY<TV,TV_INT>& gradient,
        const ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& Hessian) const;

    void Apply_BC();
    int Allocate_Projection_System_Variable();
    void Compute_Laplacian(int var);
    void Compute_Gradient(int nvar);
    void Reseeding();
    T Phase_And_Phi(const TV& X,PHASE_ID& phase) const;
    void Step(std::function<void()> func,const char* name,bool dump_substep=true,bool do_step=true);
    void Invalidate_Particle(int p);
    bool Neumann_Boundary_Condition(const FACE_INDEX<TV::m>& face,ARRAY<T,PHASE_ID>& bc) const;
//#####################################################################
};
}
#endif
