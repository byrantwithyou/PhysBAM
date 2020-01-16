//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include "FLUID_BC.h"
#include "FLUID_BOUNDARY_VECTOR.h"
#include "FLUID_SOLVER.h"
#include "FLUID_STATE.h"
#include "PARTITIONED_DRIVER.h"
#include "SOLID_BC.h"
#include "SOLID_BOUNDARY_VECTOR.h"
#include "SOLID_FLUID_INTERFACE.h"
#include "SOLID_SOLVER.h"
#include "SOLID_STATE.h"
namespace PhysBAM
{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTITIONED_DRIVER<TV>::
PARTITIONED_DRIVER()
{
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARTITIONED_DRIVER<TV>::
~PARTITIONED_DRIVER()
{
    delete fluid_new_state;
    delete fluid_old_state;
    delete fluid_prev_state;
    delete fluid_bc;

    delete solid_new_state;
    delete solid_old_state;
    delete solid_prev_state;
    delete solid_bc;
}

//#####################################################################
// Function Run
//#####################################################################
template<class TV> void PARTITIONED_DRIVER<TV>::
Run()
{
    Initialize();
    Write(0);
    for(int frame=1;frame<=last_frame;frame++)
    {
        Simulate_Frame(frame);
        Write(frame);
    }
}

//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void PARTITIONED_DRIVER<TV>::
Initialize()
{
    fluid_solver->Initialize();
    fluid_new_state=fluid_solver->Make_State();
    fluid_old_state=fluid_solver->Make_State();
    fluid_prev_state=fluid_solver->Make_State();
    fluid_bc=fluid_solver->Make_BC();

    solid_solver->Initialize();
    solid_new_state=solid_solver->Make_State();
    solid_old_state=solid_solver->Make_State();
    solid_prev_state=solid_solver->Make_State();
    solid_bc=solid_solver->Make_BC();

    time=0;
}

//#####################################################################
// Function Simulate_Frame
//#####################################################################
template<class TV> void PARTITIONED_DRIVER<TV>::
Simulate_Frame(int frame)
{
    T target_time=frame*frame_dt;
    bool done=false;
    for(int substep=1;!done;substep++)
    {
        T dt=Compute_Dt(target_time,done);
        Simulate_Time_Step(dt);
        time+=dt;
    }
    time=target_time;
}

//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> auto PARTITIONED_DRIVER<TV>::
Compute_Dt(T target_time,bool& done) const -> T
{
    T solid_dt=solid_solver->Compute_Dt(time);
    T fluid_dt=fluid_solver->Compute_Dt(time);
    T dt=std::min(solid_dt,fluid_dt);
    if(max_dt && dt>max_dt) dt=max_dt;
    if(min_dt && dt<min_dt) dt=min_dt;
    if(fixed_dt) dt=fixed_dt;
    if(time+dt*(T)1.01>=target_time){dt=target_time-time;done=true;}
    else if(time+2*dt>=target_time) dt=min(dt,(T).51*(target_time-time));
    return dt;
}

//#####################################################################
// Function Simulate_Time_Step
//#####################################################################
template<class TV> void PARTITIONED_DRIVER<TV>::
Simulate_Time_Step(T dt)
{
    solid_solver->Save(solid_old_state);
    fluid_solver->Save(fluid_old_state);

    SOLID_BOUNDARY_VECTOR<TV>* solid_vector=solid_solver->Make_Boundary_Vector();
    FLUID_BOUNDARY_VECTOR<TV>* fluid_vector=fluid_solver->Make_Boundary_Vector();

    interface->Compute_Coupling_Weights(solid_solver,fluid_solver,time,dt);
    interface->Get_Boundary(solid_vector);
    interface->Get_Boundary(fluid_vector);

    fluid_solver->Predict_Time_Step(time,dt);

    
    for(int i=0;i<max_subiterations;i++)
    {
//        interface->Compute_BC(fluid_solver,solid_bc,time,dt);
        fluid_solver->Get_Force(fluid_vector);
        interface->Distribute_Force(solid_vector,fluid_vector);
        solid_solver->Restore(solid_old_state);
        solid_solver->Simulate_Time_Step(solid_vector,time,dt);

        interface->Compute_Coupling_Weights(solid_solver,fluid_solver,time,dt);
        interface->Get_Boundary(solid_vector);
        interface->Get_Boundary(fluid_vector);

        if(use_bpp) BPP_Projection(dt); // NOTE: this chanages solid velocity

//        interface->Compute_BC(solid_solver,fluid_bc,time,dt);
        solid_solver->Get_Velocity(solid_vector);
        interface->Interpolate_Velocity(fluid_vector,solid_vector);
        fluid_solver->Restore(fluid_old_state);
        fluid_solver->Simulate_Time_Step(fluid_vector,time,dt);

        // Must wait for two states to be available to compute convergence
        if(i>0 && Is_Subiteration_Converged())
            break;

        solid_solver->Save(solid_prev_state);
        fluid_solver->Save(fluid_prev_state);
    }
}

//#####################################################################
// Function BPP_Projection
//#####################################################################
template<class TV> void PARTITIONED_DRIVER<TV>::
BPP_Projection(T dt)
{
    ARRAY<FLUID_BOUNDARY_VECTOR<TV>*> C;
    fluid_solver->Get_Constraints(C);

    if(p0.m!=C.m) p0.Resize(C.m,init_all);

    SOLID_BOUNDARY_VECTOR<TV>* solid_velocity=solid_solver->Make_Boundary_Vector();
    interface->Get_Boundary(solid_velocity);
    solid_solver->Get_Velocity(solid_velocity);

    ARRAY<SOLID_BOUNDARY_VECTOR<TV>*> Cs(C.m);
    ARRAY<SOLID_BOUNDARY_VECTOR<TV>*> Mi_Cs(C.m);
    ARRAY<T> rhs(C.m);

    for(int i=0;i<C.m;i++)
    {
        Cs(i)=solid_solver->Make_Boundary_Vector();
        interface->Distribute_Force(Cs(i),C(i));

        Mi_Cs(i)=solid_solver->Make_Boundary_Vector();
        solid_solver->Mass_Inverse(Mi_Cs(i),Cs(i));

        rhs(i)=-solid_solver->Inner_Product(Cs(i),solid_velocity);
    }

    MATRIX_MXN<T> M(C.m);
    for(int i=0;i<C.m;i++)
        for(int j=0;j<=i;j++)
            M(i,j)=M(j,i)=solid_solver->Inner_Product(Cs(i),Mi_Cs(j));

    ARRAY<T> dp0=M.Cholesky_Solve(rhs);

    p0+=dp0/dt;

    for(int i=0;i<dp0.m;i++)
        solid_solver->Apply_Velocity_Change(dp0(i),Mi_Cs(i));

    C.Delete_Pointers_And_Clean_Memory();
    Cs.Delete_Pointers_And_Clean_Memory();
    Mi_Cs.Delete_Pointers_And_Clean_Memory();
    delete solid_velocity;
}

//#####################################################################
// Function Is_Subiteration_Converged
//#####################################################################
template<class TV> bool PARTITIONED_DRIVER<TV>::
Is_Subiteration_Converged() const
{
    if(utol && fluid_solver->Diff_u(fluid_prev_state)>=utol) return false;
    if(ptol && fluid_solver->Diff_p(fluid_prev_state)>=ptol) return false;
    if(xtol && solid_solver->Diff_x(solid_prev_state)>=xtol) return false;
    if(vtol && solid_solver->Diff_v(solid_prev_state)>=vtol) return false;
    return true;
}

//#####################################################################
// Function Write
//#####################################################################
template<class TV> void PARTITIONED_DRIVER<TV>::
Write(int frame) const
{
    fluid_solver->Write(frame);
    solid_solver->Write(frame);
}

//#####################################################################
// Function Read
//#####################################################################
template<class TV> void PARTITIONED_DRIVER<TV>::
Read(int frame)
{
    fluid_solver->Read(frame);
    solid_solver->Read(frame);
}

template class PARTITIONED_DRIVER<VECTOR<float,2> >;
template class PARTITIONED_DRIVER<VECTOR<float,3> >;
template class PARTITIONED_DRIVER<VECTOR<double,2> >;
template class PARTITIONED_DRIVER<VECTOR<double,3> >;
}
