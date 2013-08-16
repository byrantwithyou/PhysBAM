//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_EVOLUTION
//#####################################################################
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Vectors/VECTOR.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_SYSTEM.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BACKWARD_EULER_EVOLUTION<TV>::
BACKWARD_EULER_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities_input)
    :SOLIDS_EVOLUTION<TV>(solids_parameters_input,solid_body_collection_input,example_forces_and_velocities_input),newtons_method(*new NEWTONS_METHOD<T>),
    minimization_system(*new BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>(solid_body_collection)),
    minimization_objective(*new BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>(solid_body_collection,minimization_system)),
    dv(static_cast<GENERALIZED_VELOCITY<TV>&>(*minimization_objective.v1.Clone_Default()))
{
    newtons_method.max_iterations=100000;
    newtons_method.max_krylov_iterations=2000;
    newtons_method.krylov_tolerance=1;
    newtons_method.fail_on_krylov_not_converged=false;
    newtons_method.tolerance=1e-5;
    newtons_method.angle_tolerance=1e-2;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BACKWARD_EULER_EVOLUTION<TV>::
~BACKWARD_EULER_EVOLUTION()
{
    delete &dv;
}
//#####################################################################
// Function Advance_One_Time_Step_Position
//#####################################################################
template<class TV> void BACKWARD_EULER_EVOLUTION<TV>::
Advance_One_Time_Step_Position(const T dt,const T time, const bool solids)
{
}
//#####################################################################
// Function Advance_One_Time_Step_Velocity
//#####################################################################
template<class TV> void BACKWARD_EULER_EVOLUTION<TV>::
Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids)
{
    LOG::SCOPE scope("Advance_One_Time_Step_Position");
    solid_body_collection.Print_Energy(time,0);

    minimization_objective.dt=dt;
    minimization_objective.time=time;
    minimization_system.dt=dt;
    minimization_system.time=time;
    minimization_objective.Reset();
    dv.Resize(minimization_objective.v1);
    dv*=0;
    minimization_objective.Test(dv,minimization_system);

    bool converged=newtons_method.Newtons_Method(minimization_objective,minimization_system,dv);
    PHYSBAM_ASSERT(converged);
// TODO for rigid bodies    R.Normalize(), update angular momentum

    solid_body_collection.Print_Energy(time+dt,1);
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class TV> void BACKWARD_EULER_EVOLUTION<TV>::
Initialize_Rigid_Bodies(const T frame_rate, const bool restart)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    // initialize kinematic object positions and velocities
    if(!restart){
        kinematic_evolution.Get_Current_Kinematic_Keyframes(1/frame_rate,time);
        kinematic_evolution.Set_External_Positions(rigid_body_collection.rigid_body_particles.frame,time);
        kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particles.twist,time,time);
        rigid_body_collection.Update_Angular_Momentum();
        for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i)){rigid_body_collection.rigid_body_particles.frame(i).r.Normalize();}}
}
//#####################################################################
// Function Use_CFL
//#####################################################################
template<class TV> bool BACKWARD_EULER_EVOLUTION<TV>::
Use_CFL() const
{
    return false;
}
//#####################################################################
namespace PhysBAM{
template class BACKWARD_EULER_EVOLUTION<VECTOR<float,1> >;
template class BACKWARD_EULER_EVOLUTION<VECTOR<float,2> >;
template class BACKWARD_EULER_EVOLUTION<VECTOR<float,3> >;
template class BACKWARD_EULER_EVOLUTION<VECTOR<double,1> >;
template class BACKWARD_EULER_EVOLUTION<VECTOR<double,2> >;
template class BACKWARD_EULER_EVOLUTION<VECTOR<double,3> >;
}
