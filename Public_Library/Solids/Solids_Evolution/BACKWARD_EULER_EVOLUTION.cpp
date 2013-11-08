//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_EVOLUTION
//#####################################################################
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Vectors/VECTOR.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/LAGGED_FORCE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
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
    minimization_system(*new BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>(solid_body_collection,&example_forces_and_velocities_input)),
    minimization_objective(*new BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>(solid_body_collection,minimization_system)),
    dv(static_cast<GENERALIZED_VELOCITY<TV>&>(*minimization_objective.v1.Clone_Default())),
    tmp0(static_cast<GENERALIZED_VELOCITY<TV>&>(*minimization_objective.v1.Clone_Default())),
    tmp1(static_cast<GENERALIZED_VELOCITY<TV>&>(*minimization_objective.v1.Clone_Default()))
{
    newtons_method.max_iterations=100000;
    newtons_method.max_krylov_iterations=2000;
    newtons_method.krylov_tolerance=1;
    newtons_method.fail_on_krylov_not_converged=false;
    newtons_method.tolerance=1e-5;
    newtons_method.angle_tolerance=1e-2;
    minimization_system.tmp=&tmp1;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BACKWARD_EULER_EVOLUTION<TV>::
~BACKWARD_EULER_EVOLUTION()
{
    delete &dv;
    delete &tmp0;
    delete &tmp1;
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
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    solid_body_collection.Print_Energy(time,0);

    Set_External_Positions(particles.X,time);
    Set_External_Velocities(particles.V,time,time);

    minimization_objective.dt=dt;
    minimization_objective.time=time;
    minimization_system.dt=dt;
    minimization_system.time=time;
    minimization_objective.Reset();
    dv.Resize(minimization_objective.v1);
    tmp0.Resize(minimization_objective.v1);
    tmp1.Resize(minimization_objective.v1);

    for(int i=0;i<solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
        if(LAGGED_FORCE<TV>* lf=dynamic_cast<LAGGED_FORCE<TV>*>(solid_body_collection.deformable_body_collection.deformables_forces(i)))
            lf->Lagged_Update_Position_Based_State(time);

    minimization_objective.Initial_Guess(dv);
    minimization_objective.Test(dv,minimization_system);

    bool converged=newtons_method.Newtons_Method(minimization_objective,minimization_system,dv);
    PHYSBAM_ASSERT(converged);
// TODO for rigid bodies    R.Normalize(), update angular momentum

    minimization_objective.Adjust_For_Collision(dv);
    minimization_objective.Compute_Unconstrained(dv,0,&tmp0,0);
    solid_body_collection.Print_Energy(time+dt,1);
    tmp1=tmp0;
    minimization_objective.Project_Gradient_And_Prune_Constraints(tmp1,true);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before friction",1,1);
    for(int i=0;i<minimization_system.collisions.m;i++){
        const typename BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::COLLISION& c=minimization_system.collisions(i);
        T normal_force=c.n.Dot(tmp0.V.array(c.p)-tmp1.V.array(c.p));
        TV& v=minimization_objective.v1.V.array(c.p);
        TV t=v.Projected_Orthogonal_To_Unit_Direction(c.n);
        T t_mag=t.Normalize();
        if(t_mag<=minimization_objective.coefficient_of_friction(c.object)*normal_force/particles.mass(c.p))
            v.Project_On_Unit_Direction(c.n);
        else v-=coefficient_of_friction/particles.mass(c.p)*normal_force*t;}
    solid_body_collection.Print_Energy(time+dt,2);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after friction",1,1);

    solid_body_collection.Print_Energy(time+dt,3);
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
