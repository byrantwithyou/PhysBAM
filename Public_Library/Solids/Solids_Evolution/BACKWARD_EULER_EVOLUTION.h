//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_EVOLUTION
//#####################################################################
#ifndef __BACKWARD_EULER_EVOLUTION__
#define __BACKWARD_EULER_EVOLUTION__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Matrices/MATRIX_POLICY.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
namespace PhysBAM{

template<class T> struct NEWTONS_METHOD;
template<class TV> class BACKWARD_EULER_MINIMIZATION_OBJECTIVE;
template<class TV> class BACKWARD_EULER_SYSTEM;
template<class TV> class BACKWARD_EULER_MINIMIZATION_SYSTEM;
template<class TV> class GENERALIZED_VELOCITY;

template<class TV>
class BACKWARD_EULER_EVOLUTION:public SOLIDS_EVOLUTION<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef SOLIDS_EVOLUTION<TV> BASE;
  public:
    using BASE::krylov_vectors;using BASE::world_space_rigid_mass;using BASE::world_space_rigid_mass_inverse;
    using BASE::time;using BASE::solid_body_collection;using BASE::solids_parameters;using BASE::kinematic_evolution;
    using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::example_forces_and_velocities;
    using BASE::Set_External_Positions;using BASE::Set_External_Velocities;

    NEWTONS_METHOD<T>& newtons_method;
    BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>& minimization_system;
    BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>& minimization_objective;
    GENERALIZED_VELOCITY<TV>& dv,&tmp0,&tmp1;
    T coefficient_of_friction;
    bool fail_on_newton_not_converged;
    bool test_diff;

    BACKWARD_EULER_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities_input);
    virtual ~BACKWARD_EULER_EVOLUTION();

//#####################################################################
    bool Use_CFL() const override;
    void Advance_One_Time_Step_Position(const T dt,const T time,const bool solids) override;
    void Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids) override;
    void Initialize_Rigid_Bodies(const T frame_rate, const bool restart) override;
//#####################################################################
};
}
#endif
