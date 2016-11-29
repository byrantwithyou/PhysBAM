//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Jonathan Su, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANGLE_JOINT
//#####################################################################
#ifndef __ANGLE_JOINT__
#define __ANGLE_JOINT__

#include <Rigids/Joints/JOINT.h>
namespace PhysBAM{

template<class TV>
class ANGLE_JOINT:public JOINT<TV>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
public:
    using JOINT<TV>::J;using JOINT<TV>::J_inverse;using JOINT<TV>::F_pj;using JOINT<TV>::Constraint_Matrix_Helper;
    bool constrain_angle;
    T angle_min,angle_max;

    ANGLE_JOINT()
        :constrain_angle(false)
    {}

    virtual ~ANGLE_JOINT();

    void Set_Angle_Constraints(const bool use_constraints=true,const T angle_min_input=-pi,const T angle_max_input=pi)
    {constrain_angle=use_constraints;if(constrain_angle){angle_min=angle_min_input;angle_max=angle_max_input;}}

    void Constrain_Relative_Angular_Velocity_Helper(const FRAME<VECTOR<T,3> >& parent_frame,VECTOR<T,3>& relative_angular_velocity) const
    {relative_angular_velocity.Project_On_Unit_Direction((parent_frame.r*F_pj().r).Rotated_X_Axis());}

    void Constrain_Relative_Angular_Velocity_Helper(const FRAME<VECTOR<T,2> >& parent_frame,VECTOR<T,1>& relative_angular_velocity) const
    {}

    void Constrain_Relative_Angular_Velocity(const FRAME<TV>& parent_frame,T_SPIN& relative_angular_velocity) const override
    {Constrain_Relative_Angular_Velocity_Helper(parent_frame,relative_angular_velocity);}

    VECTOR<bool,T_SPIN::m> Angular_Constraints() const override
    {VECTOR<bool,T_SPIN::m> constrain(VECTOR<bool,T_SPIN::m>::All_Ones_Vector());constrain.x=constrain_angle && angle_min>=angle_max;return constrain;}

//#####################################################################
    bool Has_Angular_Constraint() const override;
    void Constrain_Prismatically(TV& translation) const override;
    void Constrain_Angles(T_SPIN& angles) const override;
    void Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix=0) const override;
    void Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const override;
    void Update_State_From_Joint_Frame(const bool enforce_constraints=true) override;
//#####################################################################
};

template<class T_input>
class ANGLE_JOINT<VECTOR<T_input,1> >:public JOINT<VECTOR<T_input,1> >
{
    typedef T_input T;
public:
    ANGLE_JOINT()
    {PHYSBAM_FATAL_ERROR();}

    virtual ~ANGLE_JOINT();
};
}
#endif
