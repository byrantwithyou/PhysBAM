//#####################################################################
// Copyright 2008, Michael Lentine, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NORMAL_JOINT
//#####################################################################
#ifndef __NORMAL_JOINT__
#define __NORMAL_JOINT__

#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Rigids/Joints/JOINT.h>
namespace PhysBAM{

template<class TV>
class NORMAL_JOINT:public JOINT<TV>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef JOINT<TV> BASE;using BASE::J;using BASE::J_inverse;using BASE::F_pj;

public:
    using BASE::Constraint_Matrix_Helper;
    T x_min,x_max;

    NORMAL_JOINT()
        :x_min((T)0),x_max((T)0)
    {}

    virtual ~NORMAL_JOINT();

//#####################################################################
    void Constrain_Prismatically(TV& translation) const override; // This is in joint space
    void Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const override; // This is in world space.
    void Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix=0) const override;
    void Update_State_From_Joint_Frame(const bool enforce_constraints=true) override;
    TV Prismatic_Component_Translation() const override;
    void Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix=0) const override;
//#####################################################################
};

template<class T_input>
class NORMAL_JOINT<VECTOR<T_input,1> >:public JOINT<VECTOR<T_input,1> >
{
    typedef T_input T;
public:
    NORMAL_JOINT()
    {PHYSBAM_FATAL_ERROR();}

    virtual ~NORMAL_JOINT();
};
}
#endif
