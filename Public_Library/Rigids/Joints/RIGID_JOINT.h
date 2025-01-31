//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_JOINT
//#####################################################################
#ifndef __RIGID_JOINT__
#define __RIGID_JOINT__

#include <Rigids/Joints/JOINT.h>
namespace PhysBAM{

template<class TV>
class RIGID_JOINT:public JOINT<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef JOINT<TV> BASE;typedef typename TV::SPIN T_SPIN;
    using BASE::J;using BASE::J_inverse;

    bool prismatic_component;
    TV prismatic_translation;

    RIGID_JOINT()
        :prismatic_component(false)
    {}

    virtual ~RIGID_JOINT();

    void Set_Prismatic_Component_Translation(const TV& translation)
    {prismatic_translation=translation;prismatic_component=true;}

    VECTOR<bool,T_SPIN::m> Angular_Constraints() const override
    {VECTOR<bool,T_SPIN::m> constrain(VECTOR<bool,T_SPIN::m>::All_Ones_Vector()); return constrain;}

//#####################################################################
    bool Has_Angular_Constraint() const override;
    void Constrain_Prismatically(TV& translation) const override;
    void Constrain_Relative_Angular_Velocity(const FRAME<TV>& parent_frame,T_SPIN& relative_angular_velocity) const override;
    void Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const override;
    TV Prismatic_Component_Translation() const override;
    void Constrain_Angles(T_SPIN& angles) const override;
    void Update_State_From_Joint_Frame(const bool enforce_constraints=true) override;
    void Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix=0) const override;
//#####################################################################
};
}
#endif
