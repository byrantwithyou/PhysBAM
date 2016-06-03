//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar, Jonathan Su
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRISMATIC_TWIST_JOINT
//#####################################################################
#ifndef __PRISMATIC_TWIST_JOINT__
#define __PRISMATIC_TWIST_JOINT__

#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Rigids/Joints/ANGLE_JOINT.h>
namespace PhysBAM{

template<class TV>
class PRISMATIC_TWIST_JOINT:public ANGLE_JOINT<TV>
{
    typedef typename TV::SCALAR T;typedef ANGLE_JOINT<TV> BASE;typedef typename TV::SPIN T_SPIN;
    using BASE::J;using BASE::F_pj;

public:
    using BASE::Constraint_Matrix_Helper;
    VECTOR<bool,TV::m> constrain;
    TV prismatic_min,prismatic_max;

    PRISMATIC_TWIST_JOINT()
    {constrain.x=true;}

    virtual ~PRISMATIC_TWIST_JOINT();

    void Set_Prismatic_Constraints(const VECTOR<bool,TV::m>& constrain_input,const TV& min=TV(),const TV& max=TV())
    {constrain=constrain_input;prismatic_min=min;prismatic_max=max;}

private:
    VECTOR<bool,TV::m> Equality_Constraint() const
    {return constrain.Componentwise_And(prismatic_min.Componentwise_Greater_Equal(prismatic_max));}
public:
//#####################################################################
    bool Has_Prismatic_Constraint() const override;
    void Constrain_Prismatically(TV& translation) const override; // This is in joint space
    void Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const override; // This is in world space.
    TV Prismatic_Component_Translation() const override;
    void Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix=0) const override;
//#####################################################################
};

template<class T_input>
class PRISMATIC_TWIST_JOINT<VECTOR<T_input,1> >:public JOINT<VECTOR<T_input,1> >
{
    typedef T_input T;typedef JOINT<VECTOR<T,1> > BASE;
public:
    PRISMATIC_TWIST_JOINT()
    {PHYSBAM_FATAL_ERROR();}

    virtual ~PRISMATIC_TWIST_JOINT();
};
}
#endif
