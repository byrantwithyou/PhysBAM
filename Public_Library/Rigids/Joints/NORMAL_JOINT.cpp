//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/MATRIX_1X1.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Rigids/Joints/NORMAL_JOINT.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> NORMAL_JOINT<TV>::
~NORMAL_JOINT()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_input> NORMAL_JOINT<VECTOR<T_input,1> >::
~NORMAL_JOINT()
{}
//#####################################################################
// Function Constrain_Prismatically
//#####################################################################
template<class TV> void NORMAL_JOINT<TV>::
Constrain_Prismatically(TV& translation) const
{
    translation(0)=clamp(translation(0),x_min,x_max);
}
//#####################################################################
// Function Constrain_Relative_Linear_Velocity
//#####################################################################
template<class TV> void NORMAL_JOINT<TV>::
Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const
{
    ROTATION<TV> joint_orientation=parent_frame.r*F_pj().r;
    TV u=joint_orientation.Rotated_Axis(0);
    relative_linear_velocity-=TV::Dot_Product(relative_linear_velocity,u)*u;
}
//#####################################################################
// Function Angular_Constraint_Matrix
//#####################################################################
template<class TV> void NORMAL_JOINT<TV>::
Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix) const
{
    angular_constraint_matrix.Resize(T_SPIN::dimension,0);
    if(angular_unconstrained_matrix) (*angular_unconstrained_matrix)=MATRIX<T,T_SPIN::dimension>::Identity_Matrix();
}
//#####################################################################
// Function Update_State_From_Joint_Frame
//#####################################################################
template<class TV> void NORMAL_JOINT<TV>::
Update_State_From_Joint_Frame(const bool enforce_constraints)
{
    if(enforce_constraints){
        Constrain_Prismatically(J.t);
        J_inverse=J.Inverse();}
}
//#####################################################################
// Function Prismatic_Component_Translation
//#####################################################################
template<class TV> TV NORMAL_JOINT<TV>::
Prismatic_Component_Translation() const
{
    TV current_translation=J.t;
    Constrain_Prismatically(current_translation);
    return current_translation;
}
//#####################################################################
// Function Prismatic_Constraint_Matrix
//#####################################################################
template<class TV> void NORMAL_JOINT<TV>::
Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix) const
{
    TV_BOOL constrain;
    constrain(0)=x_min>=x_max;
    Constraint_Matrix_Helper(parent_frame.r*F_pj().r,constrained_matrix,unconstrained_matrix,constrain);
}
//#####################################################################
namespace PhysBAM{
template class NORMAL_JOINT<VECTOR<float,1> >;
template class NORMAL_JOINT<VECTOR<float,2> >;
template class NORMAL_JOINT<VECTOR<float,3> >;
template class NORMAL_JOINT<VECTOR<double,1> >;
template class NORMAL_JOINT<VECTOR<double,2> >;
template class NORMAL_JOINT<VECTOR<double,3> >;
}
