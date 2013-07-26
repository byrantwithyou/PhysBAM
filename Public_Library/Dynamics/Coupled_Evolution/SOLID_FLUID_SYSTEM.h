//#####################################################################
// Copyright 2007-2008, Avi (Snarky) Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_SYSTEM
//#####################################################################
#ifndef __SOLID_FLUID_SYSTEM__
#define __SOLID_FLUID_SYSTEM__
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Dynamics/Coupled_Evolution/PRESSURE_VELOCITY_VECTOR.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class TV> class BACKWARD_EULER_SYSTEM;
template<class TV> class GRID;

//#####################################################################
// Class SOLID_FLUID_SYSTEM
//#####################################################################
template<class TV,class T_MATRIX>
class SOLID_FLUID_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<PAIR<int,T> > FACE_WEIGHT_ELEMENTS;typedef ARRAY<FACE_WEIGHT_ELEMENTS*,FACE_INDEX<TV::m> > T_FACE_ARRAYS_FACE_WEIGHT_ELEMENTS;
    typedef typename TV::SPIN T_SPIN;
    typedef KRYLOV_SYSTEM_BASE<typename TV::SCALAR> BASE;
public:
    typedef PRESSURE_VELOCITY_VECTOR<TV> VECTOR_T;
    static const int rows_per_rigid_body=TV::dimension+T_SPIN::dimension;
    using BASE::use_preconditioner;

    BACKWARD_EULER_SYSTEM<TV>& solid_system;
    const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_deformable_array;
    const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_rigid_array;
    const ARRAY<DIAGONAL_MATRIX<T,TV::m> >& fluid_mass;
    const ARRAY<DIAGONAL_MATRIX<T,TV::m> >& rigid_body_fluid_mass;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> > modified_mass,one_over_modified_mass;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> > modified_world_space_rigid_mass,modified_world_space_rigid_mass_inverse;
    const ARRAY<SYMMETRIC_MATRIX<T,TV::SPIN::m> >& modified_world_space_rigid_inertia_tensor;
    ARRAY<SYMMETRIC_MATRIX<T,TV::SPIN::m> > modified_world_space_rigid_inertia_tensor_inverse;
    const T fluid_tolerance,solid_tolerance;
    const ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array;
    mutable ARRAY<ARRAY<T> > temp_array;

    SOLID_FLUID_SYSTEM(BACKWARD_EULER_SYSTEM<TV>& solid_system_input,const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_deformable_array_input,
        const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_rigid_array_input,const ARRAY<DIAGONAL_MATRIX<T,TV::m> >& fluid_mass_input,
        const ARRAY<DIAGONAL_MATRIX<T,TV::m> >& rigid_body_fluid_mass_input,const ARRAY<SYMMETRIC_MATRIX<T,TV::SPIN::m> >& modified_world_space_rigid_inertia_tensor_input,
        const T fluid_tolerance_input,const T solid_tolerance_input,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array_input);

    virtual ~SOLID_FLUID_SYSTEM();

    double Check_Positive_Definite(VECTOR_T& V,VECTOR_T& F,VECTOR_T& R) const
    {
        Multiply(F,R);return Inner_Product(V,R);
    }

    const T Solid_Sign() const
    {return (T)-1;}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE {} // TODO

//#####################################################################
    // void Print_Matrix(VECTOR_T& V,VECTOR_T& F);
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;  // solve MR=V
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V1,const KRYLOV_VECTOR_BASE<T>& V2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    static void Add_J_Deformable_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,const GENERALIZED_VELOCITY<TV>& V,ARRAY<T>& pressure);
    static void Add_J_Rigid_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> pressure);
    static void Add_J_Deformable_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,const ARRAY<T>& pressure,GENERALIZED_VELOCITY<TV>& V);
    static void Add_J_Rigid_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,const ARRAY<T>& pressure,GENERALIZED_VELOCITY<TV>& V);
//#####################################################################
};
}
#endif
