//#####################################################################
// Copyright 2008-2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_SYSTEM_MPI_SLIP
//#####################################################################
#ifndef __FLUID_SYSTEM_MPI_SLIP__
#define __FLUID_SYSTEM_MPI_SLIP__
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <Tools/Krylov_Solvers/PCG_SPARSE_SYSTEM.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID_SLIP.h>
namespace PhysBAM{
//#####################################################################
// Class FLUID_SYSTEM_MPI_SLIP
//#####################################################################
template<class TV>
class FLUID_SYSTEM_MPI_SLIP:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef ARRAY<T> VECTOR_T;
    typedef KRYLOV_VECTOR_WRAPPER<T,VECTOR_T&> KRYLOV_VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
public:
    using BASE::use_preconditioner;

    MPI_SOLID_FLUID_SLIP<TV>* mpi_solid_fluid;
    ARRAY<int>& coupled_deformable_particle_indices;

    const SPARSE_MATRIX_FLAT_MXN<T>& J;
    const SPARSE_MATRIX_FLAT_MXN<T>& N;
    const SPARSE_MATRIX_FLAT_MXN<T>& div;
    const SPARSE_MATRIX_FLAT_MXN<T>& PP;
    const SPARSE_MATRIX_FLAT_MXN<T>& W;
    const ARRAY<T>& M_inverse;

    SPARSE_MATRIX_FLAT_MXN<T> J_transpose;
    SPARSE_MATRIX_FLAT_MXN<T> N_transpose;
    SPARSE_MATRIX_FLAT_MXN<T> div_transpose;

    SPARSE_MATRIX_FLAT_MXN<T> C_s,C_f;
    SPARSE_MATRIX_FLAT_MXN<T> C_s_transpose,C_f_transpose;
    SPARSE_MATRIX_FLAT_MXN<T> C_f_M_inverse_C_f_transpose;
    SPARSE_MATRIX_FLAT_MXN<T> div_M_inverse_div_transpose_precondition;

    mutable ARRAY<T> fluid_velocities_size_vector;
    mutable ARRAY<T> solid_velocities_size_vector;
    mutable ARRAY<T> pressures_size_vector;
    mutable ARRAY<T> preconditioned_pressures_size_vector;
    mutable ARRAY<T> lagrange_multipliers_size_vector;

    const bool leakproof_solve;
    INTERVAL<int> interior_regions;
    GENERALIZED_VELOCITY<TV>& solid_velocity;
    INTERVAL<int> divergence_indices;

    FLUID_SYSTEM_MPI_SLIP(const bool use_preconditioner_input,
        const SPARSE_MATRIX_FLAT_MXN<T>& W_input,
        const SPARSE_MATRIX_FLAT_MXN<T>& N_input,
        const SPARSE_MATRIX_FLAT_MXN<T>& div_input,
        const ARRAY<T>& M_inverse_input,
        const SPARSE_MATRIX_FLAT_MXN<T>& J_input,
        const SPARSE_MATRIX_FLAT_MXN<T>& P_input,
        const SPARSE_MATRIX_FLAT_MXN<T>& div_precondition,
        const bool leakproof_solve_input,
        const bool using_slip,
        const INTERVAL<int>& interior_regions_input,
        const INTERVAL<int>& divergence_indices_input,
        MPI_SOLID_FLUID_SLIP<TV>* mpi_solid_fluid_input,
        ARRAY<int>& coupled_deformable_particle_indices_input,
        GENERALIZED_VELOCITY<TV>& solid_velocity_input);

    void Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF,bool transpose=false) const override;
    void Apply(const VECTOR_T& V,ARRAY<T>& result_dual_cells_size_vector) const;

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const override // only nullspace stuff for fluids - leave out for now
    {}

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BR) const override; // solve MR=V
    
    void Project(KRYLOV_VECTOR_BASE<T>& x) const override
    {}

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV0,const KRYLOV_VECTOR_BASE<T>& BV1) const override;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BV) const override;

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const override {}

    void Send_Generalized_Velocity_To_Solid(const INDIRECT_ARRAY<const ARRAY_VIEW<TV> > V_boundary,const INDIRECT_ARRAY<const ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const;
    void Send_Generalized_Velocity_To_Solid(const GENERALIZED_VELOCITY<TV>& V) const
    {Send_Generalized_Velocity_To_Solid(V.V.array.Subset(coupled_deformable_particle_indices),V.rigid_V);}

    void Get_Generalized_Velocity_From_Solid(INDIRECT_ARRAY<ARRAY_VIEW<TV> > V_boundary,INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const;
    void Get_Generalized_Velocity_From_Solid(GENERALIZED_VELOCITY<TV>& V) const
    {Get_Generalized_Velocity_From_Solid(V.V.array.Subset(coupled_deformable_particle_indices),V.rigid_V);}
//#####################################################################
};
}
#endif
