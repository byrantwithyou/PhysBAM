//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_SYSTEM_MPI
//#####################################################################
#ifndef __FLUID_SYSTEM_MPI__
#define __FLUID_SYSTEM_MPI__
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
namespace PhysBAM{
//#####################################################################
// Class FLUID_SYSTEM_MPI
//#####################################################################
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class TV> class MPI_SOLID_FLUID;
template<class TV> class GENERALIZED_VELOCITY;
template<class TV>
class FLUID_SYSTEM_MPI:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;

public:
    typedef ARRAY<T> VECTOR_T;
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_T> > KRYLOV_VECTOR_T;
    static const int rows_per_rigid_body=TV::m+T_SPIN::m;

    const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_deformable_array;
    const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_rigid_array;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array;
    ARRAY<INTERVAL<int> > interior_regions;
    MPI_SOLID_FLUID<TV>* mpi_solid_fluid;
    T tolerance_ratio;
    GENERALIZED_VELOCITY<TV>& temp;
    GENERALIZED_VELOCITY<TV>& solid_velocity;
    mutable ARRAY<ARRAY<T> > temp_array;
    ARRAY<int>& coupled_deformable_particle_indices;

    FLUID_SYSTEM_MPI(const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_deformable_array_input,const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_rigid_array_input,
        ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array_input,const ARRAY<INTERVAL<int> >& interior_regions_input,const T tolerance_ratio_input,MPI_SOLID_FLUID<TV>* mpi_solid_fluid_input,
        GENERALIZED_VELOCITY<TV>& temp_input,GENERALIZED_VELOCITY<TV>& solid_velocity_input,ARRAY<int>& coupled_deformable_particle_indices_input,bool precondition);

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const override // only nullspace stuff for fluids - leave out for now
    {}

    void Project(KRYLOV_VECTOR_BASE<T>& x) const override
    {}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const override {}

    void Send_Generalized_Velocity_To_Solid(const GENERALIZED_VELOCITY<TV>& V) const
    {Send_Generalized_Velocity_To_Solid(V.V.array.Subset(coupled_deformable_particle_indices),V.rigid_V);}

    void Get_Generalized_Velocity_From_Solid(GENERALIZED_VELOCITY<TV>& V) const
    {Get_Generalized_Velocity_From_Solid(V.V.array.Subset(coupled_deformable_particle_indices),V.rigid_V);}

//#####################################################################
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result,bool transpose=false) const override;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const override; // solve MR=V
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V0,const KRYLOV_VECTOR_BASE<T>& V1) const override;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const override;
    void Send_Generalized_Velocity_To_Solid(const INDIRECT_ARRAY<const ARRAY_VIEW<TV> > V_boundary,const INDIRECT_ARRAY<const ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const;
    void Get_Generalized_Velocity_From_Solid(INDIRECT_ARRAY<ARRAY_VIEW<TV> > V_boundary,INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const;
private:
    void Add_J_Deformable_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> pressure) const;
    void Add_J_Rigid_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> pressure) const;
    void Add_J_Deformable_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,ARRAY_VIEW<const T> pressure,GENERALIZED_VELOCITY<TV>& V) const;
    void Add_J_Rigid_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,ARRAY_VIEW<const T> pressure,GENERALIZED_VELOCITY<TV>& V) const;
//#####################################################################
};
}
#endif
