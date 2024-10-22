//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_SYSTEM_MPI
//#####################################################################
#ifndef __SOLID_SYSTEM_MPI__
#define __SOLID_SYSTEM_MPI__
#include <Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Class SOLID_SYSTEM_MPI
//#####################################################################
template<class TV> class MPI_SOLID_FLUID;
template<class TV> class NEWMARK_EVOLUTION;
template<class TV>
class SOLID_SYSTEM_MPI:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;

public:
    typedef GENERALIZED_VELOCITY<TV> VECTOR_T;

    BACKWARD_EULER_SYSTEM<TV>& solid_system;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> >& fluid_mass;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> >& rigid_body_fluid_mass;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> > modified_mass,one_over_modified_mass;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> > modified_world_space_rigid_mass,modified_world_space_rigid_mass_inverse;
    ARRAY<SYMMETRIC_MATRIX<T,TV::SPIN::m> > modified_world_space_rigid_inertia_tensor;
    ARRAY<SYMMETRIC_MATRIX<T,TV::SPIN::m> > modified_world_space_rigid_inertia_tensor_inverse;
    MPI_SOLID_FLUID<TV>* mpi_solid_fluid;
    ARRAY<ARRAY<int> >& coupled_deformable_particle_indices;
    mutable ARRAY<ARRAY<TV> > recv_fluid_V_boundary_arrays;
    mutable ARRAY<ARRAY<TWIST<TV> > > recv_fluid_rigid_V_boundary_arrays;

    NEWMARK_EVOLUTION<TV>& newmark_evolution;

    SOLID_SYSTEM_MPI(BACKWARD_EULER_SYSTEM<TV>& solid_system_input,ARRAY<DIAGONAL_MATRIX<T,TV::m> >& fluid_mass_input,ARRAY<DIAGONAL_MATRIX<T,TV::m> >& rigid_body_fluid_mass_input,
        ARRAY<SYMMETRIC_MATRIX<T,TV::SPIN::m> >& modified_world_space_rigid_inertia_tensor_input,MPI_SOLID_FLUID<TV>* mpi_solid_fluid_input,
        ARRAY<ARRAY<int> >& coupled_deformable_particle_indices_input,NEWMARK_EVOLUTION<TV>& newmark_evolution_input,const int rigid_V_size,bool precondition=true);

    ~SOLID_SYSTEM_MPI();

    const T Solid_Sign() const
    {return (T)-1;}

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const override
    {R=V;}

    void Project(KRYLOV_VECTOR_BASE<T>& V) const override
    {solid_system.Project(V);}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const override
    {return solid_system.Convergence_Norm(R);}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const override {}

//#####################################################################
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& BV) const override;
    void Send_Generalized_Velocity_To_Fluid(const GENERALIZED_VELOCITY<TV>& V) const;
    void Get_Generalized_Velocity_From_Fluid(GENERALIZED_VELOCITY<TV>& V) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F,bool transpose=false) const override;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V0,const KRYLOV_VECTOR_BASE<T>& V1) const override;
//#####################################################################
};
}
#endif
