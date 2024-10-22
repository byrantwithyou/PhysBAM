//#####################################################################
// Copyright 2008-2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_SYSTEM_MPI_SLIP
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Vectors/Dot_Product.h>
#include <Core/Vectors/TWIST.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Dynamics/Parallel_Computation/FLUID_SYSTEM_MPI_SLIP.h>
#ifdef USE_MPI
#include <Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_SYSTEM_MPI_SLIP<TV>::
FLUID_SYSTEM_MPI_SLIP(const bool use_preconditioner_input,const SPARSE_MATRIX_FLAT_MXN<T>& W_input,const SPARSE_MATRIX_FLAT_MXN<T>& N_input,const SPARSE_MATRIX_FLAT_MXN<T>& div_input,
    const ARRAY<T>& M_inverse_input,const SPARSE_MATRIX_FLAT_MXN<T>& J_input,const SPARSE_MATRIX_FLAT_MXN<T>& P_input,const SPARSE_MATRIX_FLAT_MXN<T>& div_precondition,
    const bool leakproof_solve_input,const bool using_slip,const INTERVAL<int>& interior_regions_input,const INTERVAL<int>& divergence_indices_input,
    MPI_SOLID_FLUID_SLIP<TV>* mpi_solid_fluid_input,ARRAY<int>& coupled_deformable_particle_indices_input,GENERALIZED_VELOCITY<TV>& solid_velocity_input)
    :BASE(use_preconditioner_input,!use_preconditioner_input),mpi_solid_fluid(mpi_solid_fluid_input),coupled_deformable_particle_indices(coupled_deformable_particle_indices_input),
    J(J_input),N(N_input),div(div_input),PP(P_input),W(W_input),M_inverse(M_inverse_input),leakproof_solve(leakproof_solve_input),interior_regions(interior_regions_input),
    solid_velocity(solid_velocity_input),divergence_indices(divergence_indices_input)
{
    N.Transpose(N_transpose);
    div.Transpose(div_transpose);
    
    C_f=div+PP*N*W;
    C_f.Transpose(C_f_transpose);

    pressures_size_vector.Resize(div.m);
    if(!leakproof_solve){
        J.Transpose(J_transpose);
        C_s=PP*N*J;
        C_s*=-1;
        C_s.Transpose(C_s_transpose);
        // want to zero out columns outside the local box
        for(int i=0;i<C_s_transpose.A.m;i++)
            if(interior_regions.Lazy_Outside_Half_Open(C_s_transpose.A(i).j))
                C_s_transpose.A(i).a=0;
        solid_velocities_size_vector.Resize(J.n);}
    if(use_preconditioner){
        SPARSE_MATRIX_FLAT_MXN<T> div_precondition_transpose;
        div_precondition.Transpose(div_precondition_transpose);
        div_M_inverse_div_transpose_precondition=(div_precondition*div_precondition_transpose.Scale_Rows(M_inverse));
        delete div_M_inverse_div_transpose_precondition.C;
        div_M_inverse_div_transpose_precondition.C=div_M_inverse_div_transpose_precondition.Create_Submatrix(divergence_indices);
        const bool modified_incomplete_cholesky=true;
        const T modified_incomplete_cholesky_coefficient=(T).97;
        const T preconditioner_zero_tolerance=(T)1e-8;
        const T preconditioner_zero_replacement=(T)1e-8;
        div_M_inverse_div_transpose_precondition.C->In_Place_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
            preconditioner_zero_tolerance,preconditioner_zero_replacement);
        preconditioned_pressures_size_vector.Resize(div_precondition.m);
    }
    fluid_velocities_size_vector.Resize(M_inverse.Size());
    lagrange_multipliers_size_vector.Resize(C_f.m);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI_SLIP<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF,bool transpose) const
{
    const KRYLOV_VECTOR_T& V=debug_cast<const KRYLOV_VECTOR_T&>(BV);KRYLOV_VECTOR_T& F=debug_cast<KRYLOV_VECTOR_T&>(BF);
    // get x values from solid
    solid_velocity.V.Fill(TV());solid_velocity.rigid_V.Fill(TWIST<TV>());
    C_f_transpose.Times(V.v,fluid_velocities_size_vector);
    fluid_velocities_size_vector*=M_inverse;
    C_f.Times(fluid_velocities_size_vector,F.v);

    if(!leakproof_solve){
        // Get the parts of solid V we need
        Get_Generalized_Velocity_From_Solid(solid_velocity); // MPI
        solid_velocity.Pack(solid_velocities_size_vector);
        C_s.Times(solid_velocities_size_vector,lagrange_multipliers_size_vector);

        F.v+=lagrange_multipliers_size_vector;

        C_s_transpose.Times(V.v,solid_velocities_size_vector);
        solid_velocity.Unpack(solid_velocities_size_vector);

        // Send back our contribution to F
        Send_Generalized_Velocity_To_Solid(solid_velocity);} // MPI
}
//#####################################################################
// Function Apply
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI_SLIP<TV>::
Apply(const VECTOR_T& V,ARRAY<T>& result_dual_cells_size_vector) const
{
    result_dual_cells_size_vector.Resize(div_transpose.m);
    C_f_transpose.Times(V,fluid_velocities_size_vector);
    fluid_velocities_size_vector *= M_inverse;
    result_dual_cells_size_vector=fluid_velocities_size_vector;
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI_SLIP<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BR) const // solve MR=V
{
    const KRYLOV_VECTOR_T& V=debug_cast<const KRYLOV_VECTOR_T&>(BV);KRYLOV_VECTOR_T& R=debug_cast<KRYLOV_VECTOR_T&>(BR);
    R.Copy((T)1,V);
    pressures_size_vector.Resize(divergence_indices.Size());
    div_M_inverse_div_transpose_precondition.C->Solve_Forward_Substitution(V.v.Array_View(divergence_indices),pressures_size_vector,true);
    div_M_inverse_div_transpose_precondition.C->Solve_Backward_Substitution(pressures_size_vector,R.v.Array_View(divergence_indices),false,true);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double FLUID_SYSTEM_MPI_SLIP<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV0,const KRYLOV_VECTOR_BASE<T>& BV1) const
{
    const KRYLOV_VECTOR_T& V0=debug_cast<const KRYLOV_VECTOR_T&>(BV0);const KRYLOV_VECTOR_T& V1=debug_cast<const KRYLOV_VECTOR_T&>(BV1);
    return Dot_Product_Double_Precision(V0.v.Array_View(interior_regions),V1.v.Array_View(interior_regions));
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_SYSTEM_MPI_SLIP<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BV) const
{
    const KRYLOV_VECTOR_T& V=debug_cast<const KRYLOV_VECTOR_T&>(BV);
    return V.v.Array_View(interior_regions).Maximum_Magnitude();
}
#ifdef USE_MPI
//#####################################################################
// Function void Send_Generalized_Velocity_To_Solid
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI_SLIP<TV>::
Send_Generalized_Velocity_To_Solid(const INDIRECT_ARRAY<const ARRAY_VIEW<TV> > V_boundary,const INDIRECT_ARRAY<const ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const
{
    // TODO: optimize this by only communicating the solid boundary velocities
    int tag_fluid_to_solid=mpi_solid_fluid->Get_Unique_Tag();
    int buffer_size=MPI_UTILITIES::Pack_Size(V_boundary,rigid_V_boundary,*mpi_solid_fluid->comm)+1;
    ARRAY<char> buffer_send(buffer_size);int position=0;
    MPI_UTILITIES::Pack(V_boundary,rigid_V_boundary,buffer_send,position,*mpi_solid_fluid->comm);
    mpi_solid_fluid->comm->Send(&buffer_send(0),position,MPI::PACKED,mpi_solid_fluid->solid_node,tag_fluid_to_solid);
}
//#####################################################################
// Function void Send_Generalized_Velocity_To_Solid
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI_SLIP<TV>::
Get_Generalized_Velocity_From_Solid(INDIRECT_ARRAY<ARRAY_VIEW<TV> > V_boundary,INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const
{
    // TODO: optimize this by only communicating the solid boundary velocities
    int tag_solid_to_fluid=mpi_solid_fluid->Get_Unique_Tag();
    int buffer_size=MPI_UTILITIES::Pack_Size(V_boundary,rigid_V_boundary,*mpi_solid_fluid->comm)+1;
    ARRAY<char> buffer(buffer_size);int position=0;
    mpi_solid_fluid->comm->Recv(&buffer(0),buffer_size,MPI::PACKED,mpi_solid_fluid->solid_node,tag_solid_to_fluid);
    MPI_UTILITIES::Unpack(V_boundary,rigid_V_boundary,buffer,position,*mpi_solid_fluid->comm);
}
#else
template<class TV> void FLUID_SYSTEM_MPI_SLIP<TV>::Send_Generalized_Velocity_To_Solid(const INDIRECT_ARRAY<const ARRAY_VIEW<TV> > V_boundary,const INDIRECT_ARRAY<const ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void FLUID_SYSTEM_MPI_SLIP<TV>::Get_Generalized_Velocity_From_Solid(INDIRECT_ARRAY<ARRAY_VIEW<TV> > V_boundary,INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
#endif
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,d) \
    template class FLUID_SYSTEM_MPI_SLIP<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
}
