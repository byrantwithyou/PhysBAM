//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_SYSTEM_MPI
//#####################################################################
#include <Dynamics/Parallel_Computation/FLUID_SYSTEM_MPI.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#ifdef USE_MPI
#include <Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#endif
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Vectors/Dot_Product.h>
#include <Core/Vectors/TWIST.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_SYSTEM_MPI<TV>::
FLUID_SYSTEM_MPI(const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_deformable_array_input,const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_rigid_array_input,
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array_input,const ARRAY<INTERVAL<int> >& interior_regions_input,const T tolerance_ratio_input,MPI_SOLID_FLUID<TV>* mpi_solid_fluid_input,
    GENERALIZED_VELOCITY<TV>& temp_input,GENERALIZED_VELOCITY<TV>& solid_velocity_input,ARRAY<int>& coupled_deformable_particle_indices_input,bool precondition)
    :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(precondition,!precondition),
    J_deformable_array(J_deformable_array_input),J_rigid_array(J_rigid_array_input),A_array(A_array_input),interior_regions(interior_regions_input),mpi_solid_fluid(mpi_solid_fluid_input),
    tolerance_ratio(tolerance_ratio_input),temp(temp_input),solid_velocity(solid_velocity_input),coupled_deformable_particle_indices(coupled_deformable_particle_indices_input)
{
    temp_array.Resize(A_array.m);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& bx,KRYLOV_VECTOR_BASE<T>& bresult,bool transpose) const
{
    const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);KRYLOV_VECTOR_T& result=debug_cast<KRYLOV_VECTOR_T&>(bresult);
    // get x values from solid
    Get_Generalized_Velocity_From_Solid(solid_velocity); // MPI
    temp.V.Fill(TV());temp.rigid_V.Fill(TWIST<TV>());
    for(int i=0;i<A_array.m;i++){
        const SPARSE_MATRIX_FLAT_MXN<T>& A=A_array(i);
        const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable=J_deformable_array(i);
        const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid=J_rigid_array(i);
        result.v(i).Resize(A.n);A.Times(x.v(i),result.v(i));

        Add_J_Rigid_Transpose_Times_Velocity(J_rigid,solid_velocity,result.v(i).Array_View(interior_regions(i)));
        Add_J_Deformable_Transpose_Times_Velocity(J_deformable,solid_velocity,result.v(i).Array_View(interior_regions(i)));

        Add_J_Rigid_Times_Pressure(J_rigid,x.v(i).Array_View(interior_regions(i)),temp);
        Add_J_Deformable_Times_Pressure(J_deformable,x.v(i).Array_View(interior_regions(i)),temp);}
    Send_Generalized_Velocity_To_Solid(temp); // MPI
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& bx,KRYLOV_VECTOR_BASE<T>& bresult) const // solve MR=V
{
    const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);KRYLOV_VECTOR_T& result=debug_cast<KRYLOV_VECTOR_T&>(bresult);
    for(int i=0;i<A_array.m;i++){
        const SPARSE_MATRIX_FLAT_MXN<T>& A=A_array(i);
        VECTOR_T x_i,result_i;
        temp_array(i).Resize(interior_regions(i).Size());
        A.C->Solve_Forward_Substitution(x.v(i).Array_View(interior_regions(i)),temp_array(i),true); // diagonal should be treated as the identity
        A.C->Solve_Backward_Substitution(temp_array(i),result.v(i).Array_View(interior_regions(i)),false,true);} // diagonal is inverted to save on divides
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double FLUID_SYSTEM_MPI<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV0,const KRYLOV_VECTOR_BASE<T>& BV1) const
{
    const KRYLOV_VECTOR_T& V0=debug_cast<const KRYLOV_VECTOR_T&>(BV0),&V1=debug_cast<const KRYLOV_VECTOR_T&>(BV1);
    double product=0;
    for(int i=0;i<V0.v.m;i++)
        product+=Dot_Product_Double_Precision(V0.v(i).Array_View(interior_regions(i)),V1.v(i).Array_View(interior_regions(i)));
    return product;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_SYSTEM_MPI<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bx) const
{
    const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);
    T fluid_convergence_norm=(T)0;
    for(int i=0;i<A_array.m;i++)
        fluid_convergence_norm=max(fluid_convergence_norm,x.v(i).Array_View(interior_regions(i)).Maximum_Magnitude());
    return tolerance_ratio*fluid_convergence_norm;
}
//#####################################################################
// Function Add_J_Deformable_Transpose_Times_Velocity
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI<TV>::
Add_J_Deformable_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> pressure) const
{
    assert(pressure.m==J_deformable.n && J_deformable.m==V.V.indices.Size()*TV::m);
    // computes pressure+=J*V.V
    J_deformable.For_Each(
        [&](int i,int j,T a){
            int dynamic_particle_index=i/TV::m;
            int axis=i-dynamic_particle_index*TV::m;
            pressure(j)+=a*V.V(dynamic_particle_index)(axis);});
}
//#####################################################################
// Function Add_J_Rigid_Transpose_Times_Velocity
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI<TV>::
Add_J_Rigid_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> pressure) const
{
    assert(pressure.m==J_rigid.n && J_rigid.m==V.rigid_V.indices.Size()*rows_per_rigid_body);
    // computes pressure+=J*V.V
    J_rigid.For_Each(
        [&](int i,int j,T a){
            int rigid_particle_index=i/rows_per_rigid_body;
            int component=i-rigid_particle_index*rows_per_rigid_body;
            if(component<=TV::m)
                pressure(j)+=a*V.rigid_V(rigid_particle_index).linear(component);
            else pressure(j)+=a*V.rigid_V(rigid_particle_index).angular(component-TV::m);});
}
//#####################################################################
// Function Add_J_Deformable_Times_Pressure
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI<TV>::
Add_J_Deformable_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,ARRAY_VIEW<const T> pressure,GENERALIZED_VELOCITY<TV>& V) const
{
    assert(pressure.m==J_deformable.n && J_deformable.m==V.V.indices.Size()*TV::m);
    J_deformable.For_Each(
        [&](int i,int j,T a){
            int dynamic_particle_index=i/TV::m;
            int axis=i-dynamic_particle_index*TV::m;
            V.V(dynamic_particle_index)(axis)+=pressure(j)*a;});
}
//#####################################################################
// Function Add_J_Rigid_Times_Pressure
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI<TV>::
Add_J_Rigid_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,ARRAY_VIEW<const T> pressure,GENERALIZED_VELOCITY<TV>& V) const
{
    assert(pressure.m==J_rigid.n && J_rigid.m==V.rigid_V.indices.Size()*rows_per_rigid_body);
    // computes pressure+=J*V.V
    J_rigid.For_Each(
        [&](int i,int j,T a){
            int rigid_particle_index=i/rows_per_rigid_body;
            int component=i-rigid_particle_index*rows_per_rigid_body;
            if(component<=TV::m)
                V.rigid_V(rigid_particle_index).linear(component)+=a*pressure(j);
            else V.rigid_V(rigid_particle_index).angular(component-TV::m)+=a*pressure(j);});
}
#ifdef USE_MPI
//#####################################################################
// Function void Send_Generalized_Velocity_To_Solid
//#####################################################################
template<class TV> void FLUID_SYSTEM_MPI<TV>::
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
template<class TV> void FLUID_SYSTEM_MPI<TV>::
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
template<class TV> void FLUID_SYSTEM_MPI<TV>::Send_Generalized_Velocity_To_Solid(const INDIRECT_ARRAY<const ARRAY_VIEW<TV> > V_boundary,const INDIRECT_ARRAY<const ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void FLUID_SYSTEM_MPI<TV>::Get_Generalized_Velocity_From_Solid(INDIRECT_ARRAY<ARRAY_VIEW<TV> > V_boundary,INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> > > rigid_V_boundary) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
#endif
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,d) \
    template class FLUID_SYSTEM_MPI<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
}
