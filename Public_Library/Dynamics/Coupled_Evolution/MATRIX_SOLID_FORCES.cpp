//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_FORCES
//##################################################################### 
#include <Core/Arrays/CONSTANT_ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Dynamics/Coupled_Evolution/MATRIX_SOLID_FORCES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_SOLID_FORCES<TV>::
MATRIX_SOLID_FORCES(const SOLID_BODY_COLLECTION<TV>& collection)
    :solid_body_collection(collection),dt(0),current_velocity_time(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_SOLID_FORCES<TV>::
~MATRIX_SOLID_FORCES()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Compute(const T dt_input,const T current_velocity_time_input)
{
    dt=dt_input;
    current_velocity_time=current_velocity_time_input;

    force_dof_counts.Resize(solid_body_collection.deformable_body_collection.deformables_forces.m);force_dof_counts.Fill(0);
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++){
        PHYSBAM_ASSERT(solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces);
        force_dof_counts(k)=solid_body_collection.deformable_body_collection.deformables_forces(k)->Velocity_Dependent_Forces_Size();}

    total_force_dof=FORCE_AGGREGATE_ID(force_dof_counts.Sum());
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Times_Add(const GENERALIZED_VELOCITY<TV>& solid_velocities,ARRAY<T,FORCE_AGGREGATE_ID>& force_coefficients) const
{
    //if(mpi_solids){mpi_solids->Exchange_Binding_Boundary_Data(V.V.array);/*TODO: mpi_solids->Exchange_Boundary_Data(V.rigid_V.array);*/}
    //F.V.array.Subset(dynamic_particles).Fill(TV());F.rigid_V.array.Subset(dynamic_rigid_body_particles).Fill(TWIST<TV>());
    //solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities(V.rigid_V.array);
    //solid_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(V.V.array,V.rigid_V.array);
    //if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(V.V.array);
    // solid_body_collection.Implicit_Velocity_Independent_Forces(V.V.array,F.V.array,dt,current_velocity_time+dt);
    // we don't inherit velocity dependent forces to the drifted particles (contrary to the explicit velocity independent ones) because that would compromise symmetry
    Add_Velocity_Dependent_Forces_First_Half(solid_velocities,force_coefficients,current_velocity_time+dt);
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Times(const GENERALIZED_VELOCITY<TV>& solid_velocities,ARRAY<T,FORCE_AGGREGATE_ID>& force_coefficients) const
{
    force_coefficients.Fill(T());
    Times_Add(solid_velocities,force_coefficients);
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Transpose_Times_Add(const ARRAY<T,FORCE_AGGREGATE_ID>& force_coefficients,GENERALIZED_VELOCITY<TV>& solid_velocities) const
{
    Add_Velocity_Dependent_Forces_Second_Half(force_coefficients,solid_velocities,current_velocity_time+dt);
}
//#####################################################################
// Function Transpose_Times
//#####################################################################
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Transpose_Times(const ARRAY<T,FORCE_AGGREGATE_ID>& force_coefficients,GENERALIZED_VELOCITY<TV>& solid_velocities) const
{
    // TODO: Careful to zero out enough of the solids state.
    solid_velocities.V.Fill(TV());
    solid_velocities.rigid_V.Fill(TWIST<TV>());
    solid_velocities.kinematic_and_static_rigid_V.Fill(TWIST<TV>());
    Transpose_Times_Add(force_coefficients,solid_velocities);
    //if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data(F.V.array);
    //solid_body_collection.binding_list.Distribute_Force_To_Parents(F.V.array,F.rigid_V.array);
    //solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Distribute_Force_To_Parents(F.rigid_V.array);
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
// can depend on position too
template<class TV> FORCE_AGGREGATE_ID MATRIX_SOLID_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    return total_force_dof;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
// can depend on position too
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T,FORCE_AGGREGATE_ID> aggregate,const T time) const
{
    FORCE_AGGREGATE_ID aggregate_size(0);
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++){
        PHYSBAM_ASSERT(solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces);
        ARRAY_VIEW<T> single_force_aggregate=aggregate.Array_View(aggregate_size,force_dof_counts(k));
        aggregate_size+=force_dof_counts(k);
        solid_body_collection.deformable_body_collection.deformables_forces(k)->Add_Velocity_Dependent_Forces_First_Half(V.V.array,single_force_aggregate,time);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
// can depend on position too
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T,FORCE_AGGREGATE_ID> aggregate,GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    FORCE_AGGREGATE_ID aggregate_size(0);
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++){
        PHYSBAM_ASSERT(solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces);
        //ARRAY_VIEW<const T> single_force_aggregate=aggregate.Array_View(aggregate_size,single_aggregate_size);
        ARRAY_VIEW<const T> single_force_aggregate(aggregate.Get_Array_Pointer()+Value(aggregate_size),force_dof_counts(k));
        aggregate_size+=force_dof_counts(k);
        solid_body_collection.deformable_body_collection.deformables_forces(k)->Add_Velocity_Dependent_Forces_Second_Half(single_force_aggregate,F.V.array,time);}
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Test_Matrix() const
{
    RANDOM_NUMBERS<T> random;
    FORCE_AGGREGATE_ID size=Velocity_Dependent_Forces_Size();
    ARRAY<T,FORCE_AGGREGATE_ID> aggregate(size),aggregate2(size);
    random.Fill_Uniform(aggregate,-1,1);

    GENERALIZED_VELOCITY<TV> pattern(solid_body_collection);
    GENERALIZED_VELOCITY<TV>& GV0=(GENERALIZED_VELOCITY<TV>&)*pattern.Clone_Default();
    GENERALIZED_VELOCITY<TV>& GV1=(GENERALIZED_VELOCITY<TV>&)*pattern.Clone_Default();
    GENERALIZED_VELOCITY<TV>& GV2=(GENERALIZED_VELOCITY<TV>&)*pattern.Clone_Default();
    GENERALIZED_VELOCITY<TV>& GV3=(GENERALIZED_VELOCITY<TV>&)*pattern.Clone_Default();
    GENERALIZED_VELOCITY<TV>& GV4=(GENERALIZED_VELOCITY<TV>&)*pattern.Clone_Default();

    int n=GV0.Raw_Size();
    ARRAY<T> r0(n),r1(n),r2(n),r3(n),r4(n);
    random.Fill_Uniform(r0,-1,1);
    GV0.Set(r0);

    Add_Velocity_Dependent_Forces_Second_Half(aggregate,GV1,0);
    Add_Velocity_Dependent_Forces_First_Half(GV0,aggregate2,0);
    Add_Velocity_Dependent_Forces_Second_Half(aggregate2,GV2,0);
    solid_body_collection.Add_Velocity_Dependent_Forces(GV0,GV3,0);
    solid_body_collection.Add_Velocity_Dependent_Forces(GV1,GV4,0);
    GV1.Get(r1);
    GV2.Get(r2);
    GV3.Get(r3);
    GV4.Get(r4);

    T inner_solids=r0.Dot(r1);
    T inner_aggregate=aggregate.Dot(aggregate2);
    LOG::cout<<"MATRIX_SOLID_FORCES Symmetry Test: "<<inner_solids<<"  vs  "<<inner_aggregate<<"  relative  "<<abs(inner_solids-inner_aggregate)/maxabs((T)1e-30,inner_solids,inner_aggregate)<<std::endl;

    T inner1=r0.Dot(r0);
    T inner3=r2.Dot(r2);
    T inner4=r3.Dot(r3);
    T inner_def=r0.Dot(r3);
    r2+=r3;
    T innerD=r2.Dot(r2);
    LOG::cout<<"MATRIX_SOLID_FORCES C vs D test: "<<inner3<<"  vs  "<<inner4<<"   diff "<<innerD<<"   orig "<<inner1<<"  relative  "<<abs(innerD)/maxabs((T)1e-30,inner3,inner4)<<std::endl;
    LOG::cout<<"MATRIX_SOLID_FORCES D definiteness (should be negative): "<<inner_def<<std::endl;

    T solids_inner1=r0.Dot(r4);
    T solids_inner2=r1.Dot(r3);
    LOG::cout<<"MATRIX_SOLID_FORCES D Symmetry Test: "<<solids_inner1<<"  vs  "<<solids_inner2<<"  relative  "<<abs(solids_inner1-solids_inner2)/maxabs((T)1e-30,solids_inner1,solids_inner2)<<std::endl;

    delete & GV0;
    delete & GV1;
    delete & GV2;
    delete & GV3;
    delete & GV4;
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Print_Each_Matrix(int n,const GENERALIZED_VELOCITY<TV>& V,T sqrt_dt) const
{
    OCTAVE_OUTPUT<T> oo(LOG::sprintf("C-%i.txt",n).c_str());

    GENERALIZED_VELOCITY<TV> G(V);
    int c=G.Raw_Size();
    oo.Begin_Sparse_Matrix("C",Value(total_force_dof),c);
    G*=0;
    ARRAY<T,FORCE_AGGREGATE_ID> force_coefficients(total_force_dof);
    for(int i=0;i<c;i++){
        G.Raw_Get(i)=1;
        Times(G,force_coefficients);
        G.Raw_Get(i)=0;
        force_coefficients*=sqrt_dt;
        oo.Append_Sparse_Column(ARRAY_VIEW<T,FORCE_AGGREGATE_ID>(force_coefficients));}

    oo.End_Sparse_Matrix();
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_FORCES<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    int offset=0;
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++){
        PHYSBAM_ASSERT(solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces);
        int start=data.m;
        solid_body_collection.deformable_body_collection.deformables_forces(k)->Add_Raw_Velocity_Dependent_Forces_First_Half(data);
        for(int i=start;i<data.m;i++) data(i).x+=offset;
        offset+=Value(force_dof_counts(k));}

    // TODO: Correct indices for indirect mapping.
    PHYSBAM_ASSERT(solid_body_collection.deformable_body_collection.dynamic_particles.m==solid_body_collection.deformable_body_collection.particles.V.m);
}
//#####################################################################
namespace PhysBAM{
template class MATRIX_SOLID_FORCES<VECTOR<float,1> >;
template class MATRIX_SOLID_FORCES<VECTOR<float,2> >;
template class MATRIX_SOLID_FORCES<VECTOR<float,3> >;
template class MATRIX_SOLID_FORCES<VECTOR<double,1> >;
template class MATRIX_SOLID_FORCES<VECTOR<double,2> >;
template class MATRIX_SOLID_FORCES<VECTOR<double,3> >;
}
