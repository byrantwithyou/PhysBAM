//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROTATION_SPRINGS
//#####################################################################
#include <Core/Vectors/Dot_Product.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Forces_And_Torques/ROTATION_SPRINGS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ROTATION_SPRINGS<TV>::
ROTATION_SPRINGS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,SEGMENT_MESH& mesh)
    :RIGIDS_FORCES<TV>(rigid_body_collection),mpi_solids(0),mesh(mesh)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ROTATION_SPRINGS<TV>::
~ROTATION_SPRINGS()
{}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void ROTATION_SPRINGS<TV>::
Update_Position_Based_State(const T time)
{
    if(mpi_solids) PHYSBAM_NOT_IMPLEMENTED("MPI support");

    states.Remove_All();
    for(int s=0;s<mesh.elements.m;s++){
        int parent,child;mesh.elements(s).Get(parent,child);
        const FRAME<TV> &q_world_parent=rigid_body_collection.rigid_body_particles.frame(parent),&q_world_child=rigid_body_collection.rigid_body_particles.frame(child); // q_a_b means frame from b to a
        const FRAME<TV> &q_joint_parent=object_to_joint_frames(s)[0],&q_joint_child=object_to_joint_frames(s)[1];
        QUATERNION<T> q1=(q_world_parent.r*q_joint_parent.Inverse().r).Quaternion(); // parent joint to world
        QUATERNION<T> q2=(q_world_child.r*q_joint_child.Inverse().r).Quaternion(); // child joint to world
        QUATERNION<T> q=q2.Conjugate()*q1; // parent joint to child joint (note conjugate instead of inverse, since we have raw quaternions now)
        T_SPIN v_limit=sin(clamp_max(angle_limits(s),(T)pi));
        if(!abs(q.v).All_Less(v_limit)){
            STATE state;state.edge=s;
            T_SPIN v_error=q.v-clamp(q.v,-v_limit,v_limit);
            T_SPIN stress=stiffness(s)*v_error;
            state.torque=-(T).5*((q1.s*q2.s-Dot_Product(q1.v,q2.v))*stress+Dot_Product(q2.v,stress)*q1.v+Dot_Product(q1.v,stress)*q2.v
                +T_SPIN::Cross_Product((q1.s*q2.v+q2.s*q1.v),stress));
            states.Append(state);}}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void ROTATION_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    if(mpi_solids) PHYSBAM_NOT_IMPLEMENTED("MPI support");

    for(int i=0;i<states.m;i++){const STATE& state=states(i);
        const VECTOR<int,2>& nodes=mesh.elements(state.edge);
        rigid_F(nodes[0]).angular+=state.torque;
        rigid_F(nodes[1]).angular-=state.torque;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void ROTATION_SPRINGS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    if(mpi_solids) PHYSBAM_NOT_IMPLEMENTED("MPI support");

    for(int i=0;i<states.m;i++){const STATE& state=states(i);
        const VECTOR<int,2>& nodes=mesh.elements(state.edge);
        T_SPIN torque=damping(state.edge)*(rigid_V(nodes[1]).angular-rigid_V(nodes[0]).angular);
        rigid_F(nodes[0]).angular+=torque;
        rigid_F(nodes[1]).angular-=torque;}
}
//#####################################################################
template class ROTATION_SPRINGS<VECTOR<float,3> >;
template class ROTATION_SPRINGS<VECTOR<double,3> >;
}
