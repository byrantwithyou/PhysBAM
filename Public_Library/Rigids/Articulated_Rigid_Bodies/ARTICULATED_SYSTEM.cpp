//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/UNION_FIND.h>
#include <Core/Utilities/DEBUG_CAST.h>
#include <Core/Vectors/TWIST.h>
#include <Core/Vectors/VECTOR.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_SYSTEM.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_VECTOR.h>
#include <Rigids/Joints/JOINT.h>
#include <Rigids/Joints/JOINT_MESH.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <iomanip>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ARTICULATED_SYSTEM<TV>::
ARTICULATED_SYSTEM(ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body_input)
    :BASE(false,false),articulated_rigid_body(articulated_rigid_body_input),break_loops(false),internal_x(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ARTICULATED_SYSTEM<TV>::
~ARTICULATED_SYSTEM()
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Initialize()
{
    const JOINT_MESH<TV>& joint_mesh=articulated_rigid_body.joint_mesh;
    prismatic_projection.Resize(joint_mesh.Size());
    angular_projection.Resize(joint_mesh.Size());
    location.Resize(joint_mesh.Size());
    effective_mass.Resize(joint_mesh.Size());
    intermediate_twists.Resize(articulated_rigid_body.rigid_body_collection.rigid_body_particles.Size());
    keep_joint.Remove_All();
    for(JOINT_ID j(0);j<joint_mesh.Size();j++) if(joint_mesh.Is_Active(j)){
        const RIGID_BODY<TV>& parent_body=*articulated_rigid_body.Parent(j);
        const RIGID_BODY<TV>& child_body=*articulated_rigid_body.Child(j);
        location(j)=joint_mesh(j)->Location(parent_body,child_body);
        prismatic_projection(j)=joint_mesh(j)->Prismatic_Projection_Matrix(parent_body.Frame());
        angular_projection(j)=joint_mesh(j)->Angular_Projection_Matrix(parent_body.Frame());
        MATRIX<T,TWIST<TV>::dimension> M,N;
        parent_body.Effective_Inertia_Inverse(M,location(j));
        child_body.Effective_Inertia_Inverse(N,location(j));
        effective_mass(j)=M+N;
        effective_mass(j).In_Place_Cholesky_Inverse();}

    if(break_loops){
        keep_joint.Resize(joint_mesh.Size());
        UNION_FIND<int> u(intermediate_twists.m);
        for(JOINT_ID j(0);j<joint_mesh.Size();j++) if(joint_mesh.Is_Active(j)){
            int p=articulated_rigid_body.Parent_Id(j),c=articulated_rigid_body.Child_Id(j);
            if(u.Find(p)==u.Find(c)) continue;
            keep_joint(j)=true;
            u.Union(c,p);}}
}
//#####################################################################
// Function Gather
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Gather(const ARRAY_VIEW<const TWIST<TV>,JOINT_ID> x,ARRAY_VIEW<TWIST<TV> > y) const
{
    const JOINT_MESH<TV>& joint_mesh=articulated_rigid_body.joint_mesh;
    for(JOINT_ID j(0);j<joint_mesh.Size();j++) if(joint_mesh.Is_Active(j)){
        if(break_loops && !keep_joint(j)) continue;
        y(articulated_rigid_body.Parent_Id(j))+=articulated_rigid_body.Parent(j)->Gather(x(j),location(j));
        y(articulated_rigid_body.Child_Id(j))-=articulated_rigid_body.Child(j)->Gather(x(j),location(j));}
}
//#####################################################################
// Function Inverse_Mass
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Inverse_Mass(ARRAY_VIEW<TWIST<TV> > x) const
{
    for(int i=0;i<articulated_rigid_body.rigid_body_collection.rigid_body_particles.Size();i++)
        x(i)=articulated_rigid_body.rigid_body_collection.Rigid_Body(i).Inertia_Inverse_Times(x(i));
}
//#####################################################################
// Function Scatter
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Scatter(const ARRAY_VIEW<const TWIST<TV> > x,ARRAY_VIEW<TWIST<TV>,JOINT_ID> y) const
{
    const JOINT_MESH<TV>& joint_mesh=articulated_rigid_body.joint_mesh;
    for(JOINT_ID j(0);j<joint_mesh.Size();j++) if(joint_mesh.Is_Active(j)){
        if(break_loops && !keep_joint(j)) continue;
        y(j)+=articulated_rigid_body.Parent(j)->Scatter(x(articulated_rigid_body.Parent_Id(j)),location(j));
        y(j)-=articulated_rigid_body.Child(j)->Scatter(x(articulated_rigid_body.Child_Id(j)),location(j));}
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result,bool transpose) const
{
    ARRAY<TWIST<TV>,JOINT_ID>& twist=debug_cast<ARTICULATED_VECTOR<TV>&>(result).v;
    const ARRAY<TWIST<TV>,JOINT_ID>& wrench=debug_cast<const ARTICULATED_VECTOR<TV>&>(x).v;
    intermediate_twists.Fill(TWIST<TV>());
    twist.Fill(TWIST<TV>());

    Gather(wrench,intermediate_twists);
    Inverse_Mass(intermediate_twists);
    Scatter(intermediate_twists,twist);
    if(break_loops)
        for(JOINT_ID j(0);j<articulated_rigid_body.joint_mesh.Size();j++) if(articulated_rigid_body.joint_mesh.Is_Active(j))
            if(!keep_joint(j))
                twist(j)=wrench(j);

    Kinetic_Energy();
}
//#####################################################################
// Function Kinetic_Energy
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Kinetic_Energy() const
{
    if(!internal_x) return;
    intermediate_twists.Fill(TWIST<TV>());
    Gather(internal_x->v,intermediate_twists);
    Inverse_Mass(intermediate_twists);
    intermediate_twists+=articulated_rigid_body.rigid_body_collection.rigid_body_particles.twist;
    T ke=0;
    for(int i=0;i<articulated_rigid_body.rigid_body_collection.rigid_body_particles.Size();i++){
        RIGID_BODY<TV>& rb=articulated_rigid_body.rigid_body_collection.Rigid_Body(i);
        if(rb.Has_Infinite_Inertia()) continue;
        TWIST<TV> wrench=rb.Inertia_Times(intermediate_twists(i));
        ke+=(T).5*(TV::Dot_Product(wrench.linear,intermediate_twists(i).linear)+TV::SPIN::Dot_Product(wrench.angular,intermediate_twists(i).angular));}
    LOG::cout<<std::setprecision(20)<<"system ke "<<ke<<std::endl;
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    const JOINT_MESH<TV>& joint_mesh=articulated_rigid_body.joint_mesh;
    ARRAY<TWIST<TV>,JOINT_ID>& twist=debug_cast<ARTICULATED_VECTOR<TV>&>(x).v;
    for(JOINT_ID j(0);j<joint_mesh.Size();j++) if(joint_mesh.Is_Active(j)){
        if(break_loops && !keep_joint(j)) continue;
        twist(j).linear=prismatic_projection(j)*twist(j).linear;
        twist(j).angular=angular_projection(j)*twist(j).angular;}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double ARTICULATED_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const JOINT_MESH<TV>& joint_mesh=articulated_rigid_body.joint_mesh;
    const ARRAY<TWIST<TV>,JOINT_ID>& xx=debug_cast<const ARTICULATED_VECTOR<TV>&>(x).v;
    const ARRAY<TWIST<TV>,JOINT_ID>& yy=debug_cast<const ARTICULATED_VECTOR<TV>&>(y).v;

    // TODO: Better inner product.
    double r=0;
    for(JOINT_ID j(0);j<joint_mesh.Size();j++)
        if(joint_mesh.Is_Active(j))
            r+=xx(j).Get_Vector().Dot(yy(j).Get_Vector());
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR ARTICULATED_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    return (T)sqrt(Inner_Product(x,x));
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void ARTICULATED_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
}
namespace PhysBAM{
template class ARTICULATED_SYSTEM<VECTOR<float,2> >;
template class ARTICULATED_SYSTEM<VECTOR<float,3> >;
template class ARTICULATED_SYSTEM<VECTOR<double,2> >;
template class ARTICULATED_SYSTEM<VECTOR<double,3> >;
}
