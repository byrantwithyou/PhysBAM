//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_MXN
//#####################################################################
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GENERALIZED_VELOCITY<TV>::
GENERALIZED_VELOCITY(const SOLID_BODY_COLLECTION<TV>& solid_body_collection)
    :GENERALIZED_VELOCITY(solid_body_collection.deformable_body_collection.particles.V,
        solid_body_collection.rigid_body_collection.rigid_body_particles.twist,
        solid_body_collection)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GENERALIZED_VELOCITY<TV>::
GENERALIZED_VELOCITY(ARRAY_VIEW<TV> V_full,ARRAY_VIEW<TWIST<TV> > rigid_V_full,const SOLID_BODY_COLLECTION<TV>& solid_body_collection)
    :V(V_full,solid_body_collection.deformable_body_collection.dynamic_particles),rigid_V(rigid_V_full,solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles),
    kinematic_and_static_rigid_V(rigid_V_full,solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies),deep_copy(false)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GENERALIZED_VELOCITY<TV>::
GENERALIZED_VELOCITY(ARRAY_VIEW<TV> V_full,const ARRAY<int>& dynamic_particles,ARRAY_VIEW<TWIST<TV> > rigid_V_full,
    const ARRAY<int>& dynamic_rigid_body_particles,const ARRAY<int>& static_and_kinematic_rigid_bodies)
    :V(V_full,dynamic_particles),rigid_V(rigid_V_full,dynamic_rigid_body_particles),
    kinematic_and_static_rigid_V(rigid_V_full,static_and_kinematic_rigid_bodies),deep_copy(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GENERALIZED_VELOCITY<TV>::
~GENERALIZED_VELOCITY()
{
    if(deep_copy){
        delete [] V.array.Get_Array_Pointer();
        delete [] rigid_V.array.Get_Array_Pointer();}
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> GENERALIZED_VELOCITY<TV>& GENERALIZED_VELOCITY<TV>::
operator=(const GENERALIZED_VELOCITY& gv)
{
    V=gv.V;
    rigid_V=gv.rigid_V;
    kinematic_and_static_rigid_V=gv.kinematic_and_static_rigid_V;
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& GENERALIZED_VELOCITY<TV>::
operator+=(const BASE& bv)
{
    const GENERALIZED_VELOCITY& v=debug_cast<const GENERALIZED_VELOCITY&>(bv);
    V+=v.V;
    rigid_V+=v.rigid_V;
    kinematic_and_static_rigid_V+=v.kinematic_and_static_rigid_V;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& GENERALIZED_VELOCITY<TV>::
operator-=(const BASE& bv)
{
    const GENERALIZED_VELOCITY& v=debug_cast<const GENERALIZED_VELOCITY&>(bv);
    V-=v.V;
    rigid_V-=v.rigid_V;
    kinematic_and_static_rigid_V-=v.kinematic_and_static_rigid_V;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& GENERALIZED_VELOCITY<TV>::
operator*=(const T a)
{
    V*=a;
    rigid_V*=a;
    kinematic_and_static_rigid_V*=a;
    return *this;
}
//#####################################################################
// Function Set_Zero
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Set_Zero()
{
    V.Fill(TV());
    rigid_V.Fill(TWIST<TV>());
    kinematic_and_static_rigid_V.Fill(TWIST<TV>());
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Copy(const T c,const BASE& bv)
{
    const GENERALIZED_VELOCITY& v=debug_cast<const GENERALIZED_VELOCITY&>(bv);
    V=c*v.V;
    rigid_V=c*v.rigid_V;
    kinematic_and_static_rigid_V=c*v.kinematic_and_static_rigid_V;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const GENERALIZED_VELOCITY& v1=debug_cast<const GENERALIZED_VELOCITY&>(bv1),&v2=debug_cast<const GENERALIZED_VELOCITY&>(bv2);
    V=c1*v1.V+v2.V;
    rigid_V=c1*v1.rigid_V+v2.rigid_V;
    kinematic_and_static_rigid_V=c1*v1.kinematic_and_static_rigid_V+v2.kinematic_and_static_rigid_V;
}
//#####################################################################
// Function Pack
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Pack(ARRAY<T> &velocities) const
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=0;i<V.Size();i++)
        for(int j=0;j<TV::m;j++)
            velocities(index++)=V.array(i)(j);
    for(int i=0;i<rigid_V.Size();i++){
        for(int j=0;j<TV::m;j++)
            velocities(index++)=rigid_V(i).linear(j);
        for(int j=0;j<T_SPIN::m;j++)
            velocities(index++)=rigid_V(i).angular(j);}
}
//#####################################################################
// Function Unpack
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Unpack(ARRAY<T> &velocities)
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=0;i<V.Size();i++)
        for(int j=0;j<TV::m;j++)
            V.array(i)(j)=velocities(index++);
    for(int i=0;i<rigid_V.Size();i++){
        for(int j=0;j<TV::m;j++)
            rigid_V(i).linear(j)=velocities(index++);
        for(int j=0;j<T_SPIN::m;j++)
            rigid_V(i).angular(j)=velocities(index++);}
}
//#####################################################################
// Function Unpack_And_Add
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Unpack_And_Add(ARRAY<T> &velocities)
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=0;i<V.Size();i++)
        for(int j=0;j<TV::m;j++)
            V.array(i)(j)+=velocities(index++);
    for(int i=0;i<rigid_V.Size();i++){
        for(int j=0;j<TV::m;j++)
            rigid_V(i).linear(j)+=velocities(index++);
        for(int j=0;j<T_SPIN::m;j++)
            rigid_V(i).angular(j)+=velocities(index++);}
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int GENERALIZED_VELOCITY<TV>::
Raw_Size() const
{
    return V.Size()*TV::m+rigid_V.Size()*TWIST<TV>::dimension;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& GENERALIZED_VELOCITY<TV>::
Raw_Get(int i)
{
    if(i<V.Size()*TV::m) return V(i/TV::m)(i%TV::m);
    i-=V.Size()*TV::m;
    int o=i%TWIST<TV>::dimension,n=i/TWIST<TV>::dimension;
    if(o<TV::m) return rigid_V(n).linear(o);
    return rigid_V(n).angular(o-TV::m);
}
//#####################################################################
// Function Exchange
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Exchange(GENERALIZED_VELOCITY<TV>& gv)
{
    V.array.Exchange(gv.V.array);
    rigid_V.array.Exchange(gv.rigid_V.array);
    kinematic_and_static_rigid_V.array.Exchange(gv.kinematic_and_static_rigid_V.array);
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* GENERALIZED_VELOCITY<TV>::
Clone_Default() const
{
    GENERALIZED_VELOCITY<TV>* gv=new GENERALIZED_VELOCITY<TV>(ARRAY_VIEW<TV>(new TV[V.array.m],V.array.m),V.indices,
        ARRAY_VIEW<TWIST<TV> >(new TWIST<TV>[rigid_V.array.m],rigid_V.array.m),rigid_V.indices,kinematic_and_static_rigid_V.indices);
    gv->deep_copy=true;
    return gv;
}
//#####################################################################
// Function Resize_Helper
//#####################################################################
template<class AV>
void Resize_Helper(AV& a,const AV& b)
{
    if(a.Size()==b.Size()) return;
    if(a.Size()>b.Size()){a.m=b.m;return;}
    AV t(new typename AV::ELEMENT[b.m],b.m);
    t.Fill(typename AV::ELEMENT());
    a.Exchange(t);
    delete [] t.Get_Array_Pointer();
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& v)
{
    if(!deep_copy) return;
    const GENERALIZED_VELOCITY<TV>& gv=debug_cast<const GENERALIZED_VELOCITY<TV>&>(v);
    Resize_Helper(V.array,gv.V.array);
    Resize_Helper(rigid_V.array,gv.rigid_V.array);
    kinematic_and_static_rigid_V.array.Set(rigid_V.array);
}
//#####################################################################
// Function Get
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Get(ARRAY_VIEW<T> a) const
{
    int k=0;
    for(int i=0;i<V.Size();i++)
        for(const auto& j:V(i))
            a(k++)=j;
    for(int i=0;i<rigid_V.Size();i++)
    {
        for(const auto& j:rigid_V(i).linear)
            a(k++)=j;
        for(const auto& j:rigid_V(i).angular)
            a(k++)=j;
    }
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Set(ARRAY_VIEW<const T> a)
{
    int k=0;
    for(int i=0;i<V.Size();i++)
        for(auto& j:V(i))
            j=a(k++);
    for(int i=0;i<rigid_V.Size();i++)
    {
        for(auto& j:rigid_V(i).linear)
            j=a(k++);
        for(auto& j:rigid_V(i).angular)
            j=a(k++);
    }
}
//#####################################################################
namespace PhysBAM{
template class GENERALIZED_VELOCITY<VECTOR<float,1> >;
template class GENERALIZED_VELOCITY<VECTOR<float,2> >;
template class GENERALIZED_VELOCITY<VECTOR<float,3> >;
template class GENERALIZED_VELOCITY<VECTOR<double,1> >;
template class GENERALIZED_VELOCITY<VECTOR<double,2> >;
template class GENERALIZED_VELOCITY<VECTOR<double,3> >;
}
