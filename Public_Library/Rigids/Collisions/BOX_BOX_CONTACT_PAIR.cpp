//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONTACT_PAIRS
//##################################################################### 
#include <Core/Math_Tools/RANGE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
#include <Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <Rigids/Collisions/SOLVE_CONTACT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGIDS_COLLISION_CALLBACKS;

namespace CONTACT_PAIRS
{
template<class TV>
bool Update_Box_Box_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    const int id_1,const int id_2,IMPLICIT_OBJECT<TV>* object1,IMPLICIT_OBJECT<TV>* object2,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;

    FRAME<TV> transform1,transform2;
    RIGID_BODY<TV>& body0=rigid_body_collection.Rigid_Body(id_1);
    RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id_2);
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object1)){
        transform1=*object_transformed->transform;object1=object_transformed->object_space_implicit_object;}
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object2)){
        transform2=*object_transformed->transform;object2=object_transformed->object_space_implicit_object;}
    RANGE<TV>& box1=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >&>(*object1).analytic;
    RANGE<TV>& box2=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >&>(*object2).analytic;
    ORIENTED_BOX<TV> box1_transformed(box1,body0.Frame());
    ORIENTED_BOX<TV> box2_transformed(box2,body1.Frame());

    TV collision_normal=body0.Frame().t-body1.Frame().t;collision_normal.Normalize();
    if(!box1_transformed.Intersection(box2_transformed)){rigid_body_collisions.skip_collision_check.Set_Last_Checked(id_1,id_2);return false;}
    if(TV::Dot_Product(collision_normal,body0.Twist().linear-body1.Twist().linear)>=0) return false;
    TV collision_location=(body0.Frame().t-body1.Frame().t)*.5+body0.Frame().t;
    TV collision_relative_velocity=body0.Pointwise_Object_Velocity(collision_location)-body1.Pointwise_Object_Velocity(collision_location);
    collision_callbacks.Swap_States(id_1,id_2);
    SOLVE_CONTACT::Update_Contact_Pair_Helper<TV>(rigid_body_collisions,collision_callbacks,id_1,id_2,dt,time,epsilon_scale,collision_location,collision_normal,collision_relative_velocity,
        correct_contact_energy,rigid_body_collisions.rolling_friction,mpi_one_ghost);
    return true;
}

#define INSTANTIATION_HELPER(T,d) \
    template bool Update_Box_Box_Contact_Pair(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        const int id_1,const int id_2,IMPLICIT_OBJECT<VECTOR<T,d> >* object1,IMPLICIT_OBJECT<VECTOR<T,d> >* object2,const bool correct_contact_energy,const int max_iterations, \
        const T epsilon_scale,const T dt,const T time,const bool mpi_one_ghost);

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
}
}
