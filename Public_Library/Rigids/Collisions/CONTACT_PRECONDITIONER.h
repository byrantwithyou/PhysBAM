//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONTACT_PRECONDITIONER
//##################################################################### 
#ifndef __CONTACT_PRECONDITIONER__
#define __CONTACT_PRECONDITIONER__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Core/Vectors/TWIST.h>
#include <Rigids/Collisions/PROJECTED_GAUSS_SEIDEL.h>
#include <Rigids/Collisions/SOLVE_CONTACT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>

namespace PhysBAM{

namespace CONTACT_PRECONDITIONER
{

template<class TV>
class PRECONDITIONER_RIGID_BODY
{
public:
    typedef typename TV::SCALAR T;
    FRAME<TV> frame;
    //TWIST<TV> twist;
    T mass;
    SYMMETRIC_MATRIX<T,TV::SPIN::m> inertia_tensor;
    
    //bool has_infinite_inertia;

    PRECONDITIONER_RIGID_BODY() {}
    PRECONDITIONER_RIGID_BODY(RIGID_BODY<TV>& body)
        :frame(body.Frame()),/*twist(body.Twist()),*/mass(body.Mass()),inertia_tensor(body.Inertia_Tensor())/*,has_infinite_inertia(body.Has_Infinite_Inertia())*/ {};

    VECTOR<T,0> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,1> > orientation,const VECTOR<T,0> angular_velocity) const
    {return VECTOR<T,0>();}
    VECTOR<T,0> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,1> > orientation,const VECTOR<T,0> angular_velocity) const
    {return inertia_tensor.Inverse_Times(angular_velocity);}
    VECTOR<T,1> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,2> >& orientation,const VECTOR<T,1> angular_velocity) const
    {return inertia_tensor*angular_velocity;}
    VECTOR<T,1> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,2> >& orientation,const VECTOR<T,1> angular_momentum) const
    {return inertia_tensor.Inverse_Times(angular_momentum);}
    VECTOR<T,3> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,3> >& orientation,const VECTOR<T,3>& angular_velocity) const
    {return orientation.Rotate(inertia_tensor*orientation.Inverse_Rotate(angular_velocity));}
    VECTOR<T,3> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,3> >& orientation,const VECTOR<T,3>& angular_momentum) const
    {return orientation.Rotate(inertia_tensor.Inverse_Times(orientation.Inverse_Rotate(angular_momentum)));}
    TWIST<TV> World_Space_Mass_Inverse_Times(const TWIST<TV>& momentum) const
    {
        TWIST<TV> velocity;
        velocity.linear=momentum.linear*(1.0/mass);
        velocity.angular=World_Space_Inertia_Tensor_Inverse_Times(frame.r,momentum.angular);
        return velocity;
    }
};

template<class TV>
typename SOLVE_CONTACT::CONTACT<TV> Build_Contact(const int id0,const int id1,
    const bool infinite_inertia0,const bool infinite_inertia1,
    const PRECONDITIONER_RIGID_BODY<TV>& body0,PRECONDITIONER_RIGID_BODY<TV>& body1,
    const TV& input_location,const TV& input_normal,const typename TV::SCALAR input_distance,const typename TV::SCALAR dt);

template<class TV>
PRECONDITIONER_RIGID_BODY<TV> Merge_Bodies(const PRECONDITIONER_RIGID_BODY<TV>& body0,const PRECONDITIONER_RIGID_BODY<TV>& body1);

template<class TV>
class PRECONDITIONER
{
public:
    typedef typename TV::SCALAR T;
    static const int d=TV::m;

    //needed for working with a rigid body collection with inactive bodies
    ARRAY<VECTOR<int,2> > coarse_to_fine_bodies_map;
    ARRAY<int> fine_to_coarse_bodies_map;
    ARRAY<int> fine_to_coarse_contacts_map;
    ARRAY<int> coarse_to_fine_contacts_map;
    ARRAY<PRECONDITIONER_RIGID_BODY<TV> > bodies;
    ARRAY<SOLVE_CONTACT::CONTACT<TV> > contacts;
    ARRAY<T> lambda_normal;
    ARRAY<VECTOR<T,d-1> > lambda_tangent;
    ARRAY<TWIST<TV> > velocities;
    ARRAY<bool> has_infinite_inertia;

    void Build_Fine_Problem(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,
        ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts_input);
    void Update_Fine_Velocities(RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    bool Build_Coarse_Problem(const ARRAY<PRECONDITIONER_RIGID_BODY<TV> >& fine_bodies,
        const ARRAY<bool>& fine_has_infinite_inertia,const ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts);
    void Update_Coarse_Problem(ARRAY<TWIST<TV> >& fine_velocities,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts);
    T Residual(ARRAY<TWIST<TV> >& velocities,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts,ARRAY<T>& lambda_normal);
    void Update_Coarse_Body_Velocities_Forces(
        ARRAY<PRECONDITIONER_RIGID_BODY<TV> >& fine_bodies,
        ARRAY<bool>& fine_has_infinite_inertia,
        ARRAY<TWIST<TV> >& fine_velocities,
        ARRAY<T>& fine_lambda_normal,
        ARRAY<VECTOR<T,d-1> >& fine_lambda_tangent,
        ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts);
    void Solve_Subset(ARRAY<int>& indices,T tolerance,int iteration_maximum);
    void Solve(T tolerance,int iteration_maximum,bool recursive=true);
};

}
}

#endif
