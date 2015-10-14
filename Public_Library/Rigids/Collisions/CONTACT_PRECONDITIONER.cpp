//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Collisions/CONTACT_PRECONDITIONER.h>
namespace PhysBAM{
namespace CONTACT_PRECONDITIONER
{

//#####################################################################
// Function Build_Contact
//#####################################################################
template<class TV> typename SOLVE_CONTACT::CONTACT<TV>
Build_Contact(const int id0,const int id1,const bool infinite_inertia0,const bool infinite_inertia1,
    const PRECONDITIONER_RIGID_BODY<TV>& body0,PRECONDITIONER_RIGID_BODY<TV>& body1,
    const TV& input_location,const TV& input_normal,const typename TV::SCALAR input_distance,const typename TV::SCALAR dt)
{
    typename SOLVE_CONTACT::CONTACT<TV> contact;
    contact.id(0)=id0;
    contact.id(1)=id1;
    contact.location=input_location;
    contact.normal=input_normal;
    contact.distance=input_distance;
    contact.coefficient_of_friction=0;
    
    TV r0=contact.location-body0.frame.t;
    TV r1=contact.location-body1.frame.t;
    
    contact.normal_constraint(0).linear=-contact.normal;
    contact.normal_constraint(0).angular=TV::Cross_Product(contact.normal,r0);
    
    contact.normal_constraint(1).linear=contact.normal;
    contact.normal_constraint(1).angular=TV::Cross_Product(contact.normal,r1);
    
    contact.inverse_mass_times_normal_constraint(0)=body0.World_Space_Mass_Inverse_Times(contact.normal_constraint(0));
    contact.inverse_mass_times_normal_constraint(1)=body1.World_Space_Mass_Inverse_Times(contact.normal_constraint(1));
    
    contact.normal_diagonal=0;
    if(!infinite_inertia0)
        contact.normal_diagonal+=TWIST<TV>::Dot_Product(contact.normal_constraint(0),contact.inverse_mass_times_normal_constraint(0));
    if(!infinite_inertia1)
        contact.normal_diagonal+=TWIST<TV>::Dot_Product(contact.normal_constraint(1),contact.inverse_mass_times_normal_constraint(1));
    contact.normal_relative_velocity=-contact.distance/dt;

    return contact;
}
//#####################################################################
// Function Transform_Inertia_Tensor
//#####################################################################
template<class T> static SYMMETRIC_MATRIX<T,0>
Transform_Inertia_Tensor(const SYMMETRIC_MATRIX<T,0>& inertia_tensor,const T mass,const VECTOR<T,1>& t,
    const ROTATION<VECTOR<T,1> >& r)
{
    return inertia_tensor;
}
//#####################################################################
// Function Transform_Inertia_Tensor
//#####################################################################
template<class T> static SYMMETRIC_MATRIX<T,1>
Transform_Inertia_Tensor(const SYMMETRIC_MATRIX<T,1>& inertia_tensor,const T mass,const VECTOR<T,2>& t,
    const ROTATION<VECTOR<T,2> >& r)
{
    return inertia_tensor+mass*t.Magnitude_Squared();
}
//#####################################################################
// Function Transform_Inertia_Tensor
//#####################################################################
template<class T> static SYMMETRIC_MATRIX<T,3>
Transform_Inertia_Tensor(const SYMMETRIC_MATRIX<T,3>& inertia_tensor,const T mass,const VECTOR<T,3>& t,
    const ROTATION<VECTOR<T,3> >& r)
{
    SYMMETRIC_MATRIX<T,3> transformed_inertia_tensor;
    transformed_inertia_tensor.From_Matrix(r.Rotation_Matrix()*inertia_tensor*r.Inverse().Rotation_Matrix() + mass*(VECTOR<T,3>::Dot_Product(t,t)*SYMMETRIC_MATRIX<T,3>::Identity_Matrix()-SYMMETRIC_MATRIX<T,3>::Outer_Product(t)));
    return transformed_inertia_tensor;
}
//#####################################################################
// Function Merge_Bodies
//#####################################################################
template<class TV> PRECONDITIONER_RIGID_BODY<TV>
Merge_Bodies(const PRECONDITIONER_RIGID_BODY<TV>& body0,const PRECONDITIONER_RIGID_BODY<TV>& body1)
{
    typedef typename TV::SCALAR T;
    PRECONDITIONER_RIGID_BODY<TV> merged_body;

    T mass=body0.mass+body1.mass;
    TV center=(body0.frame.t*body0.mass+body1.frame.t*body1.mass)*(1.0/mass);
    
    TV translation0=center-body0.frame.t;
    TV translation1=center-body1.frame.t;

    SYMMETRIC_MATRIX<T,TV::SPIN::m> inertia_tensor=
        Transform_Inertia_Tensor(body0.inertia_tensor,body0.mass,translation0,body0.frame.r) + 
        Transform_Inertia_Tensor(body1.inertia_tensor,body1.mass,translation1,body1.frame.r);

    merged_body.frame.t=center;
    merged_body.frame.r=ROTATION<TV>();
    merged_body.mass=mass;
    merged_body.inertia_tensor=inertia_tensor;

    return merged_body;
}
//#####################################################################
// Function Build_Fine_Problem
//#####################################################################
template<class TV> void PRECONDITIONER<TV>::
Build_Fine_Problem(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts_input)
{
    fine_to_coarse_bodies_map.Resize(rigid_body_collection.rigid_body_particles.Size());
    for(int i=0;i<fine_to_coarse_bodies_map.m;i++)
        if(rigid_body_collection.Is_Active(i))
        {
            RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
            bodies.Append(body);
            velocities.Append(body.Twist());
            has_infinite_inertia.Append(body.Has_Infinite_Inertia());
            coarse_to_fine_bodies_map.Append(VECTOR<int,2>(i,0));
            fine_to_coarse_bodies_map(i)=bodies.m;
        }
        else
        {
            fine_to_coarse_bodies_map(i)=0;
        }
        
    fine_to_coarse_contacts_map.Resize(contacts_input.m);
    for(int i=0;i<contacts_input.m;i++)
    {
        SOLVE_CONTACT::CONTACT<TV> contact=contacts_input(i);
        if(fine_to_coarse_bodies_map(contact.id(0)) && fine_to_coarse_bodies_map(contact.id(1)))
        {
            contact.id(0)=fine_to_coarse_bodies_map(contact.id(0));
            contact.id(1)=fine_to_coarse_bodies_map(contact.id(1));
            contact.coefficient_of_friction=0;
            contacts.Append(contact);
            coarse_to_fine_contacts_map.Append(i);
            fine_to_coarse_contacts_map(i)=contacts.m;
        }
        else
        {
            fine_to_coarse_contacts_map(i)=0;
        }
    }
    lambda_normal.Resize(contacts.m);
    lambda_tangent.Resize(contacts.m);
}
//#####################################################################
// Function Update_Fine_Velocities
//#####################################################################
template<class TV> void PRECONDITIONER<TV>::
Update_Fine_Velocities(RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
{
    for(int i=0;i<coarse_to_fine_bodies_map.m;i++)
        rigid_body_collection.Rigid_Body(coarse_to_fine_bodies_map(i)(0)).Twist()=velocities(i);
    LOG::cout << "velocities " << velocities << std::endl;
}
//#####################################################################
// Function Build_Coarse_Problem
//#####################################################################
template<class TV> bool PRECONDITIONER<TV>::
Build_Coarse_Problem(const ARRAY<PRECONDITIONER_RIGID_BODY<TV> >& fine_bodies,
    const ARRAY<bool>& fine_has_infinite_inertia,const ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts)
{
    bool merged=false;

    fine_to_coarse_bodies_map.Resize(fine_bodies.m);
    for(int i=0;i<fine_contacts.m;i++)
    {
        SOLVE_CONTACT::CONTACT<TV> contact=fine_contacts(i);
        int id0=contact.id(0);
        int id1=contact.id(1);
        if(!fine_to_coarse_bodies_map(id0) && !fine_has_infinite_inertia(id0) && !fine_to_coarse_bodies_map(id1) && !fine_has_infinite_inertia(id1))
        {
            bodies.Append(Merge_Bodies(fine_bodies(id0),fine_bodies(id1)));
            velocities.Append(TWIST<TV>());
            has_infinite_inertia.Append(false);
            coarse_to_fine_bodies_map.Append(contact.id);
            fine_to_coarse_bodies_map(id0)=bodies.m;
            fine_to_coarse_bodies_map(id1)=bodies.m;
            merged=true;
        }
    }
    for(int i=0;i<fine_bodies.m;i++)
    {
        if(!fine_to_coarse_bodies_map(i))
        {
            bodies.Append(fine_bodies(i));
            velocities.Append(TWIST<TV>());
            has_infinite_inertia.Append(fine_has_infinite_inertia(i));
            coarse_to_fine_bodies_map.Append(VECTOR<int,2>(i,0));
            fine_to_coarse_bodies_map(i)=bodies.m;
        }
    }

    fine_to_coarse_contacts_map.Resize(fine_contacts.m);
    for(int i=0;i<fine_contacts.m;i++)
    {
        const SOLVE_CONTACT::CONTACT<TV>& old_contact=fine_contacts(i);
        int id0=fine_to_coarse_bodies_map(old_contact.id(0));
        int id1=fine_to_coarse_bodies_map(old_contact.id(1));
        if(id0!=id1)
        {
            contacts.Append(Build_Contact(id0,id1,has_infinite_inertia(id0),has_infinite_inertia(id1),bodies(id0),bodies(id1),old_contact.location,old_contact.normal,0,1));
            coarse_to_fine_contacts_map.Append(i);
            fine_to_coarse_contacts_map(i)=contacts.m;
        }
        else
        {
            fine_to_coarse_contacts_map(i)=0;
        }
    }
    lambda_normal.Resize(contacts.m);
    lambda_tangent.Resize(contacts.m);

    return merged;
}
//#####################################################################
// Function Update_Coarse_Problem
//#####################################################################
template<class TV> void PRECONDITIONER<TV>::
Update_Coarse_Problem(ARRAY<TWIST<TV> >& fine_velocities,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts)
{
    lambda_normal.Fill(0);
    //lambda_tangent.Fill(VECTOR<T,d-1>((T)0));
    for(int i=0;i<fine_contacts.m;i++)
    {
        int id=fine_to_coarse_contacts_map(i);
        if(id>=0)
        {
            SOLVE_CONTACT::CONTACT<TV>& contact=contacts(id);
            SOLVE_CONTACT::CONTACT<TV>& fine_contact=fine_contacts(i);
            T current_residual=TWIST<TV>::Dot_Product(fine_contact.normal_constraint(0),fine_velocities(fine_contact.id(0)))+TWIST<TV>::Dot_Product(fine_contact.normal_constraint(1),fine_velocities(fine_contact.id(1)))-fine_contact.normal_relative_velocity;
            contact.normal_relative_velocity=-current_residual;
        }
    }
}
//#####################################################################
// Function Residual
//#####################################################################
template<class TV> typename TV::SCALAR PRECONDITIONER<TV>::
Residual(ARRAY<TWIST<TV> >& velocities,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts,ARRAY<T>& lambda_normal)
{
    T maximum_residual=0;
    for(int i=0;i<contacts.m;i++)
    {
        SOLVE_CONTACT::CONTACT<TV>& contact=contacts(i);
        T violation=TWIST<TV>::Dot_Product(contact.normal_constraint(0),velocities(contact.id(0)))+TWIST<TV>::Dot_Product(contact.normal_constraint(1),velocities(contact.id(1)))-contact.normal_relative_velocity;
        T complementarity_condition=violation*lambda_normal(i);
        T residual=max(-violation,(T)fabs(complementarity_condition));
        maximum_residual=max(maximum_residual,residual);
    }
    return maximum_residual;
}
//#####################################################################
// Function Update_Coarse_Body_Velocities_Forces
//#####################################################################
template<class TV> void PRECONDITIONER<TV>::
Update_Coarse_Body_Velocities_Forces(ARRAY<PRECONDITIONER_RIGID_BODY<TV> >& fine_bodies,
    ARRAY<bool>& fine_has_infinite_inertia,
    ARRAY<TWIST<TV> >& fine_velocities,
    ARRAY<T>& fine_lambda_normal,
    ARRAY<VECTOR<T,d-1> >& fine_lambda_tangent,
    ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts)
{
    LOG::SCOPE scope("CONTACT_PRECONDITIONER::Update_Coarse_Body_Velocities_Forces");

    for(int i=0;i<fine_contacts.m;i++)
    {
        int id=fine_to_coarse_contacts_map(i);
        if(id>=0)
        {
            SOLVE_CONTACT::CONTACT<TV>& fine_contact=fine_contacts(i);
            fine_lambda_normal(i)+=lambda_normal(id);
            if(!fine_has_infinite_inertia(fine_contact.id(0)))
                fine_velocities(fine_contact.id(0))+=fine_contact.inverse_mass_times_normal_constraint(0)*lambda_normal(id);
            if(!fine_has_infinite_inertia(fine_contact.id(1)))
                fine_velocities(fine_contact.id(1))+=fine_contact.inverse_mass_times_normal_constraint(1)*lambda_normal(id);
        }
    }
}
//#####################################################################
// Function Solve_Subset
//#####################################################################
template<class TV> void PRECONDITIONER<TV>::
Solve_Subset(ARRAY<int>& indices,T tolerance,int iteration_maximum)
{
    ARRAY<SOLVE_CONTACT::CONTACT<TV> > subset_contacts(indices.m);
    ARRAY<T> subset_lambda_normal(indices.m);
    ARRAY<VECTOR<T,d-1> > subset_lambda_tangent(indices.m);
    for(int i=0;i<indices.m;i++)
    {
        int index=indices(i);
        subset_contacts(i)=contacts(index);
        subset_lambda_normal(i)=lambda_normal(index);
        subset_lambda_tangent(i)=lambda_tangent(index);
    }

    PROJECTED_GAUSS_SEIDEL::Solve(velocities,has_infinite_inertia,subset_contacts,subset_lambda_normal,subset_lambda_tangent,tolerance,iteration_maximum,false);

    for(int i=0;i<indices.m;i++)
    {
        int index=indices(i);
        lambda_normal(index)=subset_lambda_normal(i);
        lambda_tangent(index)=subset_lambda_tangent(i);
    }
}
//#####################################################################
// Function Solve
//#####################################################################
template<class TV> void PRECONDITIONER<TV>::
Solve(T tolerance,int iteration_maximum,bool recursive)
{
    for(int i=0;i<contacts.m;i++)
        LOG::cout << "contact " << i << " " << contacts(i).normal << std::endl;

    if(bodies.m<2 || contacts.m<1 || !recursive)
    {
        LOG::SCOPE scope("CONTACT_PRECONDITIEOR::Solve non-recursive");
        PROJECTED_GAUSS_SEIDEL::Solve(velocities,has_infinite_inertia,contacts,lambda_normal,lambda_tangent,(T)1e-8,iteration_maximum,false);
    }
    else
    {
        char scope_string[512];
        sprintf(scope_string,"CONTACT_PRECONDITIONER::Solve recursive %d %d",bodies.m,contacts.m);
        LOG::SCOPE scope(scope_string);
            
        //PGS on entire system fixed iterations
        //Solve(tolerance,10,false);

        //build and apply ML preconditioner
        PRECONDITIONER<TV> preconditioner;
        bool merged=preconditioner.Build_Coarse_Problem(bodies,has_infinite_inertia,contacts);
        LOG::cout << "merged " << merged << std::endl;
        if(merged)
        {
            ARRAY<int> removed_indices;
            for(int i=0;i<contacts.m;i++)
                if(!preconditioner.fine_to_coarse_contacts_map(i))
                    removed_indices.Append(i);
                
            //PGS on removed contacts fixed iterations
            //Solve_Subset(removed_indices,(T)1e-5,10000);
                
            //Solve coarse problem with updated residuals
            preconditioner.Update_Coarse_Problem(velocities,contacts);
            preconditioner.Solve(tolerance,0,true);

            preconditioner.Update_Coarse_Body_Velocities_Forces(bodies,has_infinite_inertia,velocities,lambda_normal,lambda_tangent,contacts);
                
            //PGS on removed contacts fixed iterations
            Solve_Subset(removed_indices,1e-8,0);
            Solve(tolerance,10,false);

            //PGS on entire system fixed iterations
            //Solve(tolerance,10,false);
        }
        else
            Solve(tolerance,0,false);
    }

    LOG::cout << "residual 0 " << Residual(velocities,contacts,lambda_normal) << std::endl;
}
}
}
