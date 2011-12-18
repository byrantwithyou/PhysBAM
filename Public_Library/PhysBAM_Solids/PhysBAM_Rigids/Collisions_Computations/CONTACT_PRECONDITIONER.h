//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONTACT_PRECONDITIONER
//##################################################################### 
#ifndef __CONTACT_PRECONDITIONER__
#define __CONTACT_PRECONDITIONER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/PROJECTED_GAUSS_SEIDEL.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SOLVE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>

namespace PhysBAM{

namespace CONTACT_PRECONDITIONER
{

template<class TV>
class PRECONDITIONER_RIGID_BODY
{
public:
    typedef typename TV::SCALAR T;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_INERTIA_TENSOR;
    FRAME<TV> frame;
    //TWIST<TV> twist;
    T mass;
    T_INERTIA_TENSOR inertia_tensor;
    
    //bool has_infinite_inertia;

    PRECONDITIONER_RIGID_BODY() {}
    PRECONDITIONER_RIGID_BODY(RIGID_BODY<TV>& body):
        frame(body.Frame()),/*twist(body.Twist()),*/mass(body.Mass()),inertia_tensor(body.Inertia_Tensor())/*,has_infinite_inertia(body.Has_Infinite_Inertia())*/ {};

    VECTOR<T,0> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,1> > orientation,const VECTOR<T,0> angular_velocity) const
    {return VECTOR<T,0>();}
    VECTOR<T,0> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,1> > orientation,const VECTOR<T,0> angular_velocity) const
    {return inertia_tensor.Solve_Linear_System(angular_velocity);}
    VECTOR<T,1> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,2> >& orientation,const VECTOR<T,1> angular_velocity) const
    {return inertia_tensor*angular_velocity;}
    VECTOR<T,1> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,2> >& orientation,const VECTOR<T,1> angular_momentum) const
    {return inertia_tensor.Solve_Linear_System(angular_momentum);}
    VECTOR<T,3> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,3> >& orientation,const VECTOR<T,3>& angular_velocity) const
    {return orientation.Rotate(inertia_tensor*orientation.Inverse_Rotate(angular_velocity));}
    VECTOR<T,3> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,3> >& orientation,const VECTOR<T,3>& angular_momentum) const
    {return orientation.Rotate(inertia_tensor.Solve_Linear_System(orientation.Inverse_Rotate(angular_momentum)));}
    TWIST<TV> World_Space_Mass_Inverse_Times(const TWIST<TV>& momentum) const
    {
        TWIST<TV> velocity;
        velocity.linear=momentum.linear*(1.0/mass);
        velocity.angular=World_Space_Inertia_Tensor_Inverse_Times(frame.r,momentum.angular);
        return velocity;
    }
};

template<class TV>
typename SOLVE_CONTACT::CONTACT<TV> Build_Contact(
    const int id0,const int id1,
    const bool infinite_inertia0,const bool infinite_inertia1,
    const PRECONDITIONER_RIGID_BODY<TV>& body0,PRECONDITIONER_RIGID_BODY<TV>& body1,
    const TV& input_location,const TV& input_normal,const typename TV::SCALAR input_distance,const typename TV::SCALAR dt)
{
    typename SOLVE_CONTACT::CONTACT<TV> contact;
    contact.id(1)=id0;
    contact.id(2)=id1;
    contact.location=input_location;
    contact.normal=input_normal;
    contact.distance=input_distance;
    contact.coefficient_of_friction=0;
    
    TV r0=contact.location-body0.frame.t;
    TV r1=contact.location-body1.frame.t;

    //LOG::cout << "build contact " << body0.frame.t << " " << body1.frame.t << " " << r0 << " " << r1 << std::endl;
    
    contact.normal_constraint(1).linear=-contact.normal;
    contact.normal_constraint(1).angular=TV::Cross_Product(contact.normal,r0);
    
    contact.normal_constraint(2).linear=contact.normal;
    contact.normal_constraint(2).angular=TV::Cross_Product(contact.normal,r1);
    
    contact.inverse_mass_times_normal_constraint(1)=body0.World_Space_Mass_Inverse_Times(contact.normal_constraint(1));
    contact.inverse_mass_times_normal_constraint(2)=body1.World_Space_Mass_Inverse_Times(contact.normal_constraint(2));
    
    contact.normal_diagonal=0;
    if(!infinite_inertia0)
        contact.normal_diagonal+=TWIST<TV>::Dot_Product(contact.normal_constraint(1),contact.inverse_mass_times_normal_constraint(1));
    if(!infinite_inertia1)
        contact.normal_diagonal+=TWIST<TV>::Dot_Product(contact.normal_constraint(2),contact.inverse_mass_times_normal_constraint(2));
    contact.normal_relative_velocity=-contact.distance/dt;

    return contact;
}

template<class T>
typename RIGID_BODY_POLICY<VECTOR<T,1> >::WORLD_SPACE_INERTIA_TENSOR Transform_Inertia_Tensor(const typename RIGID_BODY_POLICY<VECTOR<T,1> >::WORLD_SPACE_INERTIA_TENSOR& inertia_tensor,const T mass,const VECTOR<T,1>& t,const ROTATION<VECTOR<T,1> >& r)
{
    return inertia_tensor;
}

template<class T>
typename RIGID_BODY_POLICY<VECTOR<T,2> >::WORLD_SPACE_INERTIA_TENSOR Transform_Inertia_Tensor(const typename RIGID_BODY_POLICY<VECTOR<T,2> >::WORLD_SPACE_INERTIA_TENSOR& inertia_tensor,const T mass,const VECTOR<T,2>& t,const ROTATION<VECTOR<T,2> >& r)
{
    return inertia_tensor+mass*t.Magnitude_Squared();
}

template<class T>
typename RIGID_BODY_POLICY<VECTOR<T,3> >::WORLD_SPACE_INERTIA_TENSOR Transform_Inertia_Tensor(const typename RIGID_BODY_POLICY<VECTOR<T,3> >::WORLD_SPACE_INERTIA_TENSOR& inertia_tensor,const T mass,const VECTOR<T,3>& t,const ROTATION<VECTOR<T,3> >& r)
{
    typedef typename RIGID_BODY_POLICY<VECTOR<T,3> >::WORLD_SPACE_INERTIA_TENSOR T_INERTIA_TENSOR;
    T_INERTIA_TENSOR transformed_inertia_tensor;
    transformed_inertia_tensor.From_Matrix(r.Rotation_Matrix()*inertia_tensor*r.Inverse().Rotation_Matrix() + mass*(VECTOR<T,3>::Dot_Product(t,t)*T_INERTIA_TENSOR::Identity_Matrix()-T_INERTIA_TENSOR::Outer_Product(t)));
    //LOG::cout << (VECTOR<T,3>::Dot_Product(t,t)*T_INERTIA_TENSOR::Identity_Matrix()+T_INERTIA_TENSOR::Outer_Product(t)) << std::endl;
    return transformed_inertia_tensor;
}

template<class TV>
PRECONDITIONER_RIGID_BODY<TV> Merge_Bodies(const PRECONDITIONER_RIGID_BODY<TV>& body0,const PRECONDITIONER_RIGID_BODY<TV>& body1)
{
    typedef typename TV::SCALAR T;
    typedef typename RIGID_BODY_POLICY<VECTOR<T,TV::dimension> >::WORLD_SPACE_INERTIA_TENSOR T_INERTIA_TENSOR;
    PRECONDITIONER_RIGID_BODY<TV> merged_body;

    T mass=body0.mass+body1.mass;
    TV center=(body0.frame.t*body0.mass+body1.frame.t*body1.mass)*(1.0/mass);
    
    TV translation0=center-body0.frame.t;
    TV translation1=center-body1.frame.t;

    T_INERTIA_TENSOR inertia_tensor=
        Transform_Inertia_Tensor(body0.inertia_tensor,body0.mass,translation0,body0.frame.r) + 
        Transform_Inertia_Tensor(body1.inertia_tensor,body1.mass,translation1,body1.frame.r);

    //LOG::cout << "body0 " << body0.mass << " " << translation0 << " " << body0.frame.r.Rotation_Matrix() << std::endl;
    //LOG::cout << (VECTOR<T,3>::Dot_Product(translation0,translation0)*T_INERTIA_TENSOR::Identity_Matrix()+T_INERTIA_TENSOR::Outer_Product(translation0)) << std::endl;

    //LOG::cout << "inertia tensor 0 " << body0.inertia_tensor << " transformed " << Transform_Inertia_Tensor(body0.inertia_tensor,body0.mass,translation0,body0.frame.r) << std::endl;
    //LOG::cout << "inertia tensor 1 " << body1.inertia_tensor << " transformed " << Transform_Inertia_Tensor(body1.inertia_tensor,body1.mass,translation1,body1.frame.r) << std::endl;
    //LOG::cout << "inertia tensor " << inertia_tensor << std::endl;

    merged_body.frame.t=center;
    merged_body.frame.r=ROTATION<TV>();
    merged_body.mass=mass;
    merged_body.inertia_tensor=inertia_tensor;

    return merged_body;
}

template<class TV>
class PRECONDITIONER
{
public:
    typedef typename TV::SCALAR T;
    static const int d=TV::dimension;

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

    void Build_Fine_Problem(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts_input)
    {
        fine_to_coarse_bodies_map.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size());
        for(int i=1;i<=fine_to_coarse_bodies_map.m;i++)
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
        for(int i=1;i<=contacts_input.m;i++)
        {
            SOLVE_CONTACT::CONTACT<TV> contact=contacts_input(i);
            if(fine_to_coarse_bodies_map(contact.id(1)) && fine_to_coarse_bodies_map(contact.id(2)))
            {
                contact.id(1)=fine_to_coarse_bodies_map(contact.id(1));
                contact.id(2)=fine_to_coarse_bodies_map(contact.id(2));
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

    void Update_Fine_Velocities(RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    {
        for(int i=1;i<=coarse_to_fine_bodies_map.m;i++)
            rigid_body_collection.Rigid_Body(coarse_to_fine_bodies_map(i)(1)).Twist()=velocities(i);
        //LOG::cout << "lambda_normal " << lambda_normal << std::endl;
        LOG::cout << "velocities " << velocities << std::endl;
    }

    bool Build_Coarse_Problem(const ARRAY<PRECONDITIONER_RIGID_BODY<TV> >& fine_bodies,const ARRAY<bool>& fine_has_infinite_inertia,const ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts)
    {
        //for(int i=1;i<=fine_bodies.m;i++)
        //    LOG::cout << "fine body inertia " << fine_bodies(i).mass.mass << std::endl;

        bool merged=false;

        fine_to_coarse_bodies_map.Resize(fine_bodies.m);
        for(int i=1;i<=fine_contacts.m;i++)
        {
            SOLVE_CONTACT::CONTACT<TV> contact=fine_contacts(i);
            int id0=contact.id(1);
            int id1=contact.id(2);
            if(!fine_to_coarse_bodies_map(id0) && !fine_has_infinite_inertia(id0) && !fine_to_coarse_bodies_map(id1) && !fine_has_infinite_inertia(id1))
            {
                bodies.Append(Merge_Bodies(fine_bodies(id0),fine_bodies(id1)));
                velocities.Append(TWIST<TV>());
                has_infinite_inertia.Append(false);
                coarse_to_fine_bodies_map.Append(contact.id);
                fine_to_coarse_bodies_map(id0)=bodies.m;
                fine_to_coarse_bodies_map(id1)=bodies.m;
                merged=true;
//    LOG::cout << "body inertia " << bodies(bodies.m).mass.mass << std::endl;
            }
        }
        for(int i=1;i<=fine_bodies.m;i++)
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
        for(int i=1;i<=fine_contacts.m;i++)
        {
            const SOLVE_CONTACT::CONTACT<TV>& old_contact=fine_contacts(i);
            int id0=fine_to_coarse_bodies_map(old_contact.id(1));
            int id1=fine_to_coarse_bodies_map(old_contact.id(2));
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

        //LOG::cout << "fine_to_coarse_bodies_map " << fine_to_coarse_bodies_map << std::endl;
        //LOG::cout << "has_infinite_inertia " << has_infinite_inertia << std::endl;

        return merged;
    }

    void Update_Coarse_Problem(ARRAY<TWIST<TV> >& fine_velocities,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts)
    {
        lambda_normal.Fill(0);
        //lambda_tangent.Fill(VECTOR<T,d-1>((T)0));
        for(int i=1;i<=fine_contacts.m;i++)
        {
            int id=fine_to_coarse_contacts_map(i);
            if(id)
            {
                SOLVE_CONTACT::CONTACT<TV>& contact=contacts(id);
                SOLVE_CONTACT::CONTACT<TV>& fine_contact=fine_contacts(i);
                T current_residual=TWIST<TV>::Dot_Product(fine_contact.normal_constraint(1),fine_velocities(fine_contact.id(1)))+TWIST<TV>::Dot_Product(fine_contact.normal_constraint(2),fine_velocities(fine_contact.id(2)))-fine_contact.normal_relative_velocity;
                //LOG::cout << "coarse_contact " << fine_contact.normal_relative_velocity << " " << fine_velocities(fine_contact.id(1)) << " " << fine_velocities(fine_contact.id(2)) << std::endl;
                contact.normal_relative_velocity=-current_residual;
            }
        }
    }

    T Residual(ARRAY<TWIST<TV> >& velocities,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts,ARRAY<T>& lambda_normal)
    {
        T maximum_residual=0;
        for(int i=1;i<=contacts.m;i++)
        {
            SOLVE_CONTACT::CONTACT<TV>& contact=contacts(i);
            T violation=TWIST<TV>::Dot_Product(contact.normal_constraint(1),velocities(contact.id(1)))+TWIST<TV>::Dot_Product(contact.normal_constraint(2),velocities(contact.id(2)))-contact.normal_relative_velocity;
            T complementarity_condition=violation*lambda_normal(i);
            T residual=max(-violation,(T)fabs(complementarity_condition));
            maximum_residual=max(maximum_residual,residual);
        }
        return maximum_residual;
    }

    void Update_Coarse_Body_Velocities_Forces(
        ARRAY<PRECONDITIONER_RIGID_BODY<TV> >& fine_bodies,
        ARRAY<bool>& fine_has_infinite_inertia,
        ARRAY<TWIST<TV> >& fine_velocities,
        ARRAY<T>& fine_lambda_normal,
        ARRAY<VECTOR<T,d-1> >& fine_lambda_tangent,
        ARRAY<SOLVE_CONTACT::CONTACT<TV> >& fine_contacts)
    {
        LOG::SCOPE scope("CONTACT_PRECONDITIONER::Update_Coarse_Body_Velocities_Forces");
        
        //LOG::cout << "fine_lambda_normal 0 - " << fine_lambda_normal << std::endl;
        //LOG::cout << "fine_velocities " << fine_velocities << std::endl;

        //LOG::cout << "fine_to_coarse_contacts_map " << fine_to_coarse_contacts_map << std::endl;

        /*for(int i=1;i<=contacts.m;i++)
        {
            SOLVE_CONTACT::CONTACT<TV>& contact=contacts(i);
            T residual=TWIST<TV>::Dot_Product(contact.normal_constraint(1),velocities(contact.id(1)))+TWIST<TV>::Dot_Product(contact.normal_constraint(2),velocities(contact.id(2)))-contact.normal_relative_velocity;
            LOG::cout << "coarse_residual " << residual << " " << contact.normal_relative_velocity << std::endl;
            }*/
        
        //LOG::cout << "coarse_velocities " << velocities << std::endl;
        //LOG::cout << "coarse_lambda_normal " << lambda_normal << std::endl;

        //LOG::cout << "fine_lambda_normal 1 - " << fine_lambda_normal << std::endl;

        //HASHTABLE<VECTOR<int,2>,VECTOR<TWIST<TV>,2> > pair_wrenches;

        for(int i=1;i<=fine_contacts.m;i++)
        {
            int id=fine_to_coarse_contacts_map(i);
            if(id)
            {
                SOLVE_CONTACT::CONTACT<TV>& fine_contact=fine_contacts(i);
                //SOLVE_CONTACT::CONTACT<TV>& coarse_contact=contacts(id);
                //LOG::cout << "contact " << i << " " << id << " " << fine_contact.normal_constraint(2) << " " << coarse_contact.normal_constraint(2) << std::endl;

                fine_lambda_normal(i)+=lambda_normal(id);
                //fine_lambda_tangent+=lambda_tangent(id);
                if(!fine_has_infinite_inertia(fine_contact.id(1)))
                    fine_velocities(fine_contact.id(1))+=fine_contact.inverse_mass_times_normal_constraint(1)*lambda_normal(id);
                if(!fine_has_infinite_inertia(fine_contact.id(2)))
                    fine_velocities(fine_contact.id(2))+=fine_contact.inverse_mass_times_normal_constraint(2)*lambda_normal(id);
                /*for(int i=1;i<d;i++)
                {
                    fine_velocities(fine_contact.id(1))+=fine_contact.inverse_mass_times_tangent_constraint(1)(i)*lambda_tangent(id)(i);
                    fine_velocities(fine_contact.id(2))+=fine_contact.inverse_mass_times_tangent_constraint(2)(i)*lambda_tangent(id)(i);
                }*/

                //pair_wrenches.Get_Or_Insert(fine_contact.id.Sorted())+=VECTOR<TWIST<TV>,2>(fine_contact.inverse_mass_times_normal_constraint(1)*lambda_normal(id),fine_contact.inverse_mass_times_normal_constraint(2)*lambda_normal(id));
            }
        }

        /*for(HASHTABLE_ITERATOR<VECTOR<int,2>,VECTOR<TWIST<TV>,2> > it(pair_wrenches);it.Valid();it.Next())
        {
            LOG::cout << "wrench " << it.Key() << " " << it.Data() << std::endl;
            }*/
        
        //LOG::cout << "fine_lambda_normal 2 - " << fine_lambda_normal << std::endl;
        //LOG::cout << "fine_velocities " << fine_velocities << std::endl;
    }

    void Solve_Subset(ARRAY<int>& indices,T tolerance,int iteration_maximum)
    {
      //HASHTABLE<VECTOR<int,2> > pairs;

        ARRAY<SOLVE_CONTACT::CONTACT<TV> > subset_contacts(indices.m);
        ARRAY<T> subset_lambda_normal(indices.m);
        ARRAY<VECTOR<T,d-1> > subset_lambda_tangent(indices.m);
        for(int i=1;i<=indices.m;i++)
        {
            int index=indices(i);
            subset_contacts(i)=contacts(index);
            subset_lambda_normal(i)=lambda_normal(index);
            subset_lambda_tangent(i)=lambda_tangent(index);
        //    pairs.Set(subset_contacts(i).id.Sorted());
        }

        //for(HASHTABLE_ITERATOR<VECTOR<int,2> > it(pairs);it.Valid();it.Next())
        //{
        //    LOG::cout << "subset pair " << it.Key() << std::endl;
        //}

        PROJECTED_GAUSS_SEIDEL::Solve(velocities,has_infinite_inertia,subset_contacts,subset_lambda_normal,subset_lambda_tangent,tolerance,iteration_maximum,false);

        for(int i=1;i<=indices.m;i++)
        {
            int index=indices(i);
            lambda_normal(index)=subset_lambda_normal(i);
            lambda_tangent(index)=subset_lambda_tangent(i);
        }
    }

    void Solve(T tolerance,int iteration_maximum,bool recursive=true)
    {
        for(int i=1;i<=contacts.m;i++)
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
                for(int i=1;i<=contacts.m;i++)
                    if(!preconditioner.fine_to_coarse_contacts_map(i))
                        removed_indices.Append(i);
                
                //PGS on removed contacts fixed iterations
                //Solve_Subset(removed_indices,(T)1e-5,10000);
                
                //Solve coarse problem with updated residuals
                preconditioner.Update_Coarse_Problem(velocities,contacts);
                preconditioner.Solve(tolerance,0,true);

                //LOG::cout << "pre coarse update velocities " << velocities << std::endl;
                
                preconditioner.Update_Coarse_Body_Velocities_Forces(bodies,has_infinite_inertia,velocities,lambda_normal,lambda_tangent,contacts);
                
                //LOG::cout << "pre subset solve velocities " << velocities << std::endl;

                //PGS on removed contacts fixed iterations
                Solve_Subset(removed_indices,1e-8,0);
                Solve(tolerance,10,false);

                //LOG::cout << "lambda_normal 1 - " << lambda_normal << std::endl;
                
                //PGS on entire system fixed iterations
                //Solve(tolerance,10,false);
            }
            else
                Solve(tolerance,0,false);
        }

        //LOG::cout << "velocities " << velocities << std::endl;
        LOG::cout << "residual 0 " << Residual(velocities,contacts,lambda_normal) << std::endl;

        /*for(int i=1;i<=bodies.m;i++)
        {
            LOG::cout << "mass " << i << " " << bodies(i).mass << std::endl;
            }*/

//print out force sums per interaction pair
        /*HASHTABLE<VECTOR<int,2>,T> pair_forces;
        HASHTABLE<VECTOR<int,2>,VECTOR<TWIST<TV>,2> > pair_wrenches;
        for(int i=1;i<=contacts.m;i++)
        {
            //LOG::cout << contacts(i).inverse_mass_times_normal_constraint << std::endl;
            VECTOR<int,2> pair=contacts(i).id.Sorted();
            pair_forces.Get_Or_Insert(pair)+=lambda_normal(i);
            pair_wrenches.Get_Or_Insert(pair)+=VECTOR<TWIST<TV>,2>(contacts(i).inverse_mass_times_normal_constraint(1)*lambda_normal(i),contacts(i).inverse_mass_times_normal_constraint(2)*lambda_normal(i));
        }

        for(HASHTABLE_ITERATOR<VECTOR<int,2>,T> it(pair_forces);it.Valid();it.Next())
        {
            LOG::cout << "force " << it.Key() << " " << it.Data() << std::endl;
        }
        for(HASHTABLE_ITERATOR<VECTOR<int,2>,VECTOR<TWIST<TV>,2> > it(pair_wrenches);it.Valid();it.Next())
        {
            LOG::cout << "wrench " << it.Key() << " " << it.Data() << std::endl;
            }*/
    }
};

}
}

#endif
