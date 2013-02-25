//#####################################################################
// Copyright 2010, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTED_GAUSS_SEIDEL
//##################################################################### 
#ifndef __PROJECTED_GAUSS_SEIDEL__
#define __PROJECTED_GAUSS_SEIDEL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SOLVE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>

namespace PhysBAM{;

namespace PROJECTED_GAUSS_SEIDEL
{
template<class T>
void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,ARRAY<T>& a,ARRAY<T>& x,T tolerance)
{
    assert(A.n==a.Size() && A.n==x.Size());

    int n=A.n;
    x.Fill(0);

    int iteration=0;

    T maximum_residual=1;
    while(maximum_residual>tolerance)
    {
        maximum_residual=0;
        for(int i=0;i<n;i++)
        {
            T row_sum=0;
            T diagonal=1;
            
            for(int j=A.offsets(i);j<A.offsets(i+1);j++)
            {
                int index=A.A(j).j;
                if(index==i)
                    diagonal=A.A(j).a;
                else
                    row_sum+=A.A(j).a*x(index);
            }

            T residual=row_sum+diagonal*x(i)-a(i);
            if(-residual>maximum_residual)
                maximum_residual=-residual;

            x(i)=(a(i)-row_sum)/diagonal;
            if(x(i)<0)
                x(i)=0;
        }

        iteration++;
    }

    LOG::cout << "PGS iterations = " << iteration << std::endl;
}

template<class T,int D>
void Multiply(SPARSE_MATRIX_FLAT_MXN<VECTOR<T,D> >& A,ARRAY<VECTOR<T,D> >& x,ARRAY<T>& result)
{
    assert(A.n==x.Size()&&A.m==result.Size());
    for(int i=0;i<A.m;i++)
    {
        T r=0;
        for(int j=A.offsets(i);j<A.offsets(i+1);j++)
            r+=VECTOR<T,D>::Dot_Product(A.A(j).a,x(A.A(j).j));
        result(i)=r;
    }
}

template<class T,int D>
    void Solve(ARRAY<MATRIX<T,D> >& A_block_diagonal,ARRAY<ARRAY<PAIR<int,VECTOR<T,D> > > >& C_block,ARRAY<VECTOR<T,D> >& a_block,ARRAY<T>& c,ARRAY<VECTOR<T,D> >& x,ARRAY<T>& lambda,T tolerance)
{
    int n_block_primal=A_block_diagonal.m;
    int n_dual=C_block.m;

    int iteration=0;

    for(int i=0;i<lambda.m;i++)
        lambda(i)=0;

    ARRAY<MATRIX<T,D> > A_block_diagonal_inverse(n_block_primal);
    for(int i=0;i<n_block_primal;i++)
    {
        A_block_diagonal(i).Cholesky_Inverse(A_block_diagonal_inverse(i));
        x(i)=A_block_diagonal_inverse(i)*a_block(i);
    }

    ARRAY<ARRAY<PAIR<int,VECTOR<T,D> > > > AiC=C_block;

    ARRAY<T> schur_complement_diagonal(n_dual);
    for(int i=0;i<n_dual;i++)
    {
        T diagonal=0;
        for(int j=0;j<C_block(i).m;j++)
        {
            VECTOR<T,D>& block=C_block(i)(j).y;
            VECTOR<T,D> Ai_block=A_block_diagonal_inverse(C_block(i)(j).x)*block;
            AiC(i)(j).y=Ai_block;
            diagonal+=VECTOR<T,D>::Dot_Product(block,Ai_block);
        }
        schur_complement_diagonal(i)=diagonal;
    }

    T maximum_residual=2*tolerance+1;

    while(maximum_residual>tolerance)
    {
        maximum_residual=0;
        for(int i=0;i<n_dual;i++)
        {
            T row_sum=0;
            for(int j=0;j<C_block(i).m;j++)
                row_sum+=VECTOR<T,D>::Dot_Product(C_block(i)(j).y,x(C_block(i)(j).x));

            T residual=c(i)-row_sum;
            if(residual>0 && residual>maximum_residual)
                maximum_residual=residual;

            T lambda_offset=residual/schur_complement_diagonal(i);

            T new_lambda=lambda(i)+lambda_offset;
            if(new_lambda<0)
            {
                lambda_offset=-lambda(i);
                new_lambda=0;
            }

            lambda(i)=new_lambda;

            for(int j=0;j<C_block(i).m;j++)
                x(C_block(i)(j).x)+=lambda_offset*AiC(i)(j).y;
        }

        //if(iteration%1000==0)
        //    LOG::cout << iteration << " maximum_residual = " << maximum_residual << std::endl;

        iteration++;
    }

    LOG::cout << "PGS iterations = " << iteration << std::endl;
    LOG::cout << "PGS maximum_residual = " << maximum_residual << std::endl;
}

//#####################################################################
// Function Solve_Projected_Gauss_Seidel
//#####################################################################
template<class TV>
bool Solve(
    ARRAY<TWIST<TV> >& velocities,
    ARRAY<bool>& has_infinite_inertia,
    ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts,
    ARRAY<typename TV::SCALAR>& lambda_normal,
    ARRAY<VECTOR<typename TV::SCALAR,TV::dimension-1> >& lambda_tangent,
    typename TV::SCALAR tolerance,
    int iteration_maximum,
    bool friction=false)
{
    LOG::SCOPE scope("PROJECTED_GAUSS_SEIDEL::Solve");

    typedef typename TV::SCALAR T;
    typedef TWIST<TV> T_TWIST;
    const int d=TV::dimension;

    int n_contacts=contacts.m;

    T maximum_residual=2*tolerance+1;

    /*for(int i=0;i<n_contacts;i++)
    {
        lambda_normal(i)=0;
        if(friction)
            for(int j=1;j<d;j++)
                lambda_tangent(i)(j)=0;
    }*/
    
    int iteration=0;

    while(maximum_residual>tolerance && (!iteration_maximum || iteration<iteration_maximum))
    {
        maximum_residual=0;
        for(int i=0;i<n_contacts;i++)
        {
            SOLVE_CONTACT::CONTACT<TV>& contact=contacts(i);
            
            T normal_violation=T_TWIST::Dot_Product(contact.normal_constraint(0),velocities(contact.id(0)))+T_TWIST::Dot_Product(contact.normal_constraint(1),velocities(contact.id(1)))-contact.normal_relative_velocity;

            T residual=max(-normal_violation,normal_violation*lambda_normal(i));
            if(residual>maximum_residual)
                maximum_residual=residual;
            
            T lambda_normal_delta=-normal_violation/contact.normal_diagonal;
            T lambda_normal_new=lambda_normal(i)+lambda_normal_delta;
            if(lambda_normal_new<0)
            {
                lambda_normal_new=0;
                lambda_normal_delta=-lambda_normal(i);
            }
            
            if(!has_infinite_inertia(contact.id(0)))
                velocities(contact.id(0))+=lambda_normal_delta*contact.inverse_mass_times_normal_constraint(0);
            if(!has_infinite_inertia(contact.id(1)))
                velocities(contact.id(1))+=lambda_normal_delta*contact.inverse_mass_times_normal_constraint(1);

            lambda_normal(i)=lambda_normal_new;

            if(friction)
                for(int j=1;j<d;j++)
                {
                    T tangent_violation=T_TWIST::Dot_Product(contact.tangent_constraint(0)(j),velocities(contact.id(0)))+T_TWIST::Dot_Product(contact.tangent_constraint(1)(j),velocities(contact.id(1)));
                    
                    T lambda_tangent_new=lambda_tangent(i)(j)-tangent_violation/contact.tangent_diagonal(j);
                    if(fabs(lambda_tangent_new)>contact.coefficient_of_friction*lambda_normal_new)
                    {
                        if(lambda_tangent_new>0)
                            lambda_tangent_new=contact.coefficient_of_friction*lambda_normal_new;
                        else
                            lambda_tangent_new=-contact.coefficient_of_friction*lambda_normal_new;
                    }
                    
                    if(!has_infinite_inertia(contact.id(0)))
                        velocities(contact.id(0))+=(lambda_tangent_new-lambda_tangent(i)(j))*contact.inverse_mass_times_tangent_constraint(0)(j);
                    if(!has_infinite_inertia(contact.id(1)))
                        velocities(contact.id(1))+=(lambda_tangent_new-lambda_tangent(i)(j))*contact.inverse_mass_times_tangent_constraint(1)(j);
                    
                    lambda_tangent(i)(j)=lambda_tangent_new;
                }
        }

        if(iteration%100==0)
          LOG::cout << "maximum_residual " << iteration << " " << maximum_residual << std::endl;
        
        iteration++;
    }

    LOG::cout << "iteration_maximum = " << iteration_maximum << std::endl;

    LOG::cout << "iterations = " << iteration << std::endl;
    LOG::cout << "maximum_residual = " << maximum_residual << std::endl;

    return true;
}
template<class TV>
bool Solve(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts,typename TV::SCALAR tolerance,int iteration_maximum)
{
    //LOG::SCOPE scope("Solve_Projected_Gauss_Seidel_Friction");

    typedef typename TV::SCALAR T;
    typedef TWIST<TV> T_TWIST;
    const int d=TV::dimension;

    int n_bodies=rigid_body_collection.rigid_body_particle.Size();
    int n_contacts=contacts.m;

    T maximum_residual=2*tolerance+1;

    ARRAY<T> lambda_normal(n_contacts);
    ARRAY<VECTOR<T,d-1> > lambda_tangent(n_contacts);
    for(int i=0;i<n_contacts;i++)
    {
        lambda_normal(i)=0;
        for(int j=1;j<d;j++)
            lambda_tangent(i)(j)=0;
    }
    ARRAY<TWIST<TV> > velocities(n_bodies);
    ARRAY<bool> has_infinite_inertia(n_bodies);
    
    for(int i=0;i<n_bodies;i++)
    {
        if(rigid_body_collection.Is_Active(i))
        {
            velocities(i)=rigid_body_collection.Rigid_Body(i).Twist();
            has_infinite_inertia(i)=rigid_body_collection.Rigid_Body(i).Has_Infinite_Inertia();
        }
    }

    int iteration=0;

    while(maximum_residual>tolerance && (!iteration_maximum || iteration<iteration_maximum))
    {
        maximum_residual=0;
        for(int i=0;i<n_contacts;i++)
        {
            SOLVE_CONTACT::CONTACT<TV>& contact=contacts(i);
            
            T normal_violation=T_TWIST::Dot_Product(contact.normal_constraint(0),velocities(contact.id(0)))+T_TWIST::Dot_Product(contact.normal_constraint(1),velocities(contact.id(1)))-contact.normal_relative_velocity;

            T residual=max(-normal_violation,normal_violation*lambda_normal(i));
            if(residual>maximum_residual)
                maximum_residual=residual;
            
            T lambda_normal_delta=-normal_violation/contact.normal_diagonal;
            T lambda_normal_new=lambda_normal(i)+lambda_normal_delta;
            if(lambda_normal_new<0)
            {
                lambda_normal_new=0;
                lambda_normal_delta=-lambda_normal(i);
            }
            
            if(!has_infinite_inertia(contact.id(0)))
                velocities(contact.id(0))+=lambda_normal_delta*contact.inverse_mass_times_normal_constraint(0);
            if(!has_infinite_inertia(contact.id(1)))
                velocities(contact.id(1))+=lambda_normal_delta*contact.inverse_mass_times_normal_constraint(1);

            ////////////////////////
            /*T new_normal_violation=T_TWIST::Dot_Product(contact.normal_constraint(0),velocities(contact.id(0)))+T_TWIST::Dot_Product(contact.normal_constraint(1),velocities(contact.id(1)))-contact.normal_relative_velocity;
            if(-new_normal_violation>tolerance)
            {
                LOG::cout << "constraint " << normal_violation << " " << new_normal_violation << std::endl;
                exit(0);
            }*/
            ///////////////////////
            
            lambda_normal(i)=lambda_normal_new;

            for(int j=1;j<d;j++)
            {
                T tangent_violation=T_TWIST::Dot_Product(contact.tangent_constraint(0)(j),velocities(contact.id(0)))+T_TWIST::Dot_Product(contact.tangent_constraint(1)(j),velocities(contact.id(1)));
                
                T lambda_tangent_new=lambda_tangent(i)(j)-tangent_violation/contact.tangent_diagonal(j);
                if(fabs(lambda_tangent_new)>contact.coefficient_of_friction*lambda_normal_new)
                {
                    if(lambda_tangent_new>0)
                        lambda_tangent_new=contact.coefficient_of_friction*lambda_normal_new;
                    else
                        lambda_tangent_new=-contact.coefficient_of_friction*lambda_normal_new;
                }

                if(!has_infinite_inertia(contact.id(0)))
                    velocities(contact.id(0))+=(lambda_tangent_new-lambda_tangent(i)(j))*contact.inverse_mass_times_tangent_constraint(0)(j);
                if(!has_infinite_inertia(contact.id(1)))
                    velocities(contact.id(1))+=(lambda_tangent_new-lambda_tangent(i)(j))*contact.inverse_mass_times_tangent_constraint(1)(j);

                lambda_tangent(i)(j)=lambda_tangent_new;
            }
        }

        //LOG::cout << "maximum_residual " << iteration << " " << maximum_residual << std::endl;
        
        iteration++;
    }

    //LOG::cout << "iterations = " << iteration << std::endl;
    //LOG::cout << "maximum_residual = " << maximum_residual << std::endl;

    //LOG::cout << "lambda_normal = " << lambda_normal << std::endl;
    //LOG::cout << "lambda_tangent = " << lambda_tangent << std::endl;

    /*for(int i=0;i<n_contacts;i++)
    {
        SOLVE_CONTACT::CONTACT<TV>& contact=contacts(i);
        
        int id1=contact.id(0),id2=contact.id(1);
        VECTOR<int,2> pair=VECTOR<int,2>(id1,id2).Sorted();
        if(pair(0)==1 && pair(1)==6)
        {
            LOG::cout << "contact " << i << " - " << contact.normal << " - " << contact.location << " - " << contact.distance << " " << contact.normal_relative_velocity << std::endl;
            LOG::cout << "residual " << TWIST<TV>::Dot_Product(contact.normal_constraint(0),velocities(id1)) << " " << (TWIST<TV>::Dot_Product(contact.normal_constraint(0),velocities(id1))+TWIST<TV>::Dot_Product(contact.normal_constraint(1),velocities(id2))-contact.normal_relative_velocity) << std::endl;
        }
    }*/

    LOG::cout << "pgs iterations " << iteration << " residual " << maximum_residual << std::endl;

    for(int i=0;i<n_bodies;i++)
    {
        if(rigid_body_collection.Is_Active(i) && !has_infinite_inertia(i))
        {
            RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
            //LOG::cout << "old velocity " << i << body.Twist() << std::endl;
            //LOG::cout << "new velocity " << i << velocities(i) << std::endl;
            body.Twist()=velocities(i);
            //body.Update_Angular_Momentum();
        }
    }

    //exit(0);

    return true;
}

}
}

#endif
