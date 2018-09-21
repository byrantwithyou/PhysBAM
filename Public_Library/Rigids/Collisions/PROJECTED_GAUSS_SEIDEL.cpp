//#####################################################################
// Copyright 2010, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Collisions/PROJECTED_GAUSS_SEIDEL.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>

namespace PhysBAM{

namespace PROJECTED_GAUSS_SEIDEL
{
template<class T>
void Solve(SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& a,ARRAY<T>& x,T tolerance)
{
    assert(A.n==a.Size() && A.n==x.Size());

    int n=A.n;
    x.Fill(0);

    int iteration=0;

    T maximum_residual=1;
    while(maximum_residual>tolerance)
    {
        maximum_residual=0;
        T row_sum=0;
        T diagonal=1;
        A.For_Each(
            [&](int i){row_sum=0;diagonal=1;},
            [&](int i,int j,T a){if(j==i) diagonal=a;else row_sum+=a*x(j);},
            [&](int i)
            {
                T residual=row_sum+diagonal*x(i)-a(i);
                if(-residual>maximum_residual) maximum_residual=-residual;
                x(i)=(a(i)-row_sum)/diagonal;
                if(x(i)<0) x(i)=0;
            });
        iteration++;
    }

    LOG::cout << "PGS iterations = " << iteration << std::endl;
}

template<class T,int D>
void Multiply(SPARSE_MATRIX_FLAT_MXN<VECTOR<T,D> >& A,ARRAY<VECTOR<T,D> >& x,ARRAY<T>& result)
{
    assert(A.n==x.Size()&&A.m==result.Size());

    T r=0;
    A.For_Each(
        [&](int i){r=0;},
        [&](int i,int j,const VECTOR<T,D>& a){r+=a.Dot(x(j));},
        [&](int i){result(i)=r;});
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

        iteration++;
    }

    LOG::cout << "PGS iterations = " << iteration << std::endl;
    LOG::cout << "PGS maximum_residual = " << maximum_residual << std::endl;
}

//#####################################################################
// Function Solve_Projected_Gauss_Seidel
//#####################################################################
template<class TV>
bool Solve(ARRAY<TWIST<TV> >& velocities,
    ARRAY<bool>& has_infinite_inertia,
    ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts,
    ARRAY<typename TV::SCALAR>& lambda_normal,
    ARRAY<VECTOR<typename TV::SCALAR,TV::m-1> >& lambda_tangent,
    typename TV::SCALAR tolerance,
    int iteration_maximum,
    bool friction)
{
    LOG::SCOPE scope("PROJECTED_GAUSS_SEIDEL::Solve");

    typedef typename TV::SCALAR T;
    typedef TWIST<TV> T_TWIST;
    const int d=TV::m;

    int n_contacts=contacts.m;

    T maximum_residual=2*tolerance+1;
    
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
    typedef typename TV::SCALAR T;
    typedef TWIST<TV> T_TWIST;
    const int d=TV::m;

    int n_bodies=rigid_body_collection.rigid_body_particles.Size();
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

        iteration++;
    }

    LOG::cout << "pgs iterations " << iteration << " residual " << maximum_residual << std::endl;

    for(int i=0;i<n_bodies;i++)
    {
        if(rigid_body_collection.Is_Active(i) && !has_infinite_inertia(i))
        {
            RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
            body.Twist()=velocities(i);
        }
    }

    return true;
}
template bool Solve<VECTOR<double,1> >(RIGID_BODY_COLLECTION<VECTOR<double,1> >&,
    ARRAY<SOLVE_CONTACT::CONTACT<VECTOR<double,1> >,int>&,VECTOR<double,1>::SCALAR,int);
template bool Solve<VECTOR<double,2> >(RIGID_BODY_COLLECTION<VECTOR<double,2> >&,
    ARRAY<SOLVE_CONTACT::CONTACT<VECTOR<double,2> >,int>&,VECTOR<double,2>::SCALAR,int);
template bool Solve<VECTOR<double,3> >(RIGID_BODY_COLLECTION<VECTOR<double,3> >&,
    ARRAY<SOLVE_CONTACT::CONTACT<VECTOR<double,3> >,int>&,VECTOR<double,3>::SCALAR,int);
template bool Solve<VECTOR<float,1> >(RIGID_BODY_COLLECTION<VECTOR<float,1> >&,
    ARRAY<SOLVE_CONTACT::CONTACT<VECTOR<float,1> >,int>&,VECTOR<float,1>::SCALAR,int);
template bool Solve<VECTOR<float,2> >(RIGID_BODY_COLLECTION<VECTOR<float,2> >&,
    ARRAY<SOLVE_CONTACT::CONTACT<VECTOR<float,2> >,int>&,VECTOR<float,2>::SCALAR,int);
template bool Solve<VECTOR<float,3> >(RIGID_BODY_COLLECTION<VECTOR<float,3> >&,
    ARRAY<SOLVE_CONTACT::CONTACT<VECTOR<float,3> >,int>&,VECTOR<float,3>::SCALAR,int);
}
}
