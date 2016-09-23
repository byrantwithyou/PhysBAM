//#####################################################################
// Copyright 2011, Mauricio Flores, Leah Fritter, Gina Ma.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MINRES
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <cfloat>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class T> MINRES<T>::
~MINRES()
{}
//#####################################################################
// Function Print_Diagnostic
//#####################################################################
template<class T> void MINRES<T>::
Print_Diagnostics(int iterations)
{
    if(print_diagnostics) LOG::Stat("minres iterations",iterations);
    if(iterations_used) *iterations_used=iterations;
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> bool MINRES<T>::
Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
    ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,T tolerance,const int min_iterations,const int max_iterations)
{
    Ensure_Size(av,x,9);
    KRYLOV_VECTOR_BASE<T>& u=*av(0);
    KRYLOV_VECTOR_BASE<T>& w=*av(1);
    KRYLOV_VECTOR_BASE<T>& w_hat=*av(2);
    KRYLOV_VECTOR_BASE<T>& c_k=*av(3);
    KRYLOV_VECTOR_BASE<T>& c_p=*av(4);
    KRYLOV_VECTOR_BASE<T>& c_pp=*av(5);
    KRYLOV_VECTOR_BASE<T>& vk=*av(6);
    KRYLOV_VECTOR_BASE<T>& vk_hat=*av(7);
    KRYLOV_VECTOR_BASE<T>& vp=*av(8);
    KRYLOV_VECTOR_BASE<T>& vp_hat=*av(3); // aliases c_k
    KRYLOV_VECTOR_BASE<T>& z=*av(2); // aliases w_hat
    KRYLOV_VECTOR_BASE<T>& vtemp=*av(8); // aliases vp

    static const T small_number=std::numeric_limits<T>::epsilon(); //Set small number
    system.Set_Boundary_Conditions(x); //Just leave there
    T residual;

    //variable setup
    VECTOR<T,2> G_pp,G_p;
    T b_k=0;
    VECTOR<T,2> r;

    for(int iterations=0;;iterations++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0); //This is definitely not clear

        if(restart){
            if(print_residuals) LOG::cout<<"restarting Minres"<<std::endl;
            //Include initial settings here, avoid definitions.
            G_pp=VECTOR<T,2>(1,0);
            G_p=VECTOR<T,2>(1,0);
            system.Multiply(x,vtemp);
            vtemp.Copy(-1,vtemp,b);
            system.Project(vtemp);
            const KRYLOV_VECTOR_BASE<T>& b_hat=system.Precondition(vtemp,z); //Hoping to make b_hat=M*b;
            T theta=sqrt((T) system.Inner_Product(vtemp,b_hat));
            if(relative_tolerance && iterations==0) tolerance*=theta;
            if(theta<=0){Print_Diagnostics(iterations);return true;}
            vk_hat.Copy(1/theta, b_hat);
            vk.Copy(1/theta, vtemp);
            r=VECTOR<T,2>(theta,0);
            b_k=0;

            //Make zero vectors:
            c_p*=0;
            c_pp*=0;
            vp*=0;
            vp_hat*=0;}

        //Lanczos process
        system.Multiply(vk_hat, u); //u=A*vk_hat;
        system.Project(u);
        const KRYLOV_VECTOR_BASE<T>& u_hat=system.Precondition(u, z); //u_hat=M*u;
        T a_k=system.Inner_Product(vk, u_hat);
        w.Copy(-a_k, vk, u);
        w.Copy(-b_k, vp, w);
        w_hat.Copy(-a_k,vk_hat,u_hat);
        w_hat.Copy(-b_k,vp_hat,w_hat);
        T b_k1=sqrt(system.Inner_Product(w,w_hat));

        //Obtain epsilon, delta, givens rotation matrix, and gamma
        VECTOR<T,2> bl(0,b_k);
        VECTOR<T,2> et=G_pp.Givens_Transpose_Times(bl);
        VECTOR<T,2> ta(et(1),a_k);
        VECTOR<T,2> ds=G_p.Givens_Transpose_Times(ta);
        VECTOR<T,2> G_k(ds(1),b_k1);
        T gamma=G_k.Normalize();
        if(gamma<small_number){Print_Diagnostics(iterations);LOG::cout << "gamma variable close to zero"<<std::endl;return false;}

        //Obtain p, c, and r (residual)
        VECTOR<T,2> pr=G_k.Givens_Transpose_Times(r);
        residual=abs(pr(1));

        vtemp.Copy(1/gamma,vk_hat);
        c_k.Copy(-et(0)/gamma,c_pp,vtemp);
        c_k.Copy(-ds(0)/gamma,c_p,c_k);

        x.Copy(pr(0), c_k, x);

        if(print_residuals) LOG::cout<<residual<<std::endl;
        if(residual<=tolerance || b_k1<small_number){Print_Diagnostics(iterations);return true;}
        if(iterations==max_iterations){Print_Diagnostics(iterations);break;}

        //Reassignment:
        G_pp=G_p;
        G_p=G_k;
        c_pp=c_p;
        c_p=c_k;
        vp=vk;
        vp_hat=vk_hat;
        vk.Copy(1/b_k1, w); //added 7/21
        vk_hat.Copy(1/b_k1, w_hat); //modified 7/21
        b_k=b_k1;
        r=VECTOR<T,2>(pr(1),0);}

    if(print_diagnostics) LOG::cout<<"Minres not converged after "<<max_iterations<<" iterations, error="<<residual<<std::endl;
    return false;
}
//#####################################################################
namespace PhysBAM{
template class MINRES<float>;
template class MINRES<double>;
}
