//#####################################################################
// Copyright 2011, Mauricio Flores, Leah Fritter, Gina Ma.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LANCZOS_CONJUGATE_GRADIENT
//#####################################################################
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Krylov_Solvers/LANCZOS_CONJUGATE_GRADIENT.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Vectors/VECTOR.h>
#include <cfloat>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class T> LANCZOS_CONJUGATE_GRADIENT<T>::
~LANCZOS_CONJUGATE_GRADIENT()
{}
//#####################################################################
// Function Print_Diagnostic
//#####################################################################
template<class T> void LANCZOS_CONJUGATE_GRADIENT<T>::
Print_Diagnostics(int iterations)
{
    if(print_diagnostics) LOG::Stat("cg iterations",iterations);
    if(iterations_used) *iterations_used=iterations;
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> bool LANCZOS_CONJUGATE_GRADIENT<T>::
Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
    ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,T tolerance,const int min_iterations,const int max_iterations)
{
    Ensure_Size(av,x,5);
    KRYLOV_VECTOR_BASE<T>& vk=*av(0);
    KRYLOV_VECTOR_BASE<T>& vk1=*av(1);
    KRYLOV_VECTOR_BASE<T>& vkm1=*av(4);
    KRYLOV_VECTOR_BASE<T>& ck=*av(2);
    KRYLOV_VECTOR_BASE<T>& ckm1=*av(3);

    static const T small_number=std::numeric_limits<T>::epsilon(); //Set small number
    system.Set_Boundary_Conditions(x); //Just leave there

    T bk=0, dkm1=0, lk=0, qk=bk, residual;

    for(int iterations=0;;iterations++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0); //This is definitely not clear

        if(restart){
            if(print_residuals) LOG::cout<<"restarting CG"<<std::endl;
            system.Multiply(x,vk);
            vk.Copy(-1,vk,b);
            T th=sqrt(system.Inner_Product(vk,vk));
            if(th<=0){Print_Diagnostics(iterations);return true;}
            if(relative_tolerance && iterations==0) tolerance*=th;
            vk.Copy(1/th,vk);
            vkm1*=0;
            dkm1=0;
            lk=0;
            qk=th;
            bk=0;
            ckm1*=0;}

        //Lanczos process
        system.Multiply(vk, vk1);
        T ak = system.Inner_Product(vk, vk1);
        vk1.Copy(-ak, vk, vk1);
        vk1.Copy(-bk, vkm1, vk1);
        T bk1 = sqrt(system.Inner_Product(vk1,vk1));

        T dk=ak-dkm1*sqr(lk);
        T lk1=bk1/dk;
        T pk=qk/dk;

        ck.Copy(-lk,ckm1,vk);
        x.Copy(pk, ck, x);
        residual=bk1*std::abs(pk);
        
        if(print_residuals) LOG::cout<<residual<<std::endl;
        if(bk1<small_number){Print_Diagnostics(iterations+1);return true;}
        if(residual <= tolerance || iterations == max_iterations-1){Print_Diagnostics(iterations+1);break;}

        vkm1=vk;
        vk.Copy(1/bk1, vk1);
        bk=bk1;
        dkm1=dk;
        lk=lk1;
        qk=-lk*qk;
        ckm1=ck;
    }

    if(print_diagnostics) LOG::cout<<"cg not converged after "<<max_iterations<<" iterations, error = "<<residual<<std::endl;
    return false;
}
//#####################################################################
namespace PhysBAM{
template class LANCZOS_CONJUGATE_GRADIENT<float>;
template class LANCZOS_CONJUGATE_GRADIENT<double>;
}
