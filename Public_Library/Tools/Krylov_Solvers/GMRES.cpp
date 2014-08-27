//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GMRES
//#####################################################################
#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Vectors/VECTOR.h>
#include <cfloat>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class T> GMRES<T>::
~GMRES()
{}
//#####################################################################
// Function Print_Diagnostic
//#####################################################################
template<class T> void GMRES<T>::
Print_Diagnostics(int iterations)
{
    if(print_diagnostics) LOG::Stat("gmres iterations",iterations);
    if(iterations_used) *iterations_used=iterations;
}
//#####################################################################
// Function Get_Vector
//#####################################################################
template<class T> static KRYLOV_VECTOR_BASE<T>*
Get_Vector(ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,int& next)
{
    assert(next<=av.m);
    if(next==av.m) av.Append(av(0)->Clone_Default());
    return av(next++);
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> bool GMRES<T>::
Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
    ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,T tolerance,const int min_iterations,const int max_iterations)
{
    int marker=1+system.use_preconditioner;
    int next_vector=marker;
    Ensure_Size(av,x,std::max(marker,av.m));
    KRYLOV_VECTOR_BASE<T>& temp_1=*av(0);
    KRYLOV_VECTOR_BASE<T>& temp_2=*av(system.use_preconditioner);

    T local_tolerance=tolerance;
    if(relative_tolerance){
        const KRYLOV_VECTOR_BASE<T>& prec_b=system.Precondition(b,temp_2);
        local_tolerance=tolerance*sqrt(system.Inner_Product(prec_b,prec_b));}

    T residual,phi;
    ARRAY<T> t;
    ARRAY<ARRAY<T> > h,r_2;
    ARRAY<VECTOR<T,2> > rot;
    KRYLOV_VECTOR_BASE<T> *v_cur=0;
    for(int iterations=0,system_size;;iterations++,system_size++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0);
        if(restart){
            system_size=1;
            next_vector=marker;
            v_cur=Get_Vector(av,next_vector);
            h.Remove_All();
            r_2.Remove_All();
            rot.Remove_All();
            t.Remove_All();
            system.Multiply(x,temp_1);
            temp_1.Copy(-1,temp_1,b);
            const KRYLOV_VECTOR_BASE<T>& prec=system.Precondition(temp_1,temp_2);
            phi=sqrt(system.Inner_Product(prec,prec));
            v_cur->Copy(1/phi,prec);}

        system.Multiply(*v_cur,temp_1);
        const KRYLOV_VECTOR_BASE<T>& prec2=system.Precondition(temp_1,temp_2);
        if(&prec2!=&temp_1) temp_1=prec2;

        ARRAY<T>& r_2_cur=r_2(r_2.Append(ARRAY<T>()));
        ARRAY<T>& h_cur=h(h.Append(ARRAY<T>()));
        for(int i=0;i<system_size;i++){
            T ip=system.Inner_Product(temp_1,*av(marker+i));
            h_cur.Append(ip);
            temp_1.Copy(-ip,*av(marker+i),temp_1);}

        KRYLOV_VECTOR_BASE<T> *v_next=Get_Vector(av,next_vector);
        T ip=sqrt(system.Inner_Product(temp_1,temp_1));
        h_cur.Append(ip);
        v_next->Copy(1/ip,temp_1);

        T r=h_cur(0);
        for(int i=0;i<system_size-1;i++){
            VECTOR<T,2> tmp=rot(i).Givens_Transpose_Times(VECTOR<T,2>(r,h_cur(i+1)));
            r_2_cur.Append(tmp.x);
            r=tmp.y;}
        VECTOR<T,2> tmp(r,h_cur(system_size));
        r_2_cur.Append(tmp.Normalize());
        rot.Append(tmp);

        VECTOR<T,2> ph(rot(system_size-1)*phi);
        t.Append(ph.x);
        phi=-ph.y;

        residual=abs(phi);
        residual_magnitude_squared=residual*residual;
        if(residual<=local_tolerance || iterations==max_iterations || (residual>tolerance && restart_iterations && system_size==restart_iterations)){
            for(int i=r_2.m-1;i>=0;i--){
                T rhs=t(i);
                for(int j=i+1;j<r_2.m;j++)
                    rhs-=r_2(j)(i)*t(j);
                t(i)=rhs/(r_2(i)(i));
                x.Copy(t(i),*av(marker+i),x);}
            if(residual<=local_tolerance){Print_Diagnostics(iterations);return true;}
            if(iterations==max_iterations){Print_Diagnostics(iterations);break;}}
        v_cur=v_next;}
    if(print_diagnostics){LOG::cout<<"gmres has not converged after "<<max_iterations<<" iterations, norm of residual="<<residual<<std::endl;}
    return false;
}
//#####################################################################
namespace PhysBAM{
template class GMRES<float>;
template class GMRES<double>;
}
