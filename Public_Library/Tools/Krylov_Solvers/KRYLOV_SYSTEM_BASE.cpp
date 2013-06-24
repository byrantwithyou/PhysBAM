//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/maxabs.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Utilities/DEBUG_CAST.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> KRYLOV_SYSTEM_BASE<T>::
KRYLOV_SYSTEM_BASE(const bool use_preconditioner,const bool preconditioner_commutes_with_projection)
    :use_preconditioner(use_preconditioner),preconditioner_commutes_with_projection(preconditioner_commutes_with_projection)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> KRYLOV_SYSTEM_BASE<T>::
~KRYLOV_SYSTEM_BASE()
{
}
//#####################################################################
// Function Test_System
//#####################################################################
template<class T> void KRYLOV_SYSTEM_BASE<T>::
Test_System(KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& y,KRYLOV_VECTOR_BASE<T>& z) const
{
    T tolerance=(T)1e-5;
    RANDOM_NUMBERS<T> random;
    double a,b,r;
    bool pass=true;
    int n=x.Raw_Size();
    for(int i=0;i<n;i++){
        x.Raw_Get(i)=random.Get_Uniform_Number(-1,1);
        y.Raw_Get(i)=random.Get_Uniform_Number(-1,1);}

    a=Inner_Product(x,y);
    b=Inner_Product(y,x);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Inner Product Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<r<<std::endl;}

    Multiply(x,z);
    a=Inner_Product(y,z);
    Multiply(y,z);
    b=Inner_Product(x,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"System Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z=x;
    Project(z);
    a=Inner_Product(y,z);
    z=y;
    Project(z);
    b=Inner_Product(x,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Projection Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z=x;
    Project(z);
    a=Inner_Product(z,z);
    Project(z);
    b=Inner_Product(z,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Projection Idempotence Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z=x;
    Project_Nullspace(z);
    a=Inner_Product(y,z);
    z=y;
    Project_Nullspace(z);
    b=Inner_Product(x,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Project Nullspace Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z=x;
    Project_Nullspace(z);
    a=Inner_Product(z,z);
    Project_Nullspace(z);
    b=Inner_Product(z,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Project Nullspace Idempotence Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z*=(T)0;
    a=Inner_Product(y,Precondition(x,z));
    z*=(T)0;
    b=Inner_Product(x,Precondition(y,z));
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Preconditioner Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    LOG::cout<<"Krylov System Test Result: "<<(pass?"PASS":"FAIL")<<std::endl;
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class T> void KRYLOV_SYSTEM_BASE<T>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class T> void KRYLOV_SYSTEM_BASE<T>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Precondition
//#####################################################################
template<class T> const KRYLOV_VECTOR_BASE<T>& KRYLOV_SYSTEM_BASE<T>::
Precondition(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    if(!use_preconditioner) return r;
    Apply_Preconditioner(r,z);
    if(!preconditioner_commutes_with_projection){Project(z);Project_Nullspace(z);}
    return z;
}
//#####################################################################
// Function Compute_Nullspace
//#####################################################################
template<class T> void KRYLOV_SYSTEM_BASE<T>::
Compute_Nullspace(const KRYLOV_VECTOR_BASE<T>& tmp,ARRAY<KRYLOV_VECTOR_BASE<T>*>& null,int max_null) const
{
    KRYLOV_VECTOR_BASE<T>* x=tmp.Clone_Default();
    KRYLOV_VECTOR_BASE<T>* b=tmp.Clone_Default();
    KRYLOV_VECTOR_BASE<T>* z=tmp.Clone_Default();
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    RANDOM_NUMBERS<T> random;
    MINRES<T> mr;
    for(int j=0;j<max_null;j++){
        int n=z->Raw_Size();
        for(int i=0;i<n;i++)
            z->Raw_Get(i)=random.Get_Uniform_Number(-1,1);
        Project(*z);
        for(int i=0;i<null.m;i++)
            z->Copy(-Inner_Product(*z,*null(i)),*null(i),*z);
        T mg=sqrt(Inner_Product(*z,*z));
        *z*=1/mg;

        T tol=1e-15,tol_working=sqrt(tol);
        for(int i=0;i<3;i++){
            *x*=(T)0;
            *b*=(T)0;
            Multiply(*z,*b);
            mr.Solve(*this,*x,*b,av,tol_working,0,100000);
            *z-=*x;
            Project(*z);
            for(int i=0;i<null.m;i++)
                z->Copy(-Inner_Product(*z,*null(i)),*null(i),*z);
            mg=sqrt(Inner_Product(*z,*z));
            if(!mg) break;
            *z*=1/mg;
            LOG::cout<<"nullspace iteration "<<(1-mg)<<std::endl;
            tol_working=tol;}
        Multiply(*z,*b);
        mg=Inner_Product(*z,*b);
        if(abs(mg)>1e-10) break;
        LOG::cout<<"nullspace eigenvalue "<<mg<<std::endl;
        null.Append(z);
        z=tmp.Clone_Default();}
    delete x;
    delete b;
    delete z;
}
//#####################################################################
// Function Compute_Small_Eigenvectors
//#####################################################################
template<class T> void KRYLOV_SYSTEM_BASE<T>::
Compute_Small_Eigenvectors(const KRYLOV_VECTOR_BASE<T>& tmp,ARRAY<KRYLOV_VECTOR_BASE<T>*>& null,
    ARRAY<KRYLOV_VECTOR_BASE<T>*>& eigenvectors,ARRAY<T>& eigenvalues,int max_eigen,T tol,int power_iter) const
{
//    Compute_Nullspace(tmp,null,max_eigen);

    KRYLOV_VECTOR_BASE<T>* x=tmp.Clone_Default();
    KRYLOV_VECTOR_BASE<T>* z=tmp.Clone_Default();
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    RANDOM_NUMBERS<T> random;
    MINRES<T> mr;
    for(int j=null.m;j<max_eigen;j++){
        int n=z->Raw_Size();
        for(int i=0;i<n;i++)
            z->Raw_Get(i)=random.Get_Uniform_Number(-1,1);
        Project(*z);
        for(int i=0;i<null.m;i++)
            z->Copy(-Inner_Product(*z,*null(i)),*null(i),*z);
        for(int i=0;i<eigenvectors.m;i++)
            z->Copy(-Inner_Product(*z,*eigenvectors(i)),*eigenvectors(i),*z);
        T mg=sqrt(Inner_Product(*z,*z)),last=0;
        *z*=1/mg;

        for(int i=0;i<power_iter;i++){
            *x=*z;
            mr.Solve(*this,*x,*z,av,1e-15,0,1000);
            exchange(x,z);
            Project(*z);
            for(int i=0;i<null.m;i++)
                z->Copy(-Inner_Product(*z,*null(i)),*null(i),*z);
            for(int i=0;i<eigenvectors.m;i++)
                z->Copy(-Inner_Product(*z,*eigenvectors(i)),*eigenvectors(i),*z);
            mg=sqrt(Inner_Product(*z,*z));
            *z*=1/mg;
            LOG::cout<<"power method iter "<<1/mg<<std::endl;
            if(abs(mg-last)/maxabs(mg,last,(T)1e-12)<(T)1e-6) break;
            last=mg;}
        mg=1/abs(mg);
        if(mg>tol) break;
        eigenvectors.Append(z);
        eigenvalues.Append(mg);
        LOG::cout<<"eigenvalue "<<mg<<std::endl;
        z=tmp.Clone_Default();}
    delete x;
    delete z;
}
namespace PhysBAM{
template class KRYLOV_SYSTEM_BASE<float>;
template class KRYLOV_SYSTEM_BASE<double>;
}
