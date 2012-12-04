//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/MINRES.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
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

template<class T>
struct SQUARED_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
{
public:
    const KRYLOV_SYSTEM_BASE<T>& system;
    mutable KRYLOV_VECTOR_BASE<T>* tmp;
    mutable KRYLOV_VECTOR_BASE<T>* tmp2;

    SQUARED_SYSTEM(const KRYLOV_SYSTEM_BASE<T>& system_input,const bool use_preconditioner,const bool preconditioner_commutes_with_projection)
        :KRYLOV_SYSTEM_BASE<T>(use_preconditioner,preconditioner_commutes_with_projection),system(system_input),tmp(0),tmp2(0)
    {}

    virtual ~SQUARED_SYSTEM()
    {delete tmp;delete tmp2;}

    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const PHYSBAM_OVERRIDE
    {
        if(!tmp) tmp=x.Clone_Default();
        if(!tmp2) tmp2=x.Clone_Default();
        system.Multiply(x,*tmp);
        const KRYLOV_VECTOR_BASE<T>& mr=system.Precondition(*tmp,*tmp2);
        system.Multiply(mr,result);
    }

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const PHYSBAM_OVERRIDE
    {return system.Inner_Product(x,y);}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
    {return system.Convergence_Norm(x);}

    void Project(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
    {system.Project(x);}

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
    {system.Set_Boundary_Conditions(x);}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
    {system.Project_Nullspace(x);}

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const PHYSBAM_OVERRIDE
    {system.Precondition(r,z);}

//#####################################################################
};
//#####################################################################
// Function Nullspace_Check
//#####################################################################
template<class T> bool KRYLOV_SYSTEM_BASE<T>::
Nullspace_Check(KRYLOV_VECTOR_BASE<T>& null) const
{
    KRYLOV_VECTOR_BASE<T>* x=null.Clone_Default();
    KRYLOV_VECTOR_BASE<T>* b=null.Clone_Default();
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    RANDOM_NUMBERS<T> random;
    int n=null.Raw_Size();
    for(int i=0;i<n;i++)
        null.Raw_Get(i)=random.Get_Uniform_Number(-1,1);
    MINRES<T> mr;
//    cg.print_diagnostics=false;
    Project(null);
    T mg=sqrt(Inner_Product(null,null));
    null*=1/mg;

    T tol=1e-6;
    for(int i=0;tol>1e-16;i++){
        *x*=(T)0;
        *b*=(T)0;
        Multiply(null,*b);
        mr.Solve(*this,*x,*b,av,tol,0,100000);
        tol/=1e9;
        null-=*x;
        Project(null);
        mg=sqrt(Inner_Product(null,null));
        null*=1/mg;
        if(mg<1e-6) return false;
        LOG::cout<<(1-mg)<<std::endl;}
    delete x;
    delete b;
    return mg>1e-10;
}
namespace PhysBAM{
template class KRYLOV_SYSTEM_BASE<float>;
template class KRYLOV_SYSTEM_BASE<double>;
}
