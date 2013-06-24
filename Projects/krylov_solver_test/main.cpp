//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

template<class T>
class KRYLOV_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
{
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;

public:

    MATRIX_MXN<T>& matrix;

    KRYLOV_SYSTEM(MATRIX_MXN<T>& matrix_input):BASE(false,false),matrix(matrix_input){}
    virtual ~KRYLOV_SYSTEM(){}

//#####################################################################
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const {debug_cast<VECTOR_T&>(result).v=matrix*debug_cast<const VECTOR_T&>(x).v;}
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const {return debug_cast<const VECTOR_T&>(x).v.Dot(debug_cast<const VECTOR_T&>(y).v);}
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const {return debug_cast<const VECTOR_T&>(x).v.Max_Abs();}
    void Project(KRYLOV_VECTOR_BASE<T>& x) const {}
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const {}
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const {}
//#####################################################################
};

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    typedef double T;

    MATRIX_MXN<T> M(3,3);
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > b,x;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    b.v.Resize(3);
    x.v.Resize(3);

    // M(0,0)=1;M(0,1)=3;M(0,2)=5;
    // M(1,0)=3;M(1,1)=5;M(1,2)=7;
    // M(2,0)=5;M(2,1)=7;M(2,2)=10;

    M(0,0)=14;M(0,1)=32;M(0,2)=53;
    M(1,0)=32;M(1,1)=77;M(1,2)=128;
    M(2,0)=53;M(2,1)=128;M(2,2)=213;

    b.v(0)=1;b.v(1)=2;b.v(2)=4;

    KRYLOV_SYSTEM<T> system(M);

    LOG::cout<<"CONJUGATE_RESIDUAL"<<std::endl;
    CONJUGATE_RESIDUAL<T> cr;
    cr.print_residuals=true;
    x.v.Fill(T());
    cr.Solve(system,x,b,vectors,1e-10,0,5000000);
    LOG::cout<<x.v<<std::endl;

    LOG::cout<<"MINRES"<<std::endl;
    MINRES<T> mr;
    mr.print_residuals=true;
    x.v.Fill(T());
    mr.Solve(system,x,b,vectors,1e-10,0,5000000);
    LOG::cout<<x.v<<std::endl;

    return 0;
}
