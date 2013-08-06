//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Symbolics/PROGRAM.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <climits>
#include <boost/function.hpp>
using namespace PhysBAM;

template<class T>
class MINIMIZATION_OBJECTIVE:public NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>
{
public:
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > VECTOR_T;
    typedef MATRIX_SYSTEM<MATRIX_MXN<T>,T,VECTOR_T> SYSTEM_T;

    ARRAY<int> grad_vars;
    ARRAY<int> hess_vars;
    PROGRAM<T> program;

    MINIMIZATION_OBJECTIVE(std::string function,std::string vars,char out)
    {
        grad_vars.Resize(vars.size());
        hess_vars.Resize(vars.size()*vars.size());
        char buff[2]={out,0};
        program.var_out.Append(buff);
        for(int i=0;i<grad_vars.m;i++){
            buff[0]=vars[i];
            program.var_in.Append(buff);}
        program.Parse(function.c_str(),false);
        for(int i=0;i<grad_vars.m;i++) grad_vars(i)=program.Diff(0,i);
        for(int i=0;i<grad_vars.m;i++) for(int j=0;j<grad_vars.m;j++) hess_vars(i*grad_vars.m+j)=program.Diff(grad_vars(i),j);
        program.Optimize();
        program.Finalize();
        program.Print();
    }

    virtual ~MINIMIZATION_OBJECTIVE(){}

    void Compute(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const PHYSBAM_OVERRIDE
    {
        PROGRAM_CONTEXT<T> context(program);
        context.data_in=static_cast<const VECTOR_T&>(x).v;
        program.Execute(context);
        if(e) *e=context.data_out(0);
        if(g) for(int i=0;i<grad_vars.m;i++) static_cast<VECTOR_T*>(g)->v(i)=context.data_out(grad_vars(i));
        if(h) for(int i=0;i<grad_vars.m;i++) for(int j=0;j<grad_vars.m;j++) const_cast<MATRIX_MXN<T>&>(static_cast<SYSTEM_T*>(h)->A)(i,j)=context.data_out(hess_vars(i*grad_vars.m+j));
    }
};

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    typedef double T;
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T>> VECTOR_T;
    typedef MATRIX_SYSTEM<MATRIX_MXN<T>,T,VECTOR_T> SYSTEM_T;

    const char* vars=argv[1];
    char out_var=*argv[2];
    const char* function=argv[3];

    int n=strlen(vars);
    MATRIX_MXN<T> matrix(n,n);
    SYSTEM_T system(matrix);

    MINIMIZATION_OBJECTIVE<T> obj(function, vars, out_var);

    NEWTONS_METHOD<T> nm;
    nm.max_iterations=100000;
    nm.max_krylov_iterations=2000;
    nm.krylov_tolerance=1;
    nm.fail_on_krylov_not_converged=false;
    nm.tolerance=1e-5;
    nm.angle_tolerance=1e-2;
//    nm.use_golden_section_search=true;
//    nm.use_wolfe_search=false;

    VECTOR_T x0;
    x0.v.Resize(n);
    for(int i=0;i<n;i++) x0.v(i)=atof(argv[4+i]);
    obj.Test(x0,system);
    LOG::printf("Initial guess: %.16P\n", x0.v);

    bool converged=nm.Newtons_Method(obj,system,x0);
    LOG::printf("Final solution: %.16P\n", x0.v);
    PHYSBAM_ASSERT(converged);


    return 0;
}
