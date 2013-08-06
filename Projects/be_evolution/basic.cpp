//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <climits>
#include <boost/function.hpp>
using namespace PhysBAM;

template<class TV>
class MINIMIZATION_OBJECTIVE:public NONLINEAR_FUNCTION<typename TV::SCALAR(KRYLOV_VECTOR_BASE<typename TV::SCALAR>&)>
{
    typedef typename TV::SCALAR T;
public:
    typedef KRYLOV_VECTOR_WRAPPER<T,TV> VECTOR_T;
    typedef MATRIX_SYSTEM<MATRIX<T,TV::m>,T,VECTOR_T> SYSTEM_T;

    boost::function<T(TV)> f;
    boost::function<TV(TV)> df;
    boost::function<MATRIX<T,TV::m>(TV)> ddf;

    MINIMIZATION_OBJECTIVE(boost::function<T(TV)> f,boost::function<TV(TV)> df,boost::function<MATRIX<T,TV::m>(TV)> ddf)
        :f(f),df(df),ddf(ddf)
    {}

    virtual ~MINIMIZATION_OBJECTIVE(){}

    void Compute(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const PHYSBAM_OVERRIDE
    {
        TV y=static_cast<const VECTOR_T&>(x).v;
        if(e) *e=f(y);
        if(g) static_cast<VECTOR_T*>(g)->v=df(y);
        if(h) const_cast<MATRIX<T,TV::m>&>(static_cast<SYSTEM_T*>(h)->A)=ddf(y);
    }
};

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    typedef double T;
    typedef VECTOR<T,2> TV;
    typedef KRYLOV_VECTOR_WRAPPER<T,TV> VECTOR_T;
    typedef MATRIX_SYSTEM<MATRIX<T,TV::m>,T,VECTOR_T> SYSTEM_T;

    MATRIX<T,TV::m> matrix;
    SYSTEM_T system(matrix);

    MINIMIZATION_OBJECTIVE<TV> obj([](TV x){return sqr(x.Magnitude_Squared()-1);},
        [](TV x){return 4*(x.Magnitude_Squared()-1)*x;},
        [](TV x){return MATRIX<T,TV::m>::Outer_Product(x,x)*8+4*(x.Magnitude_Squared()-1);});

    NEWTONS_METHOD<T> nm;
    nm.max_iterations=100000;
    nm.max_krylov_iterations=2000;
    nm.krylov_tolerance=1;
    nm.fail_on_krylov_not_converged=false;
    nm.tolerance=1e-5;
    nm.angle_tolerance=1e-2;

    VECTOR_T x0;
    x0.v=TV(1e-6,1+1e-5);
    obj.Test(x0,system);

    bool converged=nm.Newtons_Method(obj,system,x0);
    PHYSBAM_ASSERT(converged);


    return 0;
}
