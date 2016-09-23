//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONLINEAR_FUNCTION
//#####################################################################
#ifndef __NONLINEAR_FUNCTION__
#define __NONLINEAR_FUNCTION__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Tools/Nonlinear_Equations/PARAMETER_SPACE.h>
namespace PhysBAM{

template<class F> class NONLINEAR_FUNCTION; // F must be a function type (e.g., T(T) or T(T1,T2))
template<class T,class F> class PARAMETRIC_LINE; // F must be a two argument function type
template<class T> class KRYLOV_SYSTEM_BASE;
template<class T> class KRYLOV_VECTOR_BASE;

template<class R,class T1>
class NONLINEAR_FUNCTION<R(T1)>
{
public:
    virtual ~NONLINEAR_FUNCTION(){}
//#####################################################################
    R operator()(const T1 x) const {R r=0;Compute(x,0,0,&r);return r;}
    R Prime(const T1 x) const {R r=0;Compute(x,0,&r,0);return r;}
    R Prime_Prime(const T1 x) const {R r=0;Compute(x,&r,0,0);return r;}
    virtual void Compute(const T1 x,R* ddf,R* df,R* f) const=0;
    void Test(const T1 x,bool test_second_diff);
//#####################################################################
};

template<class R,class T1,class T2>
class NONLINEAR_FUNCTION<R(T1,T2)>
{
public:
    virtual ~NONLINEAR_FUNCTION() {}
//#####################################################################
    virtual R operator()(const T1 x,const T2 y) const=0;
    virtual R Partial_X(const T1 x,const T2 y) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual R Partial_Y(const T1 x,const T2 y) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};

template<class T>
class NONLINEAR_FUNCTION<T(PARAMETER_SPACE<T>&)>
{
public:
    typedef PARAMETER_SPACE<T> T_PARAMETER_SPACE;

    virtual ~NONLINEAR_FUNCTION(){}
//#####################################################################
    virtual T operator()(const T_PARAMETER_SPACE& x) const=0;
    virtual void Gradient(const T_PARAMETER_SPACE& x,T_PARAMETER_SPACE& g) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();} // g = grad(x)
    virtual void Times_Hessian(const T_PARAMETER_SPACE& x,const T_PARAMETER_SPACE& y,T_PARAMETER_SPACE& z) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();} // z = H(x) y
//#####################################################################
};

template<class T>
class NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>
{
public:
    virtual ~NONLINEAR_FUNCTION();
//#####################################################################
    virtual void Compute(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const=0;
    virtual T operator()(const KRYLOV_VECTOR_BASE<T>& x) const;
    void Test(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>& h) const;
    virtual void Make_Feasible(KRYLOV_VECTOR_BASE<T>& x) const;
//#####################################################################
};
}
#endif
