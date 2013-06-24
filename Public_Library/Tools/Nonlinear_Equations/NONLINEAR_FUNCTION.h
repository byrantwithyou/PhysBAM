//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONLINEAR_FUNCTION
//#####################################################################
#ifndef __NONLINEAR_FUNCTION__
#define __NONLINEAR_FUNCTION__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
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
    virtual R operator()(const T1 x) const=0;
    virtual R Prime(const T1 x) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual R Prime_Prime(const T1 x) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
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

template<class T,class R,class T1,class T2>
class PARAMETRIC_LINE<T,R(T1,T2)>:public NONLINEAR_FUNCTION<R(T)>
{
public:
    const NONLINEAR_FUNCTION<R(T1,T2)>& f;
    T1 x_not,direction_x;
    T2 y_not,direction_y;

    PARAMETRIC_LINE(const NONLINEAR_FUNCTION<R(T1,T2)>& f,const T1 x,const T2 y,const T1 a,const T2 b)
        :f(f),x_not(x),direction_x(a),y_not(y),direction_y(b)
    {}

    R operator()(const T t) const PHYSBAM_OVERRIDE
    {return f(x_not+t*direction_x,y_not+t*direction_y);}

//#####################################################################
};

template<class T>
struct PARAMETER_SPACE
{
    virtual ~PARAMETER_SPACE(){}

    virtual PARAMETER_SPACE<T>& Zero_Clone() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Op(T a,const PARAMETER_SPACE& x,T b,const PARAMETER_SPACE& y)=0; // this=a*x+b*y
    virtual void Copy(const PARAMETER_SPACE& x) {Op(1,x,0,*this);} // this=x
    virtual void Zero() {Op(0,*this,0,*this);} // this=0
    virtual T Dot(const PARAMETER_SPACE& x) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();} // this dot x
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
class PARAMETRIC_LINE<T,T(PARAMETER_SPACE<T>&)>:public NONLINEAR_FUNCTION<T(T)>
{
public:
    typedef PARAMETER_SPACE<T> T_PARAMETER_SPACE;
    typedef NONLINEAR_FUNCTION<T(T_PARAMETER_SPACE&)> F;
    const F& f;
    const T_PARAMETER_SPACE &x,&dx;
    T_PARAMETER_SPACE& tmp;

    PARAMETRIC_LINE(const F& f,const T_PARAMETER_SPACE& x,const T_PARAMETER_SPACE& dx,T_PARAMETER_SPACE& tmp)
        :f(f),x(x),dx(dx),tmp(tmp)
    {}

    T operator()(const T t) const PHYSBAM_OVERRIDE
    {tmp.Op(1,x,t,dx);return f(tmp);}
//#####################################################################
};

template<class T>
class NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>
{
public:
    virtual ~NONLINEAR_FUNCTION(){}
//#####################################################################
    virtual void Compute(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const=0;
    virtual T operator()(const KRYLOV_VECTOR_BASE<T>& x) const {T E=0;Compute(x,0,0,&E);return E;}
    void Test(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>& h) const;
//#####################################################################
};

template<class T>
class PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)>:public NONLINEAR_FUNCTION<T(T)>
{
public:
    typedef NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)> F;
    const F& f;
    const KRYLOV_VECTOR_BASE<T> &x,&dx;
    KRYLOV_VECTOR_BASE<T>& tmp;

    PARAMETRIC_LINE(const F& f,const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& dx,KRYLOV_VECTOR_BASE<T>& tmp)
        :f(f),x(x),dx(dx),tmp(tmp)
    {}

    T operator()(const T t) const PHYSBAM_OVERRIDE
    {tmp.Copy(t,dx,x);return f(tmp);}
//#####################################################################
};

}
#endif
