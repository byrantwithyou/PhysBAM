//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARAMETRIC_LINE
//#####################################################################
#ifndef __PARAMETRIC_LINE__
#define __PARAMETRIC_LINE__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T,class F> class PARAMETRIC_LINE; // F must be a two argument function type
template<class T> class KRYLOV_SYSTEM_BASE;
template<class T> class KRYLOV_VECTOR_BASE;

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

    virtual void Compute(const T1 t,R* ddg,R* dg,R* g) const PHYSBAM_OVERRIDE
    {assert(g && !dg && !ddg);*g=f(x_not+t*direction_x,y_not+t*direction_y);}

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

    virtual void Compute(const T t,T* ddg,T* dg,T* g) const PHYSBAM_OVERRIDE
    {assert(g && !dg && !ddg);tmp.Op(1,x,t,dx);*g=f(tmp);}

//#####################################################################
};

template<class T>
class PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)>:public NONLINEAR_FUNCTION<T(T)>
{
public:
    typedef NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)> F;
    const F& f;
    const KRYLOV_VECTOR_BASE<T> &x,&dx;
    KRYLOV_VECTOR_BASE<T>& tmp,*tmp2;
    KRYLOV_SYSTEM_BASE<T>* system;

    PARAMETRIC_LINE(const F& f,const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& dx,
        KRYLOV_VECTOR_BASE<T>& tmp,KRYLOV_VECTOR_BASE<T>* tmp2=0,KRYLOV_SYSTEM_BASE<T>* system=0)
        :f(f),x(x),dx(dx),tmp(tmp),tmp2(tmp2),system(system)
    {
    }

    virtual ~PARAMETRIC_LINE();
    void Compute(const T t,T* ddg,T* dg,T* g) const PHYSBAM_OVERRIDE;
//#####################################################################
};

}
#endif
