//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_VELOCITY__
#define __ANALYTIC_VELOCITY__

#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <functional>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_VELOCITY
{
    typedef typename TV::SCALAR T;
    virtual ~ANALYTIC_VELOCITY(){}
    virtual TV u(const TV& X,T t) const=0;
    virtual TV dudt(const TV& X,T t) const=0;
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const=0;
    virtual SYMMETRIC_TENSOR<T,0,TV::m> ddu(const TV& X,T t) const=0;
    virtual TV Lu(const TV& X,T t) const {return Contract<1,2>(ddu(X,t));}
};

template<class TV>
struct ANALYTIC_PRESSURE
{
    typedef typename TV::SCALAR T;
    virtual ~ANALYTIC_PRESSURE(){}
    virtual T p(const TV& X,T t) const=0;
    virtual TV dp(const TV& X,T t) const=0;
};

template<class TV,class F>
struct ANALYTIC_VELOCITY_CUSTOM:public ANALYTIC_VELOCITY<TV>
{
    F f;
    typedef typename TV::SCALAR T;
    ANALYTIC_VELOCITY_CUSTOM(F f): f(f) {}
    virtual ~ANALYTIC_VELOCITY_CUSTOM(){}
    virtual TV u(const TV& X,T t) const {return f(X,t);}
    virtual TV dudt(const TV& X,T t) const
    {
        typedef DIFF_LAYOUT<T,-1> LAYOUT;
        TV ret;
        Get<0>(ret,(Diff_From_Const<LAYOUT>(TV())+f(X,Diff_From_Var<LAYOUT,0>(t))).dx);
        return ret;
    }
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const
    {
        typedef DIFF_LAYOUT<T,TV::m> LAYOUT;
        MATRIX<T,TV::m> ret;
        Get<0>(ret,(Diff_From_Const<LAYOUT>(TV())+f(Diff_From_Var<LAYOUT,0>(X),t)).dx);
        return ret;
    }
    virtual SYMMETRIC_TENSOR<T,0,TV::m> ddu(const TV& X,T t) const
    {
        typedef DIFF_LAYOUT<T,TV::m> LAYOUT;
        SYMMETRIC_TENSOR<T,0,TV::m> ret;
        Get<0,0>(ret,(Hess_From_Const<LAYOUT>(TV())+f(Hess_From_Var<LAYOUT,0>(X),t)).ddx);
        return ret;
    }
};

template<class TV,class F>
ANALYTIC_VELOCITY_CUSTOM<TV,F>* Make_Velocity(F f)
{
    return new ANALYTIC_VELOCITY_CUSTOM<TV,F>(f);
}

template<class TV,class F>
struct ANALYTIC_PRESSURE_CUSTOM:public ANALYTIC_PRESSURE<TV>
{
    F f;
    typedef typename TV::SCALAR T;
    ANALYTIC_PRESSURE_CUSTOM(F f): f(f) {}
    virtual ~ANALYTIC_PRESSURE_CUSTOM(){}
    virtual T p(const TV& X,T t) const {return f(X,t);}
    virtual TV dp(const TV& X,T t) const
    {
        typedef DIFF_LAYOUT<T,TV::m> LAYOUT;
        TV ret;
        Get<0>(ret,(Diff_From_Const<LAYOUT>(T())+f(Diff_From_Var<LAYOUT,0>(X),t)).dx);
        return ret;
    }
};

template<class TV,class F>
ANALYTIC_PRESSURE_CUSTOM<TV,F>* Make_Pressure(F f)
{
    return new ANALYTIC_PRESSURE_CUSTOM<TV,F>(f);
}
}
#endif
