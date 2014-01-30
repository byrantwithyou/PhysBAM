//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_POLYMER_STRESS__
#define __ANALYTIC_POLYMER_STRESS__

#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <boost/function.hpp>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_POLYMER_STRESS
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    virtual ~ANALYTIC_POLYMER_STRESS(){}
    virtual T_MATRIX S(const TV& X,T t) const=0;
    virtual T_MATRIX dSdX(const TV& X,T t,int i) const=0;
    virtual T_MATRIX dSdt(const TV& X,T t) const=0;
    virtual TV divS(const TV& X,T t) const {TV val;for(int i=0;i<TV::m;i++)for(int j=0;j<TV::m;j++)val(j)+=dSdX(X,t,i)(j,i);return val;}
    virtual T_MATRIX F_S(const TV& X,T t) const=0;
    virtual void Test(const TV& X) const
    {
//        RANDOM_NUMBERS<T> rand;
//        TV dX;
//        T e=1e-6,t=rand.Get_Uniform_Number(0,1),dt=rand.Get_Uniform_Number(-e,e);
//        rand.Fill_Uniform(dX,-e,e);
//        SYMMETRIC_MATRIX<T,TV::m> S0=S(X,t),S1=S((X+dX),t);
//        SYMMETRIC_MATRIX<T,TV::m> du0=du(X,t),du1=du((X+dX),t);
//        T erru=((du0+du1)*dX/2-(u1-u0)).Magnitude()/e;
//        LOG::cout<<"analytic velocity diff test "<<erru<<std::endl;
    }
};

template<class TV>
struct ANALYTIC_POLYMER_STRESS_CONST:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    ANALYTIC_POLYMER_STRESS_CONST() {}
    virtual T_MATRIX S(const TV& X,T t) const {return T_MATRIX::Identity_Matrix();}
    virtual T_MATRIX dSdX(const TV& X,T t,int i) const {return T_MATRIX();}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX();}
    virtual T_MATRIX F_S(const TV& X,T t) const {return T_MATRIX();}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_MAGNITUDE:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    T rho;
    ANALYTIC_POLYMER_STRESS_MAGNITUDE(T rho): rho(rho){}
    virtual T_MATRIX S(const TV& X,T t) const {return T_MATRIX::Identity_Matrix()*(rho*X.Magnitude_Squared()+(T)1);}
    virtual T_MATRIX dSdX(const TV& X,T t,int i) const {return T_MATRIX::Identity_Matrix()*(rho*X(i)*(T)2);}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX();}
    virtual T_MATRIX F_S(const TV& X,T t) const {return T_MATRIX();}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_LINEAR:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    T rho;
    TV a;
    ANALYTIC_POLYMER_STRESS_LINEAR(T rho,TV a): rho(rho),a(a){}
    virtual T_MATRIX S(const TV& X,T t) const {return T_MATRIX::Identity_Matrix()*(rho*X.Dot(a)+(T)1);}
    virtual T_MATRIX dSdX(const TV& X,T t,int i) const {return T_MATRIX::Identity_Matrix()*(rho*a(i));}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX();}
    virtual T_MATRIX F_S(const TV& X,T t) const {return T_MATRIX();}
};


}
#endif
