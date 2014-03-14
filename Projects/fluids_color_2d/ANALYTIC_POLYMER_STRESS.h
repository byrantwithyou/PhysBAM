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
template<class TV>
struct ANALYTIC_POLYMER_STRESS_CONST_DECAY:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    T wi_inv;
    ANALYTIC_POLYMER_STRESS_CONST_DECAY(T wi_inv): wi_inv(wi_inv){}
    virtual T_MATRIX S(const TV& X,T t) const {return T_MATRIX::Identity_Matrix()*(exp(-wi_inv*t)+(T)1);}
    virtual T_MATRIX dSdX(const TV& X,T t,int i) const {return T_MATRIX();}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX::Identity_Matrix()*(-wi_inv*exp(-wi_inv*t));}
    virtual T_MATRIX F_S(const TV& X,T t) const {return T_MATRIX();}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_WAVES:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    T rho;
    ANALYTIC_POLYMER_STRESS_WAVES(T rho): rho(rho){}//This one does not need stress forcing if wi_inv=0
    virtual T_MATRIX S(const TV& X,T t) const {T q=X.Magnitude_Squared()-2*t;T c=cos(q);T s=sin(q);return SYMMETRIC_MATRIX<T,TV::m>(rho*s+1,rho*c,1-rho*s);}
    virtual T_MATRIX dSdX(const TV& X,T t,int i) const {T q=X.Magnitude_Squared()-2*t;T c=cos(q);T s=sin(q);return SYMMETRIC_MATRIX<T,TV::m>(c,-s,-c)*2*X(i)*rho;}
    virtual T_MATRIX dSdt(const TV& X,T t) const {T q=X.Magnitude_Squared()-2*t;T c=cos(q);T s=sin(q);return SYMMETRIC_MATRIX<T,TV::m>(-c,s,c)*2*rho;}
    virtual T_MATRIX F_S(const TV& X,T t) const {return T_MATRIX();}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_QUADRATIC:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    T rho;
    VECTOR<VECTOR<VECTOR<T,3>,3>,3> a,b,c;
    ANALYTIC_POLYMER_STRESS_QUADRATIC(T rho): rho(rho){
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){
                a(i)(j)(0)=((i*2741+j*179) % 263)/263.0 -.5;
                a(i)(j)(1)=((i*1259+j*643) % 83)/83.0 -.5;
                a(i)(j)(2)=((i*317+j*421) % 113)/113.0 -.5;
                b(i)(j)(0)=((i*1069+j*1381) % 173)/173.0 -.5;
                b(i)(j)(1)=((i*1489+j*1699) % 179)/179.0 -.5;
                b(i)(j)(2)=((i*1601+j*1811) % 43)/43.0 -.5;
                c(i)(j)(0)=((i*2017+j*2027) % 89)/89.0 -.5;
                c(i)(j)(1)=((i*2239+j*2029) % 97)/97.0 -.5;
                c(i)(j)(2)=((i*2243+j*2267) % 107)/107.0 -.5;
            }
    }//This one only works with wi_inv=0
    virtual T_MATRIX S(const TV& X,T t) const {
        T a11(0),a12(0),a22(0);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++){
                    a11+=a(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k);
                    a12+=b(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k);
                    a22+=c(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k);
                }
        return SYMMETRIC_MATRIX<T,TV::m>(a11,a12,a22);}
    virtual T_MATRIX dSdX(const TV& X,T t,int dim) const {
        T a11(0),a12(0),a22(0);        
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++){
                    T ij=(dim==0)?i:j;
                    if(ij>0){
                    a11+=a(i)(j)(k)*ij*pow(X.x,i-1+dim)*pow(X.y,j-dim)*pow(t,k);
                    a12+=b(i)(j)(k)*ij*pow(X.x,i-1+dim)*pow(X.y,j-dim)*pow(t,k);
                    a22+=c(i)(j)(k)*ij*pow(X.x,i-1+dim)*pow(X.y,j-dim)*pow(t,k);
                    }}
        return SYMMETRIC_MATRIX<T,TV::m>(a11,a12,a22);}
    virtual T_MATRIX dSdt(const TV& X,T t) const {
        T a11(0),a12(0),a22(0);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=1;k<3;k++){
                    a11+=a(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                    a12+=b(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                    a22+=c(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                }
        return SYMMETRIC_MATRIX<T,TV::m>(a11,a12,a22);}
    virtual T_MATRIX F_S(const TV& X,T t) const {return T_MATRIX();}
};


}
#endif
