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
    virtual TV divS(const TV& X,T t) const
    {TV val;for(int i=0;i<TV::m;i++) val+=dSdX(X,t,i).Row(i);return val;}
    virtual void Test(const TV& X) const
    {
        RANDOM_NUMBERS<T> rand;
        TV dX;
        T e=1e-6,t=rand.Get_Uniform_Number(e,1),dt=rand.Get_Uniform_Number(-e,e);
        rand.Fill_Uniform(dX,-e,e);
        SYMMETRIC_MATRIX<T,TV::m> S0=S(X,t),S1=S(X+dX,t),S2=S(X,t+dt),dS2a=S2-S0,dS2b=(dt/2)*(dSdt(X,t)+dSdt(X,t+dt));
        T a2=dS2a.Frobenius_Norm(),b2=dS2b.Frobenius_Norm(),d2=(dS2a-dS2b).Frobenius_Norm();
        LOG::printf("stress diff t %g %g rel %g\n",a2,b2,d2/max(a2,b2,(T)1e-30));
        SYMMETRIC_MATRIX<T,TV::m> dS1a=S1-S0,dS1b;
        for(int i=0;i<TV::m;i++) dS1b+=(dSdX(X,t,i)+dSdX(X+dX,t,i))*(dX(i)/2);
        T a1=dS1a.Frobenius_Norm(),b1=dS1b.Frobenius_Norm(),d1=(dS1a-dS1b).Frobenius_Norm();
        LOG::printf("stress diff X %g %g rel %g\n",a1,b1,d1/max(a1,b1,(T)1e-30));
    }
};

template<class TV>
struct ANALYTIC_POLYMER_STRESS_CONST:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    SYMMETRIC_MATRIX<T,TV::m> s;
    ANALYTIC_POLYMER_STRESS_CONST(const SYMMETRIC_MATRIX<T,TV::m>& s=T_MATRIX::Identity_Matrix()) {}
    virtual T_MATRIX S(const TV& X,T t) const {return s;}
    virtual T_MATRIX dSdX(const TV& X,T t,int i) const {return T_MATRIX();}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX();}
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
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_TRANSLATE:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_POLYMER_STRESS<TV>* aps;
    TV vel;
    ANALYTIC_POLYMER_STRESS_TRANSLATE(ANALYTIC_POLYMER_STRESS<TV>* aps,const TV& vel): aps(aps),vel(vel) {}
    ~ANALYTIC_POLYMER_STRESS_TRANSLATE() {delete aps;}
    virtual SYMMETRIC_MATRIX<T,TV::m> S(const TV& X,T t) const {return aps->S(X-vel*t,0);}
    virtual SYMMETRIC_MATRIX<T,TV::m> dSdX(const TV& X,T t,int i) const {return aps->dSdX(X-vel*t,0,i);}
    virtual SYMMETRIC_MATRIX<T,TV::m> dSdt(const TV& X,T t) const
    {
        TV Z=X-vel*t;
        SYMMETRIC_MATRIX<T,TV::m> m;
        for(int i=0;i<TV::m;i++) m-=aps->dSdX(Z,0,i)*vel(i);
        return m;
    }
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_ROTATION:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_POLYMER_STRESS<TV>* aps;
    TV c;
    typename TV::SPIN w;
    ANALYTIC_POLYMER_STRESS_ROTATION(ANALYTIC_POLYMER_STRESS<TV>* aps,TV cc,typename TV::SPIN ww): aps(aps),c(cc),w(ww) {}
    ~ANALYTIC_POLYMER_STRESS_ROTATION() {delete aps;}
    virtual SYMMETRIC_MATRIX<T,TV::m> S(const TV& X,T t) const
    {
        TV X0=ROTATION<TV>::From_Rotation_Vector(-t*w).Rotate(X-c)+c;
        return aps->S(X0,0);
    }
    virtual SYMMETRIC_MATRIX<T,TV::m> dSdX(const TV& X,T t,int a) const
    {
        ROTATION<TV> rot=ROTATION<TV>::From_Rotation_Vector(-t*w);
        TV X0=rot.Rotate(X-c)+c;
        MATRIX<T,TV::m> rot_mat=rot.Rotation_Matrix();
        SYMMETRIC_MATRIX<T,TV::m> m;
        for(int i=0;i<TV::m;i++) m+=aps->dSdX(X0,0,i)*rot_mat(i,a);
        return m;
    }
    virtual SYMMETRIC_MATRIX<T,TV::m> dSdt(const TV& X,T t) const
    {
        TV Z=ROTATION<TV>::From_Rotation_Vector(-t*w).Rotate(X-c),X0=Z+c;
        SYMMETRIC_MATRIX<T,TV::m> m;
        TV dX0=-w.Cross(Z);
        for(int i=0;i<TV::m;i++) m+=aps->dSdX(X0,0,i)*dX0(i);
        return m;
    }
};
}
#endif
