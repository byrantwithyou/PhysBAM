//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_POLYMER_STRESS__
#define __ANALYTIC_POLYMER_STRESS__

#include <Tools/Log/LOG.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/IDENTITY_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <functional>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_POLYMER_STRESS
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    virtual ~ANALYTIC_POLYMER_STRESS(){}
    virtual T_MATRIX S(const TV& X,T t) const=0;
    virtual T_TENSOR dSdX(const TV& X,T t) const=0;
    virtual T_MATRIX dSdt(const TV& X,T t) const=0;
    virtual void Test(const TV& X) const
    {
        RANDOM_NUMBERS<T> rand;
        TV dX;
        T e=1e-6,t=rand.Get_Uniform_Number(0,1),dt=rand.Get_Uniform_Number(-e,e);
        rand.Fill_Uniform(dX,-e,e);
        SYMMETRIC_MATRIX<T,TV::m> S0=S(X,t),S1=S(X+dX,t);
        SYMMETRIC_MATRIX<T,TV::m> dS=Contract<2>(dSdX(X,t)+dSdX(X+dX,t),dX)/2;
        T erru=(dS-(S1-S0)).Frobenius_Norm()/e;
        LOG::cout<<"analytic stress x diff test "<<erru<<std::endl;
        SYMMETRIC_MATRIX<T,TV::m> R0=S(X,t),R1=S(X,t+dt);
        SYMMETRIC_MATRIX<T,TV::m> dR=(dSdt(X,t)+dSdt(X,t+dt))/2*dt;
        T errt=(dR-(R1-R0)).Frobenius_Norm()/e;
        LOG::cout<<"analytic stress t diff test "<<errt<<std::endl;
    }
};

template<class TV>
struct ANALYTIC_POLYMER_STRESS_CONST:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    ANALYTIC_POLYMER_STRESS_CONST() {}
    virtual T_MATRIX S(const TV& X,T t) const {return T_MATRIX::Identity_Matrix();}
    virtual T_TENSOR dSdX(const TV& X,T t) const {return T_TENSOR();}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX();}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_MAGNITUDE:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    T rho;
    ANALYTIC_POLYMER_STRESS_MAGNITUDE(T rho): rho(rho){}
    virtual T_MATRIX S(const TV& X,T t) const {return T_MATRIX::Identity_Matrix()*(rho*X.Magnitude_Squared()+(T)1);}
    virtual T_TENSOR dSdX(const TV& X,T t) const {return T_TENSOR()+Tensor_Product<2>(IDENTITY_MATRIX<T,TV::m>(),rho*X*(T)2);}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX();}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_LINEAR:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    T rho;
    TV a;
    ANALYTIC_POLYMER_STRESS_LINEAR(T rho,TV a): rho(rho),a(a){}
    virtual T_MATRIX S(const TV& X,T t) const {return T_MATRIX::Identity_Matrix()*(rho*X.Dot(a)+(T)1);}
    virtual T_TENSOR dSdX(const TV& X,T t) const {return T_TENSOR()+Tensor_Product<2>(IDENTITY_MATRIX<T,TV::m>(),rho*a);}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX();}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_CONST_DECAY:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    T wi_inv;
    ANALYTIC_POLYMER_STRESS_CONST_DECAY(T wi_inv): wi_inv(wi_inv){}
    virtual T_MATRIX S(const TV& X,T t) const {return T_MATRIX::Identity_Matrix()*(exp(-wi_inv*t)+(T)1);}
    virtual T_TENSOR dSdX(const TV& X,T t) const {return T_TENSOR();}
    virtual T_MATRIX dSdt(const TV& X,T t) const {return T_MATRIX::Identity_Matrix()*(-wi_inv*exp(-wi_inv*t));}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_WAVES:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    T rho;
    ANALYTIC_POLYMER_STRESS_WAVES(T rho): rho(rho){}//This one does not need stress forcing if wi_inv=0
    virtual T_MATRIX S(const TV& X,T t) const {T q=X.Magnitude_Squared()-2*t;T c=cos(q);T s=sin(q);return SYMMETRIC_MATRIX<T,TV::m>(rho*s+1,rho*c,1-rho*s);}
    virtual T_TENSOR dSdX(const TV& X,T t) const {T q=X.Magnitude_Squared()-2*t;T c=cos(q);T s=sin(q);return Tensor_Product<2>(SYMMETRIC_MATRIX<T,TV::m>(c,-s,-c),X*rho*2);}
    virtual T_MATRIX dSdt(const TV& X,T t) const {T q=X.Magnitude_Squared()-2*t;T c=cos(q);T s=sin(q);return SYMMETRIC_MATRIX<T,TV::m>(-c,s,c)*2*rho;}
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_PERIODIC:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    RANGE<TV> domain;
    ANALYTIC_POLYMER_STRESS_PERIODIC(const RANGE<TV>& domain): domain(domain){}
    virtual T_MATRIX S(const TV& X,T t) const
    {
        TV Z=(X-domain.min_corner)/domain.Edge_Lengths()*(T)(2*pi);
        TV W=(cos(Z)+sin(Z))*(1+exp(-t));
        return SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(W)+1;
    }
    virtual T_TENSOR dSdX(const TV& X,T t) const
    {
        TV Z=(X-domain.min_corner)/domain.Edge_Lengths()*(T)(2*pi);
        DIAGONAL_MATRIX<T,TV::m> dZ((T)(2*pi)/domain.Edge_Lengths());
        TV W=(cos(Z)+sin(Z))*(1+exp(-t));
        DIAGONAL_MATRIX<T,TV::m> dWdZ((cos(Z)-sin(Z))*(1+exp(-t)));
        return T_TENSOR()+Symmetric_Tensor_Product<0,1>(SYMMETRIC_MATRIX<T,TV::m>(dWdZ*dZ),W);
    }
    virtual T_MATRIX dSdt(const TV& X,T t) const
    {
        TV Z=(X-domain.min_corner)/domain.Edge_Lengths()*(T)(2*pi);
        TV W=(cos(Z)+sin(Z))*(1+exp(-t));
        TV dW=-(cos(Z)+sin(Z))*exp(-t);
        return SYMMETRIC_MATRIX<T,TV::m>::Symmetric_Outer_Product(W,dW);
    }
};
template<class TV>
struct ANALYTIC_POLYMER_STRESS_QUADRATIC:public ANALYTIC_POLYMER_STRESS<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
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
    virtual T_TENSOR dSdX(const TV& X,T t) const {
        T_TENSOR ten;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                for(int k=0;k<3;k++){
                    for(int dim=0;dim<TV::m;dim++){
                        T ij=(dim==0)?i:j;
                        if(ij>0){
                            ten.x(dim)(0,0)+=a(i)(j)(k)*ij*pow(X.x,i-1+dim)*pow(X.y,j-dim)*pow(t,k);
                            ten.x(dim)(0,1)+=b(i)(j)(k)*ij*pow(X.x,i-1+dim)*pow(X.y,j-dim)*pow(t,k);
                            ten.x(dim)(1,1)+=c(i)(j)(k)*ij*pow(X.x,i-1+dim)*pow(X.y,j-dim)*pow(t,k);
                        }}}
        return ten;}
    virtual T_MATRIX dSdt(const TV& X,T t) const {
        T a11(0),a12(0),a22(0);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=1;k<3;k++){
                    a11+=a(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                    a12+=b(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                    a22+=c(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                }
        return SYMMETRIC_MATRIX<T,TV::m>(a11,a12,a22);}
};


}
#endif
