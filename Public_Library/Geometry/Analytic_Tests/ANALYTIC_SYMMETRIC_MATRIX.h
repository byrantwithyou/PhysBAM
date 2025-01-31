//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_SYMMETRIC_MATRIX__
#define __ANALYTIC_SYMMETRIC_MATRIX__

#include <Core/Log/LOG.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/IDENTITY_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <functional>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_SYMMETRIC_MATRIX
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    virtual ~ANALYTIC_SYMMETRIC_MATRIX(){}
    virtual T_MATRIX S(const TV& X,T t) const=0;
    virtual T_TENSOR dX(const TV& X,T t) const=0;
    virtual T_MATRIX dt(const TV& X,T t) const=0;
    virtual void Test(const TV& X,T t) const
    {
        RANDOM_NUMBERS<T> rand;
        TV dX_;
        T e=1e-6,dt_=rand.Get_Uniform_Number(-e,e);
        rand.Fill_Uniform(dX_,-e,e);
        SYMMETRIC_MATRIX<T,TV::m> S0=S(X,t),S1=S(X+dX_,t);
        SYMMETRIC_MATRIX<T,TV::m> ds=Contract<2>(dX(X,t)+dX(X+dX_,t),dX_)/2;
        T erru=(ds-(S1-S0)).Frobenius_Norm()/e;
        LOG::cout<<"analytic stress x diff test "<<erru<<std::endl;
        SYMMETRIC_MATRIX<T,TV::m> R0=S(X,t),R1=S(X,t+dt_);
        SYMMETRIC_MATRIX<T,TV::m> dR=(dt(X,t)+dt(X,t+dt_))/2*dt_;
        T errt=(dR-(R1-R0)).Frobenius_Norm()/e;
        LOG::cout<<"analytic stress t diff test "<<errt<<std::endl;
    }
};

template<class TV,class G>
struct ANALYTIC_SYMMETRIC_MATRIX_CUSTOM:public ANALYTIC_SYMMETRIC_MATRIX<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    G g;
    ANALYTIC_SYMMETRIC_MATRIX_CUSTOM(G g): g(g) {}
    virtual ~ANALYTIC_SYMMETRIC_MATRIX_CUSTOM(){}
    virtual T_MATRIX S(const TV& X,T t) const {T_MATRIX m;Fill_From(m,g(X,t));return m;}
    virtual T_MATRIX dt(const TV& X,T t) const
    {
        typedef VECTOR<T,1> V1;
        return AUTO_DIFF<T_MATRIX,V1>(g(X,AUTO_DIFF<T,V1>::From_Var(V1(t),0))).dx.x.x;
    }
    virtual T_TENSOR dX(const TV& X,T t) const
    {return AUTO_DIFF<T_MATRIX,TV>(g(AUTO_DIFF<TV,TV>::From_Var(X),t)).dx;}
};

template<class TV,class G>
ANALYTIC_SYMMETRIC_MATRIX_CUSTOM<TV,G>* Make_Analytic_Symmetric_Matrix(G g)
{
    return new ANALYTIC_SYMMETRIC_MATRIX_CUSTOM<TV,G>(g);
}

template<class TV>
struct ANALYTIC_SYMMETRIC_MATRIX_QUADRATIC:public ANALYTIC_SYMMETRIC_MATRIX<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    T rho;
    VECTOR<VECTOR<VECTOR<T,3>,3>,3> a,b,c;
    ANALYTIC_SYMMETRIC_MATRIX_QUADRATIC(T rho): rho(rho){
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
    virtual T_TENSOR dX(const TV& X,T t) const {
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
    virtual T_MATRIX dt(const TV& X,T t) const {
        T a11(0),a12(0),a22(0);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=1;k<3;k++){
                    a11+=a(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                    a12+=b(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                    a22+=c(i)(j)(k)*pow(X.x,i)*pow(X.y,j)*pow(t,k-1)*k;
                }
        return SYMMETRIC_MATRIX<T,TV::m>(a11,a12,a22);}
};

template<class TV>
struct ANALYTIC_SYMMETRIC_MATRIX_PROGRAM:public ANALYTIC_SYMMETRIC_MATRIX<TV>
{
    typedef typename TV::SCALAR T;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MATRIX;
    typedef SYMMETRIC_TENSOR<T,2,TV::m,TV::m> T_TENSOR;
    PROGRAM<T> prog;
    mutable PROGRAM_CONTEXT<T> context;
    enum {n=TV::m*(TV::m+1)/2};

    ANALYTIC_SYMMETRIC_MATRIX_PROGRAM(std::string& str)
    {
        const char* axes[]={"x","y","z"};
        for(int i=0;i<TV::m;i++) prog.var_in.Append(axes[i]);
        prog.var_in.Append("t");
        char buff[10]="s00";
        for(int j=0;j<TV::m;j++)
            for(int i=j;i<TV::m;i++){
                buff[1]=i+'0';
                buff[2]=j+'0';
                LOG::cout<<"var "<<buff<<std::endl;
                prog.var_out.Append(buff);}
        prog.Parse(str.c_str(),false);
        ARRAY<int> out(IDENTITY_ARRAY<>(prog.var_out.m));
        for(int j=0;j<prog.var_in.m;j++)
            prog.Diff(out,j);
        prog.Optimize();
        prog.Finalize();
        context.Initialize(prog);
        LOG::cout<<prog.var_in<<"   "<<prog.var_out<<std::endl;
    }
    void Run(const TV& X,T t) const
    {
        for(int i=0;i<TV::m;i++) context.data_in(i)=X(i);
        context.data_in(TV::m)=t;
        prog.Execute(context);
    }
    virtual ~ANALYTIC_SYMMETRIC_MATRIX_PROGRAM(){}

    virtual T_MATRIX S(const TV& X,T t) const
    {
        Run(X,t);
        T_MATRIX out;
        for(int i=0;i<n;i++) out.array[i]=context.data_out(i);
        return out;
    }
    virtual T_MATRIX dt(const TV& X,T t) const
    {
        Run(X,t);
        T_MATRIX out;
        for(int i=0;i<n;i++) out.array[i]=context.data_out(i+(TV::m+1)*n);
        return out;
    }
    virtual T_TENSOR dX(const TV& X,T t) const
    {
        Run(X,t);
        T_TENSOR m;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<n;j++)
                m.x(i).array[j]=context.data_out(n*(i+1)+j);
        return m;
    }
};
}
#endif
