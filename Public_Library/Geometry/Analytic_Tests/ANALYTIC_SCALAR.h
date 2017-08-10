//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_SCALAR__
#define __ANALYTIC_SCALAR__

#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Tools/Symbolics/PROGRAM.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <functional>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_SCALAR
{
    typedef typename TV::SCALAR T;
    virtual ~ANALYTIC_SCALAR(){}
    virtual T f(const TV& X,T t) const=0;
    virtual TV dX(const TV& X,T t) const=0;
    virtual SYMMETRIC_MATRIX<T,TV::m> ddX(const TV& X,T t) const=0;
    virtual T L(const TV& X,T t) const{return ddX(X,t).Trace();}
    virtual T dt(const TV& X,T t) const=0;
    virtual void Test(const TV& X,T t) const
    {
        RANDOM_NUMBERS<T> rand;
        TV dX_;
        T e=1e-6,dt_=rand.Get_Uniform_Number(0,e);
        rand.Fill_Uniform(dX_,-e,e);
        T u0=f(X,t),u1=f(X+dX_,t);
        TV du0=dX(X,t),du1=dX(X+dX_,t);
        T erru=((du0+du1).Dot(dX_)/2-(u1-u0))/e;
        LOG::cout<<"analytic scalar diff test "<<erru<<std::endl;
        SYMMETRIC_MATRIX<T,TV::m> ddu0=ddX(X,t),ddu1=ddX(X+dX_,t);
        T errdu=((ddu0+ddu1)*dX_/2-(du1-du0)).Magnitude()/e;
        LOG::cout<<"analytic scalar hess test "<<errdu<<std::endl;
        T dt0=dt(X,t),dt1=dt(X,t+dt_),u2=f(X,t+dt_);
        T errut=((dt0+dt1)*dt_/2-(u2-u0))/e;
        LOG::cout<<"analytic scalar time tiff test "<<errut<<std::endl;
    }
};

template<class TV,class G>
struct ANALYTIC_SCALAR_CUSTOM:public ANALYTIC_SCALAR<TV>
{
    G g;
    typedef typename TV::SCALAR T;
    ANALYTIC_SCALAR_CUSTOM(G g): g(g) {}
    virtual ~ANALYTIC_SCALAR_CUSTOM(){}
    virtual T f(const TV& X,T t) const {return g(X,t);}
    virtual TV dX(const TV& X,T t) const
    {return AUTO_DIFF<T,TV>(g(AUTO_DIFF<TV,TV>::From_Var(X),t)).dx;}
    virtual SYMMETRIC_MATRIX<T,TV::m> ddX(const TV& X,T t) const
    {return AUTO_HESS<T,TV>(g(AUTO_HESS<TV,TV>::From_Var(X),t)).ddx;}
    virtual T dt(const TV& X,T t) const
    {
        typedef VECTOR<T,1> V1;typedef AUTO_DIFF<T,V1> AD;
        return AD(g(X,AD::From_Var(V1(t),0))).dx.x;
    }
};

template<class TV,class G>
ANALYTIC_SCALAR_CUSTOM<TV,G>* Make_Analytic_Scalar(G g)
{
    return new ANALYTIC_SCALAR_CUSTOM<TV,G>(g);
}

template<class TV>
struct ANALYTIC_SCALAR_PROGRAM:public ANALYTIC_SCALAR<TV>
{
    typedef typename TV::SCALAR T;
    PROGRAM<T> prog;
    mutable PROGRAM_CONTEXT<T> context;

    ANALYTIC_SCALAR_PROGRAM(std::string& str)
    {
        const char* axes[]={"x","y","z"};
        for(int i=0;i<TV::m;i++) prog.var_in.Append(axes[i]);
        prog.var_out.Append("p");
        prog.var_in.Append("t");
        prog.Parse(str.c_str(),false);
        ARRAY<int> out(1);
        out(0)=0;
        for(int j=0;j<TV::m+1;j++)
            prog.Diff(out,j);
        out=IDENTITY_ARRAY<>{TV::m}+1;
        for(int j=0;j<TV::m;j++)
            prog.Diff(out,j);
        prog.Optimize();
        prog.Finalize();
        context.Initialize(prog);
    }
    void Run(const TV& X,T t) const
    {
        for(int i=0;i<TV::m;i++) context.data_in(i)=X(i);
        context.data_in(TV::m)=t;
        prog.Execute(context);
    }
    virtual ~ANALYTIC_SCALAR_PROGRAM(){}
    virtual T f(const TV& X,T t) const
    {
        Run(X,t);
        return context.data_out(0);
    }
    virtual TV dX(const TV& X,T t) const
    {
        Run(X,t);
        TV out;
        for(int i=0;i<TV::m;i++) out(i)=context.data_out(i+1);
        return out;
    }
    virtual SYMMETRIC_MATRIX<T,TV::m> ddX(const TV& X,T t) const
    {
        Run(X,t);
        SYMMETRIC_MATRIX<T,TV::m> out;
        for(int j=0;j<TV::m;j++)
            for(int k=j;k<TV::m;k++)
                out(j,k)=context.data_out(TV::m+2+j+TV::m*k);
        return out;
    }
    virtual T dt(const TV& X,T t) const
    {
        Run(X,t);
        return context.data_out(TV::m+1);
    }
};
}
#endif
