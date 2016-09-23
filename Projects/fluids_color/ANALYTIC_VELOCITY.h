//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_VELOCITY__
#define __ANALYTIC_VELOCITY__

#include <Core/Arrays/IDENTITY_ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Symbolics/PROGRAM.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
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
    virtual void Test(const TV& X) const
    {
        RANDOM_NUMBERS<T> rand;
        TV dX;
        T e=1e-6,t=rand.Get_Uniform_Number(0,1),dt=rand.Get_Uniform_Number(0,e);
        rand.Fill_Uniform(dX,-e,e);
        TV u0=u(X,t),u1=u(X+dX,t);
        MATRIX<T,TV::m> du0=du(X,t),du1=du(X+dX,t);
        T erru=((du0+du1)*dX/2-(u1-u0)).Magnitude()/e;
        LOG::cout<<"analytic velocity diff test "<<erru<<std::endl;
        SYMMETRIC_TENSOR<T,0,TV::m> ddu0=ddu(X,t),ddu1=ddu(X+dX,t);
        T errdu=(Contract<2>(ddu0+ddu1,dX/2)-(du1-du0)).Frobenius_Norm()/e;
        LOG::cout<<"analytic velocity hess test "<<errdu<<std::endl;
        LOG::cout<<"analytic velocity Laplacian test "<<(Lu(X,t)-Contract<1,2>(ddu0))<<std::endl;
        TV dudt0=dudt(X,t),dudt1=dudt(X,t+dt),u2=u(X,t+dt);
        T errut=((dudt0+dudt1)*dt/2-(u2-u0)).Magnitude()/e;
        LOG::cout<<"analytic velocity time tiff test "<<errut<<std::endl;
    }
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

template<class TV>
struct ANALYTIC_VELOCITY_PROGRAM:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    PROGRAM<T> prog;
    mutable PROGRAM_CONTEXT<T> context;

    ANALYTIC_VELOCITY_PROGRAM(std::string& str)
    {
        const char* axes[]={"x","y","z"};
        const char* vel[]={"u","v","w"};
        for(int i=0;i<TV::m;i++){
            prog.var_in.Append(axes[i]);
            prog.var_out.Append(vel[i]);}
        prog.var_in.Append("t");
        prog.Parse(str.c_str(),false);
        ARRAY<int> out{IDENTITY_ARRAY<>{TV::m}};
        for(int j=0;j<prog.var_in.m;j++)
            prog.Diff(out,j);
        out.Resize(TV::m*TV::m);
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                out(i*TV::m+j)=i*TV::m+TV::m+j;
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
    virtual ~ANALYTIC_VELOCITY_PROGRAM(){}
    virtual TV u(const TV& X,T t) const
    {
        Run(X,t);
        TV out;
        for(int i=0;i<TV::m;i++) out(i)=context.data_out(i);
        return out;
    }
    virtual TV dudt(const TV& X,T t) const
    {
        Run(X,t);
        TV out;
        for(int i=0;i<TV::m;i++) out(i)=context.data_out(i+TV::m*(TV::m+1));
        return out;
    }
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const
    {
        Run(X,t);
        MATRIX<T,TV::m> m;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                m(i,j)=context.data_out(j*TV::m+i+TV::m);
        return m;
    }
    virtual SYMMETRIC_TENSOR<T,0,TV::m> ddu(const TV& X,T t) const
    {
        Run(X,t);
        SYMMETRIC_TENSOR<T,0,TV::m> out;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                for(int k=j;k<TV::m;k++)
                    out.x(i)(j,k)=context.data_out(i+TV::m*(TV::m*j+k+TV::m+2));
        return out;
    }
    virtual TV Lu(const TV& X,T t) const
    {
        Run(X,t);
        TV out;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                out(i)+=context.data_out(i+TV::m*((TV::m+1)*j+TV::m+2));
        return out;
    }
};
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

template<class TV>
struct ANALYTIC_PRESSURE_PROGRAM:public ANALYTIC_PRESSURE<TV>
{
    typedef typename TV::SCALAR T;
    PROGRAM<T> prog;
    mutable PROGRAM_CONTEXT<T> context;

    ANALYTIC_PRESSURE_PROGRAM(std::string& str)
    {
        const char* axes[]={"x","y","z"};
        for(int i=0;i<TV::m;i++) prog.var_in.Append(axes[i]);
        prog.var_in.Append("t");
        prog.var_out.Append("p");
        prog.Parse(str.c_str(),false);
        ARRAY<int> out(1);
        out(0)=0;
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
    virtual ~ANALYTIC_PRESSURE_PROGRAM(){}
    virtual T p(const TV& X,T t) const
    {
        Run(X,t);
        return context.data_out(0);
    }
    virtual TV dp(const TV& X,T t) const
    {
        Run(X,t);
        TV out;
        for(int i=0;i<TV::m;i++) out(i)=context.data_out(i+1);
        return out;
    }
};
}
#endif
