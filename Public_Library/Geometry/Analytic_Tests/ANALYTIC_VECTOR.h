//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_VECTOR__
#define __ANALYTIC_VECTOR__

#include <Core/Arrays/IDENTITY_ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Tools/Symbolics/PROGRAM.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <functional>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_VECTOR
{
    typedef typename TV::SCALAR T;
    virtual ~ANALYTIC_VECTOR(){}
    virtual TV v(const TV& X,T t) const=0;
    virtual TV dt(const TV& X,T t) const=0;
    virtual MATRIX<T,TV::m> dX(const TV& X,T t) const=0;
    virtual SYMMETRIC_TENSOR<T,0,TV::m> ddX(const TV& X,T t) const=0;
    virtual TV L(const TV& X,T t) const {return Contract<1,2>(ddX(X,t));}
    virtual void Test(const TV& X,T t) const
    {
        RANDOM_NUMBERS<T> rand;
        TV dX_;
        T e=1e-6,dt_=rand.Get_Uniform_Number(0,e);
        rand.Fill_Uniform(dX_,-e,e);
        TV u0=v(X,t),u1=v(X+dX_,t);
        MATRIX<T,TV::m> du0=dX(X,t),du1=dX(X+dX_,t);
        T erru=((du0+du1)*dX_/2-(u1-u0)).Magnitude()/e;
        LOG::cout<<"analytic velocity diff test "<<erru<<std::endl;
        SYMMETRIC_TENSOR<T,0,TV::m> ddu0=ddX(X,t),ddu1=ddX(X+dX_,t);
        T errdu=(Contract<2>(ddu0+ddu1,dX_/2)-(du1-du0)).Frobenius_Norm()/e;
        LOG::cout<<"analytic velocity hess test "<<errdu<<std::endl;
        LOG::cout<<"analytic velocity Laplacian test "<<(L(X,t)-Contract<1,2>(ddu0))<<std::endl;
        TV dt0=dt(X,t),dt1=dt(X,t+dt_),u2=v(X,t+dt_);
        T errut=((dt0+dt1)*dt_/2-(u2-u0)).Magnitude()/e;
        LOG::cout<<"analytic velocity time tiff test "<<errut<<std::endl;
    }
};

template<class TV,class G>
struct ANALYTIC_VECTOR_CUSTOM:public ANALYTIC_VECTOR<TV>
{
    G g;
    typedef typename TV::SCALAR T;
    ANALYTIC_VECTOR_CUSTOM(G g): g(g) {}
    virtual ~ANALYTIC_VECTOR_CUSTOM(){}
    virtual TV v(const TV& X,T t) const {return g(X,t);}
    virtual TV dt(const TV& X,T t) const
    {
        typedef VECTOR<T,1> V1;
        return AUTO_DIFF<TV,V1>(g(X,AUTO_DIFF<T,V1>::From_Var(V1(t),0))).dx.Column(0);
    }
    virtual MATRIX<T,TV::m> dX(const TV& X,T t) const
    {return AUTO_DIFF<TV,TV>(g(AUTO_DIFF<TV,TV>::From_Var(X),t)).dx;}
    virtual SYMMETRIC_TENSOR<T,0,TV::m> ddX(const TV& X,T t) const
    {return AUTO_HESS<TV,TV>(g(AUTO_HESS<TV,TV>::From_Var(X),t)).ddx;}
};

template<class TV,class G>
ANALYTIC_VECTOR_CUSTOM<TV,G>* Make_Analytic_Vector(G g)
{
    return new ANALYTIC_VECTOR_CUSTOM<TV,G>(g);
}

template<class TV>
struct ANALYTIC_VECTOR_PROGRAM:public ANALYTIC_VECTOR<TV>
{
    typedef typename TV::SCALAR T;
    PROGRAM<T> prog;
    mutable PROGRAM_CONTEXT<T> context;

    ANALYTIC_VECTOR_PROGRAM(std::string& str)
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
    virtual ~ANALYTIC_VECTOR_PROGRAM(){}
    virtual TV v(const TV& X,T t) const
    {
        Run(X,t);
        TV out;
        for(int i=0;i<TV::m;i++) out(i)=context.data_out(i);
        return out;
    }
    virtual TV dt(const TV& X,T t) const
    {
        Run(X,t);
        TV out;
        for(int i=0;i<TV::m;i++) out(i)=context.data_out(i+TV::m*(TV::m+1));
        return out;
    }
    virtual MATRIX<T,TV::m> dX(const TV& X,T t) const
    {
        Run(X,t);
        MATRIX<T,TV::m> m;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                m(i,j)=context.data_out(j*TV::m+i+TV::m);
        return m;
    }
    virtual SYMMETRIC_TENSOR<T,0,TV::m> ddX(const TV& X,T t) const
    {
        Run(X,t);
        SYMMETRIC_TENSOR<T,0,TV::m> out;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                for(int k=j;k<TV::m;k++)
                    out.x(i)(j,k)=context.data_out(i+TV::m*(TV::m*j+k+TV::m+2));
        return out;
    }
    virtual TV L(const TV& X,T t) const
    {
        Run(X,t);
        TV out;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                out(i)+=context.data_out(i+TV::m*((TV::m+1)*j+TV::m+2));
        return out;
    }
};
}
#endif
