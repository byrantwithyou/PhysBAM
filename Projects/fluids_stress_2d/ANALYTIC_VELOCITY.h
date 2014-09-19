//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_VELOCITY__
#define __ANALYTIC_VELOCITY__

#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <boost/function.hpp>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_VELOCITY
{
    typedef typename TV::SCALAR T;
    virtual ~ANALYTIC_VELOCITY(){}
    virtual TV u(const TV& X,T t) const=0;
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const=0;
    virtual void Test(const TV& X) const
    {
        RANDOM_NUMBERS<T> rand;
        TV dX;
        T e=1e-6,t=rand.Get_Uniform_Number(0,1);
        rand.Fill_Uniform(dX,-e,e);
        TV u0=u(X,t),u1=u((X+dX),t);
        MATRIX<T,TV::m> du0=du(X,t),du1=du((X+dX),t);
        T erru=((du0+du1)*dX/2-(u1-u0)).Magnitude()/e;
        LOG::cout<<"analytic velocity diff test "<<erru<<std::endl;
    }
};

template<class TV>
struct ANALYTIC_VELOCITY_CONST:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV au;
    ANALYTIC_VELOCITY_CONST(TV v): au(v) {}
    virtual TV u(const TV& X,T t) const {return au;}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return MATRIX<T,TV::m>();}
};

template<class TV>
struct ANALYTIC_VELOCITY_ROTATION:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV c;
    typename TV::SPIN w;
    ANALYTIC_VELOCITY_ROTATION(TV cc,typename TV::SPIN ww,T rho): c(cc),w(ww){}
    virtual TV u(const TV& X,T t) const {return w.Cross(X-c);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return MATRIX<T,TV::m>::Cross_Product_Matrix(w);}
};

template<class TV>
struct ANALYTIC_VELOCITY_RAREFACTION:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_VELOCITY_RAREFACTION(){}
    virtual TV u(const TV& X,T t) const {return X.x/(t+1)*TV::Axis_Vector(0);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {MATRIX<T,TV::m> M;M(0,0)=1/(t+1);return M;}
};

template<class TV>
struct ANALYTIC_VELOCITY_AFFINE:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV v0;
    MATRIX<T,TV::m> du0;
    ANALYTIC_VELOCITY_AFFINE(const TV& x0,const TV& v0,const MATRIX<T,TV::m>& du0,T rho): v0(v0-du0*x0),du0(du0) {}
    virtual TV u(const TV& X,T t) const {return du0*X+v0;}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return du0;}
};

template<class TV>
struct ANALYTIC_VELOCITY_QUADRATIC_X:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    T a,b,c;
    ANALYTIC_VELOCITY_QUADRATIC_X(T a,T b,T c): a(a),b(b),c(c) {}
    virtual TV u(const TV& X,T t) const {return TV::Axis_Vector(1)*((a*X.x+b)*X.x+c);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {MATRIX<T,TV::m> M;M(1,0)=2*a*X.x+b;return M;}
};

template<class TV>
struct ANALYTIC_VELOCITY_TRANSLATE:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_VELOCITY<TV>* av;
    TV vel;
    ANALYTIC_VELOCITY_TRANSLATE(ANALYTIC_VELOCITY<TV>* av,const TV& vel): av(av),vel(vel) {}
    ~ANALYTIC_VELOCITY_TRANSLATE() {delete av;}
    virtual TV u(const TV& X,T t) const {return av->u(X-vel*t,t)+vel;}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return av->du(X-vel*t,t);}
};

template<class TV>
struct ANALYTIC_VELOCITY_SHIFT_PRESSURE:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_VELOCITY<TV>* av;
    T dp;
    ANALYTIC_VELOCITY_SHIFT_PRESSURE(ANALYTIC_VELOCITY<TV>* av,T dp): av(av),dp(dp) {}
    ~ANALYTIC_VELOCITY_SHIFT_PRESSURE() {delete av;}
    virtual TV u(const TV& X,T t) const {return av->u(X,t);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return av->du(X,t);}
};

template<class TV>
struct ANALYTIC_VELOCITY_VORTEX:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_VELOCITY_VORTEX(){}
    virtual TV u(const TV& X,T t) const
    {return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y))*exp(-t);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const
    {T c=cos(X.x)*cos(X.y),s=sin(X.x)*sin(X.y);return MATRIX<T,TV::m>(c,s,-s,-c)*exp(-t);}
};
template<class TV>
struct ANALYTIC_VELOCITY_NEST:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    
    ANALYTIC_LEVELSET<TV>* al;
    ARRAY<ANALYTIC_VELOCITY<TV>*> sub_vel;
    ANALYTIC_VELOCITY_NEST(ANALYTIC_LEVELSET<TV>* ls): al(ls) {}
    ANALYTIC_VELOCITY_NEST* Add(ANALYTIC_VELOCITY<TV>* vel){sub_vel.Append(vel);return this;}
    virtual TV u(const TV& X,T t) const
    {int c=0;al->phi(X,t,c);return sub_vel(c)->u(X,t);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const
    {int c=0;al->phi(X,t,c);return sub_vel(c)->du(X,t);}
};
}
#endif
