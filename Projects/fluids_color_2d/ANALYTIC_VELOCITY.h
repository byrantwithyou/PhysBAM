//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_VELOCITY__
#define __ANALYTIC_VELOCITY__

#include <PhysBAM_Tools/Matrices/MATRIX.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_VELOCITY
{
    typedef typename TV::SCALAR T;
    virtual ~ANALYTIC_VELOCITY(){}
    virtual TV u(const TV& X,T t) const=0;
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const=0;
    virtual T p(const TV& X,T t) const=0;
    virtual TV F(const TV& X,T t) const=0;
};

template<class TV>
struct ANALYTIC_VELOCITY_CONST:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV au;
    T const_p;
    ANALYTIC_VELOCITY_CONST(TV v): au(v),const_p(0) {}
    virtual TV u(const TV& X,T t) const {return au;}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return MATRIX<T,TV::m>();}
    virtual T p(const TV& X,T t) const {return const_p;}
    virtual TV F(const TV& X,T t) const {return TV();}
};

template<class TV>
struct ANALYTIC_VELOCITY_ROTATION:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV c;
    typename TV::SPIN w;
    T rho;
    ANALYTIC_VELOCITY_ROTATION(TV cc,typename TV::SPIN ww,T rho): c(cc),w(ww),rho(rho){}
    virtual TV u(const TV& X,T t) const {return w.Cross(X-c);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return MATRIX<T,TV::m>::Cross_Product_Matrix(w);}
    virtual T p(const TV& X,T t) const {return (T).5*rho*u(X,t).Magnitude_Squared();}
    virtual TV F(const TV& X,T t) const {return TV();}
};

template<class TV>
struct ANALYTIC_VELOCITY_RAREFACTION:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_VELOCITY_RAREFACTION(){}
    virtual TV u(const TV& X,T t) const {return X.x/(t+1)*TV::Axis_Vector(0);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {MATRIX<T,TV::m> M;M(0,0)=1/(t+1);return M;}
    virtual T p(const TV& X,T t) const {return 0;}
    virtual TV F(const TV& X,T t) const {return TV();}
};

template<class TV>
struct ANALYTIC_VELOCITY_AFFINE:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV v0;
    MATRIX<T,TV::m> du0;
    T rho;
    ANALYTIC_VELOCITY_AFFINE(const TV& x0,const TV& v0,const MATRIX<T,TV::m>& du0,T rho): v0(v0-du0*x0),du0(du0),rho(rho) {}
    virtual TV u(const TV& X,T t) const {return du0*X+v0;}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return du0;}
    virtual T p(const TV& X,T t) const {return 0;}
    virtual TV F(const TV& X,T t) const {return rho*du0*(du0*X+v0);}
};

template<class TV>
struct ANALYTIC_VELOCITY_QUADRATIC_X:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    T a,b,c,mu;
    ANALYTIC_VELOCITY_QUADRATIC_X(T a,T b,T c,T mu): a(a),b(b),c(c),mu(mu) {}
    virtual TV u(const TV& X,T t) const {return TV::Axis_Vector(1)*((a*X.x+b)*X.x+c);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {MATRIX<T,TV::m> M;M(1,0)=2*a*X.x+b;return M;}
    virtual T p(const TV& X,T t) const {return 0;}
    virtual TV F(const TV& X,T t) const {return -2*a*mu*TV::Axis_Vector(1);}
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
    virtual T p(const TV& X,T t) const {return av->p(X-vel*t,t);}
    virtual TV F(const TV& X,T t) const {return av->F(X-vel*t,t);}
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
    virtual T p(const TV& X,T t) const {return av->p(X,t)+dp;}
    virtual TV F(const TV& X,T t) const {return av->F(X,t);}
};

template<class TV>
struct ANALYTIC_VELOCITY_VORTEX:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    T nu,rho;
    ANALYTIC_VELOCITY_VORTEX(T mu,T rho): nu(mu/rho),rho(rho){}
    virtual TV u(const TV& X,T t) const
    {return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y))*exp(-2*nu*t);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const
    {T c=cos(X.x)*cos(X.y),s=sin(X.x)*sin(X.y);return MATRIX<T,TV::m>(c,s,-s,-c)*exp(-2*nu*t);}
    virtual T p(const TV& X,T t) const {return (T).25*rho*(cos(2*X.x)+cos(2*X.y))*exp(-4*nu*t);}
    virtual TV F(const TV& X,T t) const {return TV();}
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
    virtual T p(const TV& X,T t) const {int c=0;al->phi(X,t,c);return sub_vel(c)->p(X,t);}
    virtual TV F(const TV& X,T t) const {int c=0;al->phi(X,t,c);return sub_vel(c)->F(X,t);}
};

template<class TV> //Reverses direction with cos(t), adds constant velocity
struct ANALYTIC_VELOCITY_VORTEX_NEW:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    T nu,rho,l; //lambda is frequency of reversing, 0 doesn't reverse at all
    TV au; //shift vector
    T const_p;
    ANALYTIC_VELOCITY_VORTEX_NEW(T mu,T rho,T lambda,TV v): nu(mu/rho),rho(rho),l(lambda),au(v),const_p(0){}
    virtual TV u(const TV& X,T t) const
    {return TV(sin(X.x)*cos(X.y)*cos(l*t),-cos(X.x)*sin(X.y)*cos(l*t))+au;}
    virtual MATRIX<T,2> du(const TV& X,T t) const
    {T c=cos(X.x)*cos(X.y)*cos(l*t),s=sin(X.x)*sin(X.y)*cos(l*t);return MATRIX<T,2>(c,s,-s,-c);}
    virtual T p(const TV& X,T t) const {return const_p+(T).25*rho*cos(l*t)*cos(l*t)*(cos(2*X.x)+cos(2*X.y));}
    virtual TV F(const TV& X,T t) const {return TV(sin(X.x)*cos(X.y)*((T)2*nu*cos(l*t)-l*sin(l*t))+cos(l*t)*cos(l*t)*(au(0)*cos(X.x)*cos(X.y)-au(1)*sin(X.x)*sin(X.y))-au(0)*l*sin(l*t),-cos(X.x)*sin(X.y)*((T)2*nu*cos(l*t)-l*sin(l*t))+cos(l*t)*cos(l*t)*(au(0)*sin(X.x)*sin(X.y)-au(1)*cos(X.x)*cos(X.y))-au(1)*l*sin(l*t));}
};

template<class TV>
struct ANALYTIC_VELOCITY_GRADED_ROTATION:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    T r,nu,rho,scale;
    ANALYTIC_VELOCITY_GRADED_ROTATION(T rad,T mu,T rho,T scale): r(rad),nu(mu/rho),rho(rho),scale(scale){}
    virtual TV u(const TV& X,T t) const {T q=X.Magnitude_Squared()-sqr(r);return sqr(q)*scale*X.Rotate_Clockwise_90();}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {T txy=4*X.x*X.y,x2=sqr(X.x),y2=sqr(X.y),r2=sqr(r),q=x2+y2-r2;return scale*q*MATRIX<T,TV::m>(txy,-4*x2-q,4*y2+q,-txy);}
    virtual T p(const TV& X,T t) const {return 0;}
    virtual TV F(const TV& X,T t) const {return (T)8*nu*scale*(2*sqr(r)-3*X.Magnitude_Squared())*X.Rotate_Clockwise_90();}
};
}

#endif
