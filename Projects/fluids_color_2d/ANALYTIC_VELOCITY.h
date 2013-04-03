//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_VELOCITY__
#define __ANALYTIC_VELOCITY__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <boost/function.hpp>

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

template<class TV>
struct ANALYTIC_VELOCITY_F:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;

    T rho,mu;
    bool use_advection;

    ANALYTIC_VELOCITY_F(T rho,T mu,bool use_advection): rho(rho),mu(mu),use_advection(use_advection) {}
    virtual TV Lu(const TV& X,T t) const=0;
    virtual TV ut(const TV& X,T t) const=0;
    virtual TV dp(const TV& X,T t) const=0;
    virtual TV F(const TV& X,T t) const {TV f=rho*ut(X,t)+dp(X,t)-mu*Lu(X,t);if(use_advection) f+=rho*du(X,t)*u(X,t);return f;}
    virtual void Test(const TV& X) const
    {
        ANALYTIC_VELOCITY<TV>::Test(X);
        RANDOM_NUMBERS<T> rand;
        TV dX;
        T e=1e-6,t=rand.Get_Uniform_Number(0,1),dt=rand.Get_Uniform_Number(-e,e);
        rand.Fill_Uniform(dX,-e,e);
        TV Lun=-(T)2*TV::m*u(X,t),Lue=Lu(X,t);
        for(int i=0;i<TV::m;i++){
            TV dX;
            dX(i)=e;
            Lun+=u(X-dX,t)+u(X+dX,t);}

        T errLu=(Lun/sqr(e)-Lue).Magnitude();
        LOG::cout<<"analytic L u diff test "<<errLu<<std::endl;

        TV ua=u(X,t),ub=u(X,t+dt),dua=ut(X,t),dub=ut(X,t+dt);
        T errut=((dua+dub)*dt/2-(ub-ua)).Magnitude()/e;
        LOG::cout<<"analytic u_t diff test "<<errut<<std::endl;

        T p0=p(X,t),p1=p((X+dX),t);
        TV dp0=dp(X,t),dp1=dp((X+dX),t);
        T errp=abs((dp0+dp1).Dot(dX)/2-(p1-p0))/e;
        LOG::cout<<"analytic p diff test "<<errp<<std::endl;
    }
};

template<class TV>
struct ANALYTIC_VELOCITY_ELLIPSE_FLOW:public ANALYTIC_VELOCITY_F<TV>
{
    typedef typename TV::SCALAR T;
    T st,mu_j;
    bool sub_pj;
    boost::function<T(T t)> a,da,dda;

    ANALYTIC_VELOCITY_ELLIPSE_FLOW(T rho,T mu,bool use_advection,T st,T mu_j,bool sub_pj,boost::function<T(T t)> a,boost::function<T(T t)> da,boost::function<T(T t)> dda):
        ANALYTIC_VELOCITY_F<TV>(rho,mu,use_advection),st(st),mu_j(mu_j),sub_pj(sub_pj),a(a),da(da),dda(dda) {}
    virtual TV u(const TV& X,T t) const {return da(t)/a(t)*TV(X.x,-X.y);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return da(t)/a(t)*MATRIX<T,TV::m>(1,0,0,-1);}
    virtual TV Lu(const TV& X,T t) const {return TV();}
    virtual TV ut(const TV& X,T t) const {T b=a(t),c=da(t),d=dda(t),e=(d*b-c*c)/sqr(b);return e*TV(X.x,-X.y);}
    virtual T p(const TV& X,T t) const
    {
        if(!sub_pj) return 0;
        TV n=N(X,t);
        return -st*K(X,t)-2*mu_j*n.Dot(du(X,t)*n);
    }
    virtual TV dp(const TV& X,T t) const
    {
        if(!sub_pj) return TV();
        return -st*dK(X,t)-4*mu_j*dN(X,t).Transpose_Times(du(X,t).Symmetric_Part()*N(X,t));
    }
    T K(const TV& X,T t) const
    {
        T a4=sqr(sqr(a(t))),x2=sqr(X.x),y2=sqr(X.y),z=a4*y2,w=a4*z+x2,e=a4*pow(w,-(T)1.5);
        TV n=N(X,t);
        return e*(x2+z);
    }
    TV dK(const TV& X,T t) const
    {
        T a4=sqr(sqr(a(t))),x2=sqr(X.x),y2=sqr(X.y),z=a4*y2,w=a4*z+x2,e=a4*pow(w,-(T)2.5);
        return e*TV((2*w-3*(z+x2))*X.x,-(w+3*(a4-1)*x2)*a4*X.y);
    }
    TV N(const TV& X,T t) const {return TV(X.x,sqr(sqr(a(t)))*X.y).Normalized();}
    MATRIX<T,TV::m> dN(const TV& X,T t) const
    {
        T a4=sqr(sqr(a(t))),z=a4*X.y,x2=X.x*X.x,z2=z*z,xz=X.x*z,d2=x2+z2,d3=sqrt(d2)*d2;
        return MATRIX<T,TV::m>(z2,-xz,-a4*xz,x2*a4)/d3;
    }
    virtual void Test(const TV& X) const
    {
        ANALYTIC_VELOCITY_F<TV>::Test(X);
        RANDOM_NUMBERS<T> rand;
        TV dX;
        T e=1e-6,t=rand.Get_Uniform_Number(0,1);
        rand.Fill_Uniform(dX,-e,e);
        TV N0=N(X,t),N1=N((X+dX),t);
        MATRIX<T,TV::m> dN0=dN(X,t),dN1=dN((X+dX),t);
        T errN=((dN0+dN1)*dX/2-(N1-N0)).Magnitude()/e;
        LOG::cout<<"analytic normal diff test "<<errN<<std::endl;
        LOG::cout<<"analytic normal orthogonality test "<<dN0.Transpose_Times(N0)<<std::endl;
    }
};
}
#endif
