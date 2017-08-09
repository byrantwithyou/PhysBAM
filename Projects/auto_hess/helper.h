#ifndef __HELPER__
#define __HELPER__

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

#define TEST(...) Hess_Test([](auto u,auto v,auto w,auto a,auto b,auto z){return __VA_ARGS__;},#__VA_ARGS__)

#define TEST2(...) Diff_Test([](auto u,auto v,auto w,auto a,auto b,auto z){return __VA_ARGS__;},#__VA_ARGS__)

template<int Q> AUTO_HESS<T,TV,Q> Diff_Helper(const AUTO_HESS<T,TV,Q>& a) {return a;}
template<int Q> AUTO_HESS<T,TV,Q> Diff_Helper(const AUTO_HESS<VECTOR<T,2>,TV,Q>& a) {return a.Dot(VECTOR<T,2>(1.1,2.3));}
template<int Q> AUTO_HESS<T,TV,Q> Diff_Helper(const AUTO_HESS<VECTOR<T,3>,TV,Q>& a) {return a.Dot(TV(1.1,2.3,-1.2));}
template<int Q,int m,int n> AUTO_HESS<T,TV,Q> Diff_Helper(const AUTO_HESS<MATRIX<T,m,n>,TV,Q>& a) {return a.Transpose_Times(a).Trace();}
template<int Q,int m> AUTO_HESS<T,TV,Q> Diff_Helper(const AUTO_HESS<SYMMETRIC_MATRIX<T,m>,TV,Q>& a) {return a.Transpose_Times(a).Trace();}
inline AUTO_HESS<T,TV,3> Diff_Helper(const VECTOR<T,2>& a) {return AUTO_HESS<T,TV,3>(a.Dot(VECTOR<T,2>(1.1,2.3)));}
inline AUTO_HESS<T,TV,3> Diff_Helper(const VECTOR<T,3>& a) {return AUTO_HESS<T,TV,3>(a.Dot(TV(1.1,2.3,-1.2)));}
inline AUTO_HESS<T,TV,3> Diff_Helper(T a) {return AUTO_HESS<T,TV,3>(a);}
template<int m,int n> AUTO_HESS<T,TV,3> Diff_Helper(const MATRIX<T,m,n>& a) {return AUTO_HESS<T,TV,3>(a.Transpose_Times(a).Trace());}
template<int m> AUTO_HESS<T,TV,3> Diff_Helper(const SYMMETRIC_MATRIX<T,m>& a) {return AUTO_HESS<T,TV,3>(a.Transpose_Times(a).Trace());}

template<class F>
void Hess_Test(F f,const char* str)
{
    T eps=1e-6;
    TV du,u,v;
    
    RANDOM_NUMBERS<T> rand;
    rand.Fill_Uniform(du,-eps,eps);
    rand.Fill_Uniform(u,-1,1);
    v=u+du;

    auto uv=u,vv=exp(uv),wv=sin(uv);
    auto av=vv.Magnitude_Squared(),bv=wv.Magnitude_Squared();
    VECTOR<T,2> zv(av,bv);
    auto v0=f(uv,vv,wv,av,bv,zv);

    auto un=AUTO_NO_DIFF<TV,TV>::From_Var(u),vn=exp(un),wn=sin(un);
    auto an=vn.Magnitude_Squared(),bn=wn.Magnitude_Squared();
    AUTO_NO_DIFF<VECTOR<T,2>,TV> zn(an,bn);
    auto n0=f(un,vn,wn,an,bn,zn);

    auto ud=AUTO_DIFF<TV,TV>::From_Var(u),vd=exp(ud),wd=sin(ud);
    auto ad=vd.Magnitude_Squared(),bd=wd.Magnitude_Squared();
    AUTO_DIFF<VECTOR<T,2>,TV> zd(ad,bd);
    auto d0=f(ud,vd,wd,ad,bd,zd);

    auto uh=AUTO_HESS<TV,TV>::From_Var(u),vh=exp(uh),wh=sin(uh);
    auto ah=vh.Magnitude_Squared(),bh=wh.Magnitude_Squared();
    AUTO_HESS<VECTOR<T,2>,TV> zh(ah,bh);
    auto h0=f(uh,vh,wh,ah,bh,zh);

    auto ui=AUTO_HESS<TV,TV>::From_Var(v),vi=exp(ui),wi=sin(ui);
    auto ai=vi.Magnitude_Squared(),bi=wi.Magnitude_Squared();
    AUTO_HESS<VECTOR<T,2>,TV> zi(ai,bi);
    auto h1=f(ui,vi,wi,ai,bi,zi);

    auto Sv0=Diff_Helper(v0);
    auto Sn0=Diff_Helper(n0);
    auto Sd0=Diff_Helper(d0);
    auto Sh0=Diff_Helper(h0);
    auto Sh1=Diff_Helper(h1);
    T x0=Sh0.x;
    T x1=Sh1.x;
    T y0=Sd0.x;
    T z0=Sv0.x;
    T w0=Sn0.x;

    T dx=abs(z0-x0)/maxabs(z0,x0,(T)1e-30);
    T dy=abs(z0-y0)/maxabs(z0,y0,(T)1e-30);
    T dw=abs(z0-w0)/maxabs(z0,w0,(T)1e-30);
    if(dw>1e-5) printf("value (NONE) %8.4f %8.4f %10.4e  %s\n",z0,w0,dw,str);
    if(dy>1e-5) printf("value (DIFF) %8.4f %8.4f %10.4e  %s\n",z0,y0,dy,str);
    if(dx>1e-5) printf("value (HESS) %8.4f %8.4f %10.4e  %s\n",z0,x0,dx,str);

    TV Ga=Sh0.dx,Gb=Sh1.dx;
    TV Da=Sd0.dx;

    T mG=Ga.Magnitude();
    T mD=Da.Magnitude();
    T mS=(Da-Ga).Magnitude();
    T mm=abs(mS)/maxabs(mG,mD,(T)1e-30);
    if(mm<1e-5){}
    else printf("gradient %8.4f %8.4f %10.4e  %s\n",mG,mD,mm,str);
    
    MATRIX<T,3> Ha=Sh0.ddx,Hb=Sh1.ddx;

    T dx0=x1-x0,dx1=(Ga+Gb).Dot(du)/2;
    TV dv0=Gb-Ga,dv1=(Ha+Hb)*du/2;
    T ddx0=dv0.Magnitude(),ddx1=dv1.Magnitude(),ddx2=(dv0-dv1).Magnitude();
    T r0=abs(dx0-dx1)/maxabs(dx0,dx1,(T)1e-30);
    T r1=ddx2/maxabs(ddx0,ddx1,(T)1e-30);
    if(r0<1e-4 && r1<1e-4){}
    else printf("%8.4f %8.4f %10.4e   %8.4f %8.4f %10.4e  %s\n",abs(dx0/eps),abs(dx1/eps),r0,ddx0/eps,ddx1/eps,r1,str);
}

template<class F>
void Diff_Test(F f,const char* str)
{
    T eps=1e-6;
    TV du,u,v;
    
    RANDOM_NUMBERS<T> rand;
    rand.Fill_Uniform(du,-eps,eps);
    rand.Fill_Uniform(u,-1,1);
    v=u+du;

    auto uv=u,vv=exp(uv),wv=sin(uv);
    auto av=vv.Magnitude_Squared(),bv=wv.Magnitude_Squared();
    VECTOR<T,2> zv(av,bv);
    auto v0=f(uv,vv,wv,av,bv,zv);

    auto un=AUTO_NO_DIFF<TV,TV>::From_Var(u),vn=exp(un),wn=sin(un);
    auto an=vn.Magnitude_Squared(),bn=wn.Magnitude_Squared();
    AUTO_NO_DIFF<VECTOR<T,2>,TV> zn(an,bn);
    auto n0=f(un,vn,wn,an,bn,zn);

    auto ud=AUTO_DIFF<TV,TV>::From_Var(u),vd=exp(ud),wd=sin(ud);
    auto ad=vd.Magnitude_Squared(),bd=wd.Magnitude_Squared();
    AUTO_DIFF<VECTOR<T,2>,TV> zd(ad,bd);
    auto d0=f(ud,vd,wd,ad,bd,zd);

    auto ui=AUTO_DIFF<TV,TV>::From_Var(v),vi=exp(ui),wi=sin(ui);
    auto ai=vi.Magnitude_Squared(),bi=wi.Magnitude_Squared();
    AUTO_DIFF<VECTOR<T,2>,TV> zi(ai,bi);
    auto d1=f(ui,vi,wi,ai,bi,zi);

    auto Sv0=Diff_Helper(v0);
    auto Sn0=Diff_Helper(n0);
    auto Sd0=Diff_Helper(d0);
    auto Sd1=Diff_Helper(d1);
    T y0=Sd0.x;
    T y1=Sd1.x;
    T z0=Sv0.x;
    T w0=Sn0.x;

    T dy=abs(z0-y0)/maxabs(z0,y0,(T)1e-30);
    T dw=abs(z0-w0)/maxabs(z0,w0,(T)1e-30);
    if(dw>1e-5) printf("value (NONE) %8.4f %8.4f %10.4e  %s\n",z0,w0,dw,str);
    if(dy>1e-5) printf("value (DIFF) %8.4f %8.4f %10.4e  %s\n",z0,y0,dy,str);

    TV Ga=Sd0.dx,Gb=Sd1.dx;
    T dx0=y1-y0,dx1=(Ga+Gb).Dot(du)/2;
    T r0=abs(dx0-dx1)/maxabs(dx0,dx1,(T)1e-30);
    if(r0>1e-4) printf("%8.4f %8.4f %10.4e  %s\n",abs(dx0/eps),abs(dx1/eps),r0,str);
}

#endif
