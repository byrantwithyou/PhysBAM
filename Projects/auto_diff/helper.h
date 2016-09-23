#ifndef __HELPER__
#define __HELPER__

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Auto_Diff/DIFF_LAYOUT.h>
using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

#define TEST(...) Hess_Test([](auto u,auto v,auto a,auto b){return __VA_ARGS__;},#__VA_ARGS__)

template<class A> auto Diff_Helper(const A& a) -> typename enable_if<IS_SCALAR<decltype(a.x)>::value && !IS_VECTOR<A>::value,A>::type {return a;}
template<class A> auto Diff_Helper(const A& a) -> typename enable_if<IS_VECTOR<decltype(a.x)>::value,decltype(a.Dot(TV(1.1,2.3,-1.2)))>::type {return a.Dot(TV(1.1,2.3,-1.2));}

template<class A> auto Diff_Helper(const A& a) -> typename enable_if<IS_SCALAR<A>::value,A>::type {return a;}
template<class A> auto Diff_Helper(const A& a) -> typename enable_if<IS_VECTOR<A>::value,decltype(a.Dot(TV(1.1,2.3,-1.2)))>::type {return a.Dot(TV(1.1,2.3,-1.2));}

template<class F>
void Hess_Test(F f,const char* str)
{
    typedef DIFF_LAYOUT<T,TV::m,TV::m,-1,-1> LAYOUT;

    T eps=1e-6;
    TV du[2],u[2],v[2];
    T da[2],a[2],b[2];
    
    RANDOM_NUMBERS<T> rand;
    for(int i=0;i<2;i++){
        rand.Fill_Uniform(du[i],-eps,eps);
        rand.Fill_Uniform(da[i],-eps,eps);
        rand.Fill_Uniform(u[i],-1,1);
        rand.Fill_Uniform(a[i],-1,1);
        v[i]=u[i]+du[i];
        b[i]=a[i]+da[i];}


    auto v0=f(u[0],u[1],a[0],a[1]);
    auto d0=f(Diff_From_Var<LAYOUT,0>(u[0]),Diff_From_Var<LAYOUT,1>(u[1]),Diff_From_Var<LAYOUT,2>(a[0]),Diff_From_Var<LAYOUT,3>(a[1]));
    auto h0=f(Hess_From_Var<LAYOUT,0>(u[0]),Hess_From_Var<LAYOUT,1>(u[1]),Hess_From_Var<LAYOUT,2>(a[0]),Hess_From_Var<LAYOUT,3>(a[1]));
    auto h1=f(Hess_From_Var<LAYOUT,0>(v[0]),Hess_From_Var<LAYOUT,1>(v[1]),Hess_From_Var<LAYOUT,2>(b[0]),Hess_From_Var<LAYOUT,3>(b[1]));

    auto Sv0=Diff_Helper(v0);
    auto Sd0=Diff_Helper(d0);
    auto Sh0=Diff_Helper(h0);
    auto Sh1=Diff_Helper(h1);
    T x0=Sh0.x;
    T x1=Sh1.x;
    T y0=Sd0.x;
    T z0=Sv0;

    T dy=abs(x0-y0)/maxabs(x0,y0,(T)1e-30);
    T dz=abs(x0-z0)/maxabs(x0,z0,(T)1e-30);
    if(dy<1e-5 && dz<1e-5){}
    else printf("value %8.4f %8.4f %8.4f %10.4e %10.4e  %s\n",x0,y0,z0,dy,dz,str);

    VECTOR<VECTOR<T,3>,2> G0a,G0b;
    VECTOR<T,2> G2a,G2b;
    Extract<0>(G0a,Sh0.dx);
    Extract<0>(G0b,Sh1.dx);
    Extract<2>(G2a,Sh0.dx);
    Extract<2>(G2b,Sh1.dx);

    VECTOR<VECTOR<T,3>,2> D0a;
    VECTOR<T,2> D2a;
    Extract<0>(D0a,Sd0.dx);
    Extract<2>(D2a,Sd0.dx);

    T mG=sqrt(G0a.Flattened().Magnitude_Squared()+G2a.Magnitude_Squared());
    T mD=sqrt(D0a.Flattened().Magnitude_Squared()+D2a.Magnitude_Squared());
    T mS=sqrt((D0a-G0a).Flattened().Magnitude_Squared()+(D2a-G2a).Magnitude_Squared());
    T mm=abs(mS)/maxabs(mG,mD,(T)1e-30);
    if(mm<1e-5){}
    else printf("gradient %8.4f %8.4f %10.4e  %s\n",mG,mD,mm,str);
    
    MATRIX<MATRIX<T,3>,2,2> H00a,H00b;
    MATRIX<VECTOR<T,3>,2,2> H02a,H20a,H02b,H20b;
    MATRIX<T,2,2> H22a,H22b;
    Extract<0,0>(H00a,Sh0.ddx);
    Extract<0,0>(H00b,Sh1.ddx);
    Extract<2,0>(H20a,Sh0.ddx);
    Extract<2,0>(H20b,Sh1.ddx);
    Extract<0,2>(H02a,Sh0.ddx);
    Extract<0,2>(H02b,Sh1.ddx);
    Extract<2,2>(H22a,Sh0.ddx);
    Extract<2,2>(H22b,Sh1.ddx);

    T dx0=x1-x0,dx1=0;
    for(int i=0;i<2;i++){
        dx1+=(G0a(i)+G0b(i)).Dot(du[i])/2;
        dx1+=(G2a(i)+G2b(i))*da[i]/2;}

    T ddx0=0,ddx1=0,ddx2=0;
    for(int i=0;i<2;i++){
        ddx0+=(G0b(i)-G0a(i)).Magnitude_Squared();
        ddx0+=sqr(G2b(i)-G2a(i));
        TV dv;
        T ds=0;
        for(int j=0;j<2;j++){
            dv+=H00a(i,j)*du[j];
            dv+=H02a(i,j)*da[j];
            ds+=H20a(i,j).Dot(du[j]);
            ds+=H22a(i,j)*da[j];}
        ddx1+=dv.Magnitude_Squared();
        ddx1+=sqr(ds);
        ddx2+=(G0b(i)-G0a(i)-dv).Magnitude_Squared();
        ddx2+=sqr(G2b(i)-G2a(i)-ds);}
    ddx0=sqrt(ddx0);
    ddx1=sqrt(ddx1);
    ddx2=sqrt(ddx2);

    T r0=abs(dx0-dx1)/maxabs(dx0,dx1,(T)1e-30);
    T r1=ddx2/maxabs(ddx0,ddx1,(T)1e-30);
    if(r0<1e-5 && r1<1e-5){}
    else printf("%8.4f %8.4f %10.4e   %8.4f %8.4f %10.4e  %s\n",abs(dx0/eps),abs(dx1/eps),r0,ddx0/eps,ddx1/eps,r1,str);
}

#endif
