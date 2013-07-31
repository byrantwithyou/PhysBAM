//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/exchange.h>
#include <Tools/Math_Tools/min.h>
#include <Tools/Nonlinear_Equations/LINE_SEARCH.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <limits>
#include <cmath>
using std::abs;
using namespace PhysBAM;
//#####################################################################
// Function Update_Interval
//#####################################################################
template<class T> void LINE_SEARCH<T>::
Update_Interval(BRACKET& s,T d,T Fd)
{
    if(s.m>d){exchange(s.m,d);exchange(s.Fm,Fd);} // a < m < d < b
    T am=min(s.Fa,s.Fm),bd=min(s.Fb,Fd);
    if(am<bd){s.b=d;s.Fb=Fd;} // a m d
    else{s.a=s.m;s.m=d;s.Fa=s.Fm;s.Fm=Fd;} // m d b
}
//#####################################################################
// Function Compute_Quadratic_Minimum
//#####################################################################
template<class T> int LINE_SEARCH<T>::
Compute_Quadratic_Minimum(const BRACKET& s,T& x,T tolerance) // 0=terminate, 1=fail, 2=good
{
    T m=(s.Fa-s.Fm)*(s.m-s.b),n=(s.Fm-s.Fb)*(s.a-s.m),p=2*(n-m);
    if(p<0) return 1;
    if(p<=tolerance*(abs(s.Fa)+abs(s.Fm)+abs(s.Fb))*(abs(s.a)+abs(s.b))) return 0;
    T r=(s.a+s.m)*n-(s.m+s.b)*m;
    x=r/p;
    if(x<=s.a || x>=s.m || x==s.m) return 1;
    return 2;
}
//#####################################################################
// Function Best_Value
//#####################################################################
template<class T> T LINE_SEARCH<T>::
Best_Value(const BRACKET& s)
{
    if(s.Fa<=s.Fm) return (s.Fa<=s.Fb)?s.a:s.b;
    return (s.Fm<=s.Fb)?s.m:s.b;
}
//#####################################################################
// Function Line_Search_Quadratic_Golden_Section
//#####################################################################
template<class T> bool LINE_SEARCH<T>::
Line_Search_Quadratic_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance,T quadratic_tolerance)
{
    T m=(T).5*(a+b),t;
    BRACKET s={a,m,b,F(a),F(m),F(b)};
    int i;
    for(i=0;i<max_iterations;i++){
        if(abs(s.b-s.a)<interval_tolerance*(abs(s.a)+abs(s.b))) break;
        int r=Compute_Quadratic_Minimum(s,t,quadratic_tolerance);
        if(!r) break;
        if(r==1){
            if(s.m-s.a>s.b-s.m) t=(T).5*(s.a+s.m);
            else t=(T).5*(s.m+s.b);}
        Update_Interval(s,t,F(t));}
    x=Best_Value(s);
    return i<=max_iterations;
}
//#####################################################################
// Function Line_Search_Quadratic_Golden_Section
//#####################################################################
template<class T> bool LINE_SEARCH<T>::
Line_Search_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance)
{
    PHYSBAM_ASSERT(a<=b);
    T tau=(T).5*(sqrt((T)5)-1),m=a+tau*(b-a),t;
    BRACKET s={a,m,b,F(a),F(m),F(b)};
    int i;
    interval_tolerance*=abs(s.a)+abs(s.b);
    for(i=0;i<max_iterations;i++){
        if(abs(s.b-s.a)<interval_tolerance) break;
        t=s.a+s.b-s.m;
        Update_Interval(s,t,F(t));}
    x=Best_Value(s);
    return i<=max_iterations;
}
//#####################################################################
// Function Line_Search_Wolfe_Conditions
//#####################################################################
template<class T> bool LINE_SEARCH<T>::
Line_Search_Wolfe_Conditions(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,T c1,T c2,T x_max)
{
    PHYSBAM_ASSERT(0<c1 && c1<c2 && c2<1);
    WOLFE_HELPER a0={0};
    F.Compute(0,0,&a0.dfa,&a0.fa);
    PHYSBAM_ASSERT(a0.dfa<0);
    WOLFE_HELPER x0=a0,x1={1};
    for(;x1.a<=x_max;x1.a*=2){
        F.Compute(x1.a,0,&x1.dfa,&x1.fa);
        if(x1.fa>a0.fa+c1*x1.a*a0.dfa || x1.fa>=x0.fa)
            return Line_Search_Wolfe_Conditions_Zoom(F,x0,x1,a0,x,c1,c2);
        if(abs(x1.dfa)<=-c2*a0.dfa){
            x=x1.a;
            return true;}
        if(x1.dfa>=0)
            return Line_Search_Wolfe_Conditions_Zoom(F,x1,x0,a0,x,c1,c2);
        x0=x1;}
    x=0;
    return false;
}
//#####################################################################
// Function Line_Search_Wolfe_Conditions_Zoom
//#####################################################################
template<class T> bool LINE_SEARCH<T>::
Line_Search_Wolfe_Conditions_Zoom(NONLINEAR_FUNCTION<T(T)>& F,WOLFE_HELPER lo,WOLFE_HELPER hi,WOLFE_HELPER x0,T& x,T c1,T c2)
{
    T min_interval=abs(hi.a-lo.a)*std::numeric_limits<T>::epsilon();
    while(abs(hi.a-lo.a)>min_interval)
    {
        WOLFE_HELPER xj={(lo.a+hi.a)/2};
        F.Compute(xj.a,0,&xj.dfa,&xj.fa);
        if(xj.fa>x0.fa+c1*xj.a*x0.dfa || xj.fa>=x0.fa) hi=xj;
        else if(abs(xj.dfa)<=-c2*x0.dfa){
            x=xj.a;
            return true;}
        else{
            if(((hi.a>lo.a) && xj.dfa>=0) || ((hi.a<lo.a) && xj.dfa<=0)) hi=lo;
            lo=xj;}}
    if(lo.a>0) x=lo.a;
    x=1;
    return true;
}
namespace PhysBAM{
template struct LINE_SEARCH<float>;
template struct LINE_SEARCH<double>;
}
