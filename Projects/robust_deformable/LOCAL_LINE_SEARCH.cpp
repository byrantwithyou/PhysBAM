//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/exchange.h>
#include <Tools/Math_Tools/min.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <cmath>
#include <fstream>
#include "LOCAL_LINE_SEARCH.h"
using std::abs;
using namespace PhysBAM;
//#####################################################################
// Function Update_Interval
//#####################################################################
template<class T> void LOCAL_LINE_SEARCH<T>::
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
template<class T> int LOCAL_LINE_SEARCH<T>::
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
template<class T> T LOCAL_LINE_SEARCH<T>::
Best_Value(const BRACKET& s)
{
    if(s.Fa<=s.Fm) return (s.Fa<=s.Fb)?s.a:s.b;
    return (s.Fm<=s.Fb)?s.m:s.b;
}
//#####################################################################
// Function Line_Search_Quadratic_Golden_Section
//#####################################################################
template<class T> bool LOCAL_LINE_SEARCH<T>::
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
template<class T> bool LOCAL_LINE_SEARCH<T>::
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
//####################################################################
// Function Line_Search_For_Wolfe_Conditions
//####################################################################
template<class T> bool LOCAL_LINE_SEARCH<T>::
Line_Search_For_Wolfe_Conditions(NONLINEAR_FUNCTION<T(T)>& F,T gx,T a,T c1,T c2,T& x,int max_iterations,T a_tolerance)
{    
    T a_i=a;
    T aBefore=0;
    int i;
    T F_i;
    T F_max=0;
    T F_0=F(0);
    T gx_ai;

    PHYSBAM_ASSERT(gx<0);
    for(i=0;i<max_iterations;i++){
        F_i=F(a_i);
        if(F_i-F_0>-c1*gx*a_i){
            gx_ai=F.Prime(aBefore);
            a_i=Zoom_Interval(F,gx,gx_ai,F_0,aBefore,a_i,c1,c2,max_iterations);
            break;}
        else{         
            gx_ai=F.Prime(a_i);
            if(fabs(gx_ai)<=-c2*gx) break;
            if(gx_ai>=0.0){
                a_i=Zoom_Interval(F,gx,gx_ai,F_0,a_i,aBefore,c1,c2,max_iterations);
                break;}
            if(i==0) F_max=F(a_tolerance);
            aBefore=a_i;
            a_i=Interpolate_Quadratic_Alpha(F_i,F_max,gx_ai,a_i,a_tolerance);}}

    x = a_i;

    PHYSBAM_ASSERT(x>=0); 
//    PHYSBAM_ASSERT(F(x)<=F_0);
 
    return i<=max_iterations;
}
//####################################################################
// Function Zoom_Interval
//####################################################################
template<class T> T LOCAL_LINE_SEARCH<T>::
Zoom_Interval(NONLINEAR_FUNCTION<T(T)>& F,T gx,T gx_ai,T F_0,T a_lo,T a_hi,T c1,T c2,int max_iterations)
{
    T a_j;
    int j;
    T F_temp;
    T F_lo=F(a_lo);
    T F_hi=F(a_hi);
    T gx_lo=gx_ai;
    T gx_hi=F.Prime(a_hi);
    T gx_aj;
    for(j=0;j<max_iterations;j++){
        a_j=Interpolate_Quadratic_Alpha(F_hi,F_lo,gx_hi,a_hi,a_lo);
        //a_j=Interpolate_Cubic_Alpha(F_lo,F_hi,gx_lo,gx_hi,a_lo,a_hi);
        F_temp=F(a_j);
        if(F_temp>=F_lo){
            a_hi=a_j;
            F_hi=F_temp;
            gx_hi=F.Prime(a_hi);}
        else{
            gx_aj=F.Prime(a_j);
            if(fabs(gx_aj)<=-c2*gx){
                return a_j;}
            if(gx_aj*(a_hi-a_lo)>=0){
                a_hi=a_lo;
                F_hi=F_lo;
                gx_hi=gx_lo;}
            a_lo=a_j;
            F_lo=F_temp;
            gx_lo=gx_aj;}
    }
    return 1e-3;
}
//####################################################################
// Function Interpolate_Quadratic_Alpha
//####################################################################
template<class T> T LOCAL_LINE_SEARCH<T>::
Interpolate_Quadratic_Alpha(T Fa1,T Fa2,T gx_a1,T a1,T a2) 
 {
     T d1=(T)2*a1*(Fa1-Fa2)+(a2*a2-a1*a1)*gx_a1;
     T d2=(T)2*(Fa1-Fa2+(a2-a1)*gx_a1);

     return d1/d2; 
 }

//####################################################################
// Function Interpolate_Cubic_Alpha
//####################################################################
template<class T> T LOCAL_LINE_SEARCH<T>::
Interpolate_Cubic_Alpha(T Fa1,T Fa2,T gx_a1,T gx_a2,T a1,T a2)
 {
     printf("gx_lo: %lf, gx_hi: %lf\n",Fa1,Fa2);
     T d1=gx_a1+gx_a2-(T)3*(Fa1-Fa2)/(a1-a2);
     if(d1*d1-gx_a1*gx_a2<0) return -1.0;
     T d2=sqrt(d1*d1-gx_a1*gx_a2);

     return a2-(a2-a1)*(gx_a2+d2-d1)/(gx_a2-gx_a1+(T)2*d2);
 }
namespace PhysBAM{
template struct LOCAL_LINE_SEARCH<float>;
template struct LOCAL_LINE_SEARCH<double>;
}
