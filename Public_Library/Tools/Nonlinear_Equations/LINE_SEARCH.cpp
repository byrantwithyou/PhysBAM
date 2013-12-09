//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/exchange.h>
#include <Tools/Math_Tools/maxabs.h>
#include <Tools/Math_Tools/min.h>
#include <Tools/Nonlinear_Equations/LINE_SEARCH.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <cmath>
#include <limits>
using std::abs;
using namespace PhysBAM;
//#####################################################################
// Function Dump_Line_Log
//#####################################################################
template<class T> void
Dump_Line_Log(NONLINEAR_FUNCTION<T(T)>& F,T a,T b)
{
    static int dump_id=0;
    T mn=1e-12,mx=1e4;
    int n=50;
    ARRAY<T> pts;
    ARRAY<T> dpts;
    T f,df,f0,df0;
    F.Compute(a,0,&df0,&f0);
    for(int i=0;i<=n;i++){
        F.Compute(a-pow((T)2,-(T)i)*(b-a),0,&df,&f);
        pts.Append(f);
        dpts.Append(df);}
    pts.Append(f0);
    dpts.Append(df0);
    for(int i=n;i>=0;i--){
        F.Compute(a+pow((T)2,-(T)i)*(b-a),0,&df,&f);
        pts.Append(f);
        dpts.Append(df);}

    FILE* red=fopen("red.txt","w");
    FILE* green=fopen("green.txt","w");
    for(int i=0;i<pts.m;i++){
        T x=pts(i)-f0,p=i-n-1;
        if(x>=0) fprintf(red,"%.16g %.16g\n",p,log10(clamp(x,mn,mx)));
        else fprintf(green,"%.16g %.16g\n",p,log10(clamp(-x,mn,mx)));}
    fclose(red);
    fclose(green);
    char buff[1000];
    sprintf(buff,"gnuplot -e \"set terminal postscript eps color ; set output 'func-%i.eps' ; plot 'red.txt' u 1:2 , 'green.txt' u 1:2\"",dump_id);
    if(system(buff)==0)
        LOG::printf("dump func-%i.eps\n",dump_id);
    else
        LOG::printf("dump func-%i.eps failed\n",dump_id);


    dump_id++;
}
//#####################################################################
// Function Dump_Line
//#####################################################################
template<class T> void
Dump_Line(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T c1)
{
    static int dump_id=0;
    int n=50,max_r=14;
    b-=a;
    T f,df,f0,df0;
    F.Compute(a,0,&df0,&f0);
    for(int r=0;r<max_r;r++){
        FILE* red=fopen("red.txt","w");
        FILE* green=fopen("green.txt","w");
        FILE* initial=fopen("initial.txt","w");
        T dx = 0.5*b/n;
        fprintf(initial,"%.16g %.16g %.16g %.16g\n",a,(T)0,dx,df0*dx);
        fclose(initial);
        for(int i=-n;i<=n;i++){
            T x=a+b*i/n;
            F.Compute(x,0,&df,&f);
            df *= dx;
            fprintf(f>=f0?red:green,"%.16g %.16g %.16g %.16g\n",x,f-f0,dx,df);}
        fclose(red);
        fclose(green);
        char buff[1000];
        sprintf(buff,"gnuplot -e \"set terminal postscript eps color ; set output 'func-%03i-%03i.eps' ; set st ar 1 nohead lc 1 ; set st ar 2 nohead lc 2 ; set st ar 3 nohead lc 3 ; plot 'red.txt' u 1:2:3:4 w vec as 1, 'green.txt' u 1:2:3:4 w vec as 2, 'initial.txt' u 1:2:3:4 w vec as 3 \"",dump_id,r);
        if(system(buff)==0)
            LOG::printf("dump func-%03i-%03i.eps\n",dump_id,r);
        else
            LOG::printf("dump func-%03i-%03i.eps failed\n",dump_id,r);

        b/=10;}
    dump_id++;
}
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
// Function Line_Search_Golden_Section
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
        if(s.a > a && abs(s.b-s.a)<interval_tolerance) break;
        t=s.a+s.b-s.m;
        Update_Interval(s,t,F(t));}
    x=Best_Value(s);
    return i<=max_iterations;
}
//#####################################################################
// Function Line_Search_Backtracking
//#####################################################################
template<class T> bool LINE_SEARCH<T>::
Line_Search_Backtracking(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,T c)
{
    PHYSBAM_ASSERT(0<c && c<1);
    WOLFE_HELPER x0={a};
    F.Compute(x0.a,0,&x0.dfa,&x0.fa);
    PHYSBAM_ASSERT(x0.dfa<0);
    WOLFE_HELPER x1={b};
    T min_interval=abs(a-b)*std::numeric_limits<T>::epsilon();
    while(abs(x1.a-x0.a)>min_interval){
        F.Compute(x1.a,0,&x1.dfa,&x1.fa);
        if(x1.fa<x0.fa+c*x1.a*x0.dfa){
            x=x1.a;
            return true;}
        x1.a=New_Point_Interpolation(x0,x1);
    }
    x=a;
    return false;
}
//#####################################################################
// Function Line_Search_Wolfe_Conditions
//#####################################################################
template<class T> bool LINE_SEARCH<T>::
Line_Search_Wolfe_Conditions(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,T c1,T c2,T x_max)
{
    PHYSBAM_ASSERT(0<c1 && c1<c2 && c2<1);
    WOLFE_HELPER a0={a};
    F.Compute(a0.a,0,&a0.dfa,&a0.fa);
    PHYSBAM_ASSERT(a0.dfa<0);
    WOLFE_HELPER x0=a0,x1={b};
    for(;x1.a<=x_max;x1.a=2*x1.a-x0.a){
        F.Compute(x1.a,0,&x1.dfa,&x1.fa);
        if(x1.fa>a0.fa+c1*x1.a*a0.dfa || x1.fa>=x0.fa){
            if(!Line_Search_Wolfe_Conditions_Zoom(F,x0,x1,a0,x,c1,c2))
                Line_Search_Derivative_Bisection(F,a,b,x,pow(std::numeric_limits<T>::epsilon(),2./3),x_max);
            return true;}
        if(abs(x1.dfa)<=-c2*a0.dfa){
            x=x1.a;
            return true;}
        if(x1.dfa>=0){
            if(!Line_Search_Wolfe_Conditions_Zoom(F,x1,x0,a0,x,c1,c2))
                Line_Search_Derivative_Bisection(F,a,b,x,pow(std::numeric_limits<T>::epsilon(),2./3),x_max);
            return true;}
        x0=x1;}
    x=a;
    return false;
}
//#####################################################################
// Function New_Point_Interpolation
//#####################################################################
template<class T> T LINE_SEARCH<T>::
New_Point_Interpolation(const WOLFE_HELPER& lo,const WOLFE_HELPER& hi)
{
    T k=hi.fa-lo.fa,a=-2*k+lo.dfa+hi.dfa,b=3*k-2*lo.dfa-hi.dfa,c=lo.dfa;
    T r=b*b-3*a*c;
    if(r>0){
        T s=sqrt(r),num=s-b,den=3*a;
        if(den<0){num=-num;den=-den;}
        if(num>(T).1*den && num<(T).9*den)
            return lo.a+num/den*(hi.a-lo.a);}
    return (lo.a+hi.a)/2;
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
        WOLFE_HELPER xj={New_Point_Interpolation(lo,hi)};
        F.Compute(xj.a,0,&xj.dfa,&xj.fa);
        if(xj.fa>x0.fa+c1*xj.a*x0.dfa || xj.fa>=x0.fa) hi=xj;
        else if(abs(xj.dfa)<=-c2*x0.dfa){
            x=xj.a;
            return true;}
        else{
            if(((hi.a>lo.a) && xj.dfa>=0) || ((hi.a<lo.a) && xj.dfa<=0)) hi=lo;
            lo=xj;}}
    if(lo.a-x0.a>=min_interval*100){
        LOG::printf("take decrease on zoom with %g\n",lo.a);
        x=lo.a;
        return true;}
    LOG::printf("exit zoom on %g\n",lo.a);
    return false;
}
//#####################################################################
// Function Line_Search_Derivative_Bisection
//#####################################################################
template<class T> void LINE_SEARCH<T>::
Line_Search_Derivative_Bisection(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,T allowed_relative_increase,T x_max)
{
    WOLFE_HELPER x0={a},x1={b};
    F.Compute(x0.a,0,&x0.dfa,&x0.fa);
    WOLFE_HELPER a0=x0;
    PHYSBAM_ASSERT(a0.dfa<0);
    F.Compute(x1.a,0,&x1.dfa,&x1.fa);
    if(x1.fa-a0.fa<=allowed_relative_increase*abs(a0.fa)){
        x=b;
        return;}

    Line_Search_Upper_Trust_Interval(F,x0,x1,a0,allowed_relative_increase);

    if(x1.dfa>0){
        T min_interval_width=maxabs(x1.a-a0.a,x1.a,a0.a)*std::numeric_limits<T>::epsilon()*2;
        while(x1.a-x0.a>min_interval_width){
            WOLFE_HELPER xj={(x1.a+x0.a)/2};
            F.Compute(xj.a,0,&xj.dfa,&xj.fa);
            if(xj.fa-a0.fa>allowed_relative_increase*abs(a0.fa)){
                x1=xj;
                Line_Search_Upper_Trust_Interval(F,x0,x1,a0,allowed_relative_increase);
                if(x1.dfa<=0) break;}
            else if(xj.dfa<=0) x0=xj;
            else x1=xj;}}
    x=x1.a;
}
//#####################################################################
// Function Line_Search_Upper_Trust_Interval
//#####################################################################
template<class T> void LINE_SEARCH<T>::
Line_Search_Upper_Trust_Interval(NONLINEAR_FUNCTION<T(T)>& F,WOLFE_HELPER x0,WOLFE_HELPER& x1,WOLFE_HELPER a0,T allowed_relative_increase)
{
    T min_interval_width=maxabs(x1.a-a0.a,x1.a,a0.a)*std::numeric_limits<T>::epsilon()*2;
    while(x1.a-x0.a>min_interval_width){
//        LOG::printf("%.16g %.16g  %.16g %.16g  %.16g %.16g\n",x0.a,x1.a,x0.fa,x1.fa,x0.dfa,x1.dfa);
        WOLFE_HELPER xj={(x1.a+x0.a)/2};
        F.Compute(xj.a,0,&xj.dfa,&xj.fa);
        if(xj.fa-a0.fa>allowed_relative_increase*abs(a0.fa)) x1=xj;
        else x0=xj;}
    x1=x0;
}
namespace PhysBAM{
template struct LINE_SEARCH<float>;
template struct LINE_SEARCH<double>;
}
