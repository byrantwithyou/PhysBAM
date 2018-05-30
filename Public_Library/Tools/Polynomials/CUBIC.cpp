//#####################################################################
// Copyright 2002, 2003, Zhaosheng Bao, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Math_Tools/exchange_sort.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Math_Tools/max.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Polynomials/CUBIC.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <cassert>
#include <cmath>
using namespace PhysBAM;

using ::std::abs;
using ::std::pow;
using ::std::sqrt;

//#####################################################################
// Constructor
//#####################################################################
template<class T> CUBIC<T>::
CUBIC(const T c3_input,const T c2_input,const T c1_input,const T c0_input)
    :c3(c3_input),c2(c2_input),c1(c1_input),c0(c0_input),roots(0),error_tolerance((T)1e-14)
{}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void CUBIC<T>::
Compute(const T x,T* ddf,T* df,T* f) const
{
    if(f) *f=((c3*x+c2)*x+c1)*x+c0;
    if(df) *df=(3*c3*x+2*c2)*x+c1;
    if(ddf) *ddf=6*c3*x+2*c2;
}
//#####################################################################
// Function Prime
//#####################################################################
template<class T> T CUBIC<T>::
Prime(const T x) const
{
    return (3*c3*x+2*c2)*x+c1;
}
//#####################################################################
// Function Compute_Roots_Noniterative_In_Interval
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Roots_Noniterative_In_Interval(const T xmin,const T xmax)
{
    Compute_Roots_Noniterative();
    if(roots == 1){if(root[0] < xmin || root[0] > xmax) roots=0;}
    else{ // roots=3
        if(root[2] < xmin || root[0] > xmax){roots=0;return;}
        if(root[1] < xmin){if(root[2] > xmax){roots=0;return;}else{roots=1;root[0]=root[2];return;}}
        if(root[1] > xmax){if(root[0] < xmin){roots=0;return;}else{roots=1;return;}}
        if(root[2] < xmax){if(root[0] > xmin) return;else{roots=2;root[0]=root[1];root[1]=root[2];}}
        if(root[0] > xmin){roots=2;return;}else{roots=1;root[0]=root[1];}}
}
//#####################################################################
// Function Compute_Roots_Noniterative
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Roots_Noniterative()
{
    if(c3 == 0){
        QUADRATIC<T> quadratic(c2,c1,c0);
        quadratic.Compute_Roots();
        roots=quadratic.roots;
        root[0]=quadratic.root[0];
        root[1]=quadratic.root[1];}
    else{
        T one_over_c3=1/c3,a1=c2*one_over_c3,a2=c1*one_over_c3,a3=c0*one_over_c3;
        T Q=((T)1/9)*a1*a1-((T)1/3)*a2,R=a1*((T)1/27*a1*a1-((T)1/6)*a2)+(T).5*a3;
        T Q_Q_Q=Q*Q*Q,R2_minus_Q3=R*R-Q_Q_Q;       
        if(R2_minus_Q3 <= 0){ // three real roots
            roots=3;
            T theta=acos(R/sqrt(Q_Q_Q)),theta_over_3=((T)1/3)*theta,minus_two_sqrt_Q=(T)(-2)*sqrt(Q),minus_a1_over_3=(-(T)1/3)*a1;
            root[0]=minus_two_sqrt_Q*cos(theta_over_3)+minus_a1_over_3;
            root[1]=minus_two_sqrt_Q*cos(theta_over_3+(T)pi*2/3)+minus_a1_over_3;
            root[2]=minus_two_sqrt_Q*cos(theta_over_3+(T)pi*4/3)+minus_a1_over_3;
            exchange_sort(root[0],root[1],root[2]);}      
        else{ // one real root
            roots=1;
            root[0]=pow(sqrt(R2_minus_Q3)+abs(R),((T)1/3));
            root[0]+=(T)Q/root[0];
            root[0]*=(R<(T)0)?(T)1:(T)-1;
            root[0]-=((T)1/3)*a1;}}
}
//#####################################################################
// Function Compute_Roots
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Roots()
{
    root[0]=0;root[1]=0;root[2]=0; // initialize
    if(c3 == 0){QUADRATIC<T> quadratic(c2,c1,c0);quadratic.Compute_Roots();roots=quadratic.roots;root[0]=quadratic.root[0];root[1]=quadratic.root[1];return;}
    else{ // c3 != 0 - cubic - bound the sign changes, i.e. roots
        T bound=(T)1.01*max((T)3*abs(c2/c3),sqrt((T)3*abs(c1/c3)),pow((T)3*abs(c0/c3),(T)1/3));
        Compute_Roots_In_Interval(-bound,bound);}
}
//#####################################################################
// Function Compute_Roots_In_Interval
//#####################################################################
template<class T> void CUBIC<T>::
Insert_Root_In_Extrema_Interval(const T xmin,const T xmax)
{
    T y0=Value(xmin),y1=Value(xmax);
    if(!y0) root[roots++]=xmin;
    else if(!y1) root[roots++]=xmax;
    else if((y0>0)!=(y1>0)){
        ITERATIVE_SOLVER<T> iterative_solver;
        iterative_solver.tolerance=error_tolerance;
        root[roots++]=iterative_solver.Bisection_Secant_Root(*this,xmin,xmax);}
}
//#####################################################################
// Function Compute_Roots_In_Interval
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Roots_In_Interval(const T xmin,const T xmax)
{
    if(c3==0){
        QUADRATIC<T> quadratic(c2,c1,c0);
        quadratic.Compute_Roots_In_Interval(xmin,xmax);
        roots=quadratic.roots;
        root[0]=quadratic.root[0];
        root[1]=quadratic.root[1];
        return;}

    INTERVAL<T> ivals[3];
    int intervals=0;
    Compute_Intervals(xmin,xmax,intervals,ivals);
    roots=0;
    for(int i=0;i<intervals;i++)
        Insert_Root_In_Extrema_Interval(ivals[i].min_corner,ivals[i].max_corner);
}
//#####################################################################
// Function Compute_Intervals
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Intervals(const T& a,const T& b,int& intervals,INTERVAL<T> interval[3])
{
    T ba=b-a,ba2=ba*ba,ba3=ba2*ba;
    T a2=a*a,a3=a2*a;
    CUBIC<T> new_cubic(c3*ba3,(3*c3*a+c2)*ba2,(3*c3*a2+2*c2*a+c1)*ba,c3*a3+c2*a2+c1*a+c0);
    new_cubic.Compute_Intervals(intervals,interval);
    for(int i=0;i<intervals;i++)
        interval[i]=interval[i]*(b-a)+a;
}
//#####################################################################
// Function Compute_Intervals
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Intervals(int& intervals,INTERVAL<T> interval[3])
{
    // Assume [a,b] are [0,1] and scale appropriately before calling this function
    T g2=c2,g1=c1,g0=c0;
    T h2=3*c3+c2;
    T h1=3*c3+2*c2+c1;
    T h0=c3+c2+c1+c0;
    
    int index=((g0>=0)!=(h0>=0));
    index|=(((g1>=0)!=(h1>=0))<<1);
    index|=(((g2>=0)!=(h2>=0))<<2);
    
    switch(index){
        case 0:
            intervals=0;
            break;
        case 1:
            intervals=1;
            interval[0]=INTERVAL<T>(0,1);
            break;
        case 2:
        case 6:
            {QUADRATIC<T> quadratic(3*c3,2*c2,c1);
            quadratic.Compute_Roots_In_Interval(0,1); // quadratic roots are the extrema of the cubic
            if(quadratic.roots!=1){LOG::cout << "2/6 more than one root" << std::endl;return;}
            if((Value(quadratic.root[0])>=0)==(g0>=0)) intervals=0;
            else{
                intervals=2;
                interval[0]=INTERVAL<T>(0,quadratic.root[0]);
                interval[1]=INTERVAL<T>(quadratic.root[0],1);}
            break;}
        case 3:
        case 7:
            {QUADRATIC<T> quadratic(3*c3,2*c2,c1);
            quadratic.Compute_Roots_In_Interval(0,1); // quadratic roots are the extrema of the cubic
            if(quadratic.roots!=1){LOG::cout << "3/7 more than one root" << std::endl;return;}
            intervals=1;
            if((Value(quadratic.root[0])>=0)==(g0>=0)) interval[0]=INTERVAL<T>(quadratic.root[0],1);
            else interval[0]=INTERVAL<T>(0,quadratic.root[0]);
            break;}
        case 4:
            {T inflection_point=-c2/(3*c3);
            QUADRATIC<T> quadratic(3*c3,2*c2,c1);
            if((quadratic.Value(inflection_point)>=0)==(g1>=0)) intervals=0;
            else{
                quadratic.Compute_Roots_In_Interval(0,1); // quadratic roots are the extrema of the cubic
                if(quadratic.roots!=2){LOG::cout << "4 should be two roots" << std::endl;return;}
                else if((Value(quadratic.root[1])>=0)!=(h0>=0)){
                    intervals=2;
                    interval[0]=INTERVAL<T>(quadratic.root[0],quadratic.root[1]);
                    interval[1]=INTERVAL<T>(quadratic.root[1],1);}
                else if((Value(quadratic.root[0])>=0)!=(g0>=0)){
                    intervals=2;
                    interval[0]=INTERVAL<T>(0,quadratic.root[0]);
                    interval[1]=INTERVAL<T>(quadratic.root[0],quadratic.root[1]);}
                else intervals=0;}
            break;}
        case 5:
            {T inflection_point=-c2/(3*c3);
            QUADRATIC<T> quadratic(3*c3,2*c2,c1);
            if((quadratic.Value(inflection_point)>=0)==(g1>=0)){
                intervals=1;
                interval[0]=INTERVAL<T>(0,1);}
            else{
                quadratic.Compute_Roots_In_Interval(0,1); // quadratic roots are the extrema of the cubic
                T e=Value(quadratic.root[0]);
                T f=Value(quadratic.root[1]);
                if(quadratic.roots!=2){
                    intervals=1;
                    interval[0]=INTERVAL<T>(0,1);}
                else if((e>=0)!=(g0>=0)){
                    if((e>=0)!=(f>=0)){
                        assert((f>=0)!=(h0>=0));
                        intervals=3;
                        interval[0]=INTERVAL<T>(0,quadratic.root[0]);
                        interval[1]=INTERVAL<T>(quadratic.root[0],quadratic.root[1]);
                        interval[2]=INTERVAL<T>(quadratic.root[1],1);}
                    else{
                        intervals=1;
                        interval[0]=INTERVAL<T>(0,quadratic.root[0]);}}
                else if((f>=0)!=(h0>=0)){
                    intervals=1;
                    interval[0]=INTERVAL<T>(quadratic.root[1],1);}
                else{
                    assert((e>=0)!=(f>=0));
                    intervals=1;
                    interval[0]=INTERVAL<T>(quadratic.root[0],quadratic.root[1]);}}
            break;}}
}
//#####################################################################
// Function Coefficients_From_Interpolation
//#####################################################################
template<class T> void CUBIC<T>::
Coefficients_From_Interpolation(T x0,T y0,T x1,T y1,T x2,T y2,T x3,T y3)
{
    T dx1=x0-x3,dx2=x1-x3,dx3=x2-x3,dy0=y0-y3,dy1=y1-y3,dy2=y2-y3;
    T dx12=dx1-dx2,dx13=dx1-dx3,dx23=dx2-dx3,den=1/(dx1*dx2*dx3*dx12*dx13*dx23);
    T sx1=dx1*dx1,sx2=dx2*dx2,sx3=dx3*dx3,n1=dx2*dx3,n2=dx1*dx3,n3=dx1*dx2;
    T k1=n1*dy0,k2=n2*dy1,k3=n3*dy2,m1=dx23*k1,m2=-dx13*k2,m3=dx12*k3;
    T p3=m1+m2+m3;
    T p2=(sx3-sx2)*k1+(sx1-sx3)*k2+(sx2-sx1)*k3;
    T p1=n1*m1+n2*m2+n3*m3;
    c3=den*p3;
    c2=den*(p2-3*p3*x3);
    c1=den*(p1+x3*(3*x3*p3-2*p2));
    c0=y3+den*x3*(x3*(p2-x3*p3)-p1);
}
//#####################################################################
namespace PhysBAM{
template class CUBIC<float>;
template class CUBIC<double>;
}
