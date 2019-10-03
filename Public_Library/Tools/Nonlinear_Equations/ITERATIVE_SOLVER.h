//#####################################################################
// Copyright 2003-2008, Ronald Fedkiw, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ITERATIVE_SOLVER
//##################################################################### 
#ifndef __ITERATIVE_SOLVER__
#define __ITERATIVE_SOLVER__
#include <cassert>
#include <climits>
#include <cmath>
#include <functional>
#include <limits>

namespace PhysBAM{

using std::abs;

//#####################################################################
// Function Bisection_Root
//#####################################################################
template<class T,class FUNC>
T Bisection_Root(FUNC F,T a,T b,
    T tolerance=45*std::numeric_limits<T>::epsilon(),int max_iterations=INT_MAX)
{
    int iterations=0;
    T Fa=F(a);
    assert(Fa*F(b)<=0);
    while(b-a>tolerance && iterations++<max_iterations){
        T m=(T).5*(a+b),Fm=F(m);
        if(Fa*Fm<=0){b=m;}
        else{a=m;Fa=Fm;}}
    return (T).5*(a+b);
}
//#####################################################################
// Function Secant_Root
//#####################################################################
template<class T,class FUNC>
T Secant_Root(FUNC F,T x0,T x1,
    T tolerance=45*std::numeric_limits<T>::epsilon(),int max_iterations=INT_MAX)
{
    int iterations=0;T x_old=x0,x=x1,Fx_old=F(x_old),Fx=F(x);
    while(abs(Fx)>tolerance && iterations++<max_iterations){
        assert(abs(Fx-Fx_old)>tolerance);
        T x_temp=x;x-=Fx*(x-x_old)/(Fx-Fx_old);x_old=x_temp;Fx_old=Fx;Fx=F(x);}
    return x;
}
//#####################################################################
// Function Bisection_Secant_Root
//#####################################################################
template<class T,class FUNC>
T Bisection_Secant_Root(FUNC F,T a,T b,
    T tolerance=45*std::numeric_limits<T>::epsilon(),int max_iterations=INT_MAX)
{
    int iterations=0;
    T Fa=F(a),Fb=F(b),x_old=a,x=b,Fx_old=Fa,Fx=Fb;
    assert(Fa*Fb<=0);
    while(b-a>tolerance && iterations++<max_iterations){
        if(abs(Fx-Fx_old)<tolerance){ // bisection method
            T m=(T).5*(a+b),Fm=F(m);if(Fa*Fm<=0){b=m;Fb=Fm;}else{a=m;Fa=Fm;}
            x_old=a;x=b;Fx_old=Fa;Fx=Fb;} // update secant points
        else{ // secant method
            T x_temp=x;x-=Fx*(x-x_old)/(Fx-Fx_old);
            if(a<x && x<b){x_old=x_temp;Fx_old=Fx;Fx=F(x);T m=x,Fm=Fx;if(Fa*Fm<=0){b=m;Fb=Fm;}else{a=m;Fa=Fm;}} // update bisection points
            else{ // throw out secant root - do bisection instead
                T m=(T).5*(a+b),Fm=F(m);if(Fa*Fm<=0){b=m;Fb=Fm;}else{a=m;Fa=Fm;}
                x_old=a;x=b;Fx_old=Fa;Fx=Fb;}}} // update secant points
    return x;
}
//#####################################################################
// Function Golden_Minimum
//#####################################################################
template<class T,class FUNC>
T Golden_Minimum(FUNC F,T a,T b,
    T tolerance=45*std::numeric_limits<T>::epsilon(),int max_iterations=INT_MAX)
{
    T tau=(T).5*(sqrt((T)5)-1);
    int iterations=0;T x1=a+(1-tau)*(b-a),x2=a+tau*(b-a),Fx1=F(x1),Fx2=F(x2);
    while(b-a>tolerance && iterations++<max_iterations){
        if(Fx1>Fx2){a=x1;x1=x2;x2=a+tau*(b-a);Fx1=Fx2;Fx2=F(x2);}
        else{b=x2;x2=x1;x1=a+(1-tau)*(b-a);Fx2=Fx1;Fx1=F(x1);}}
    return (T).5*(a+b);
}
//#####################################################################
// Function Parabolic_Minimum
//#####################################################################
template<class T,class FUNC>
T Parabolic_Minimum(FUNC F,T x0,T x1,T x2,
    T tolerance=45*std::numeric_limits<T>::epsilon(),int max_iterations=INT_MAX)
{
    int iterations=0;T d0=x2-x1,d1=x0-x2,d2=x1-x0,Fx0=F(x0),Fx1=F(x1),Fx2=F(x2);
    while(abs(d0)>tolerance && abs(d1)>tolerance && abs(d2)>tolerance && iterations++<max_iterations){
        T d0F=d0*Fx0,d1F=d1*Fx1,d2F=d2*Fx2,denominator=d0F+d1F+d2F;
        if(abs(denominator)<tolerance) return x2;
        T x=(d0F*(x2+x1)+d1F*(x0+x2)+d2F*(x1+x0))/(2*denominator);
        x0=x1;x1=x2;x2=x;d0=x2-x1;d1=x0-x2;d2=x1-x0;Fx0=Fx1;Fx1=Fx2;Fx2=F(x2);}
    return x2;
}
//#####################################################################
// Function Golden_Parabolic_Minimum
//#####################################################################
template<class T,class FUNC>
T Golden_Parabolic_Minimum(FUNC F,T a,T b,
    T tolerance=45*std::numeric_limits<T>::epsilon(),int max_iterations=INT_MAX)
{
    int iterations=0;T m=(T).5*(a+b),Fa=F(a),Fb=F(b),Fm=F(m),x0=a,x1=b,x2=m,Fx0=Fa,Fx1=Fb,Fx2=Fm;
    while(b-a>tolerance && iterations++<max_iterations){ // parabolic interpolation
        T d0=(x2-x1)*Fx0,d1=(x0-x2)*Fx1,d2=(x1-x0)*Fx2,denominator=d0+d1+d2;
        if(abs(denominator)<tolerance) return x2;
        T x=(d0*(x2+x1)+d1*(x0+x2)+d2*(x1+x0))/(2*denominator);
        if(a<x && x<m){
            x0=x1;x1=x2;x2=x;Fx0=Fx1;Fx1=Fx2;Fx2=F(x2); // update parabolic points
            if(Fx2>Fm){a=x2;Fa=Fx2;}else{b=m;m=x2;Fb=Fm;Fm=Fx2;}} // update golden points
        else if(m<x && x<b){
            x0=x1;x1=x2;x2=x;Fx0=Fx1;Fx1=Fx2;Fx2=F(x2); // update parabolic points
            if(Fx2>Fm){b=x2;Fb=Fx2;}else{a=m;m=x2;Fa=Fm;Fm=Fx2;}} // update golden points
        else{ // golden section
            if(b-m<m-a){
                T x=(T).5*(a+m),Fx=F(x);
                if(Fx>Fm){a=x;Fa=Fx;}else{b=m;m=x;Fb=Fm;Fm=Fx;} // update golden points
                x0=a;x1=b;x2=m;Fx0=Fa;Fx1=Fb;Fx2=Fm;} // update parabolic points
            else{
                T x=(T).5*(m+b),Fx=F(x);
                if(Fx>Fm){b=x;Fb=Fx;}else{a=m;m=x;Fa=Fm;Fm=Fx;} // update golden points
                x0=a;x1=b;x2=m;Fx0=Fa;Fx1=Fb;Fx2=Fm;}}} // update parabolic points
    return x2;
}



}
#endif
