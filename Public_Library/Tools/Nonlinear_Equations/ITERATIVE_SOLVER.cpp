//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ITERATIVE_SOLVER
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/sqr.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Nonlinear_Equations/LINE_SEARCH.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Nonlinear_Equations/PARAMETRIC_LINE.h>
#include <cfloat>
#include <cmath>
using namespace PhysBAM;
using std::abs;
//#####################################################################
// Function Bisection_Root
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Bisection_Root(NONLINEAR_FUNCTION<T(T)>& F,T a,T b)
{
    iterations=0;T Fa=F(a);assert(Fa*F(b)<=0);
    while(b-a>tolerance && iterations++<max_iterations){
        T m=(T).5*(a+b),Fm=F(m);if(Fa*Fm<=0){b=m;}else{a=m;Fa=Fm;}}

    if(iterations==max_iterations){
        LOG::cout<<"Bisection_Root failed with max iterations"<<std::endl;
        return 0;}
    return (T).5*(a+b);
}
//#####################################################################
// Function Newton_Root
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Newton_Root(NONLINEAR_FUNCTION<T(T)>& F,T x0)
{
    iterations=0;T x=x0,Fx=F(x),Fprime=F.Prime(x);
    while(abs(Fx)>tolerance && iterations++<max_iterations){
        assert(abs(Fprime)>tolerance);x-=Fx/Fprime;Fx=F(x);Fprime=F.Prime(x);}
    if(iterations==max_iterations){
        LOG::cout<<"Newton_Root failed with max iterations"<<std::endl;
        return 0;}
    return x;
}
//#####################################################################
// Function Secant_Root
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Secant_Root(NONLINEAR_FUNCTION<T(T)>& F,T x0,T x1)
{
    iterations=0;T x_old=x0,x=x1,Fx_old=F(x_old),Fx=F(x);
    while(abs(Fx)>tolerance && iterations++<max_iterations){
        assert(abs(Fx-Fx_old)>tolerance);T x_temp=x;x-=Fx*(x-x_old)/(Fx-Fx_old);x_old=x_temp;Fx_old=Fx;Fx=F(x);}
    if(iterations==max_iterations){
        LOG::cout<<"Secant_Root failed with max iterations"<<std::endl;
        return 0;}
    return x;
}
//#####################################################################
// Function Bisection_Secant_Root
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Bisection_Secant_Root(NONLINEAR_FUNCTION<T(T)>& F,T a,T b)
{
    iterations=0;T Fa=F(a),Fb=F(b),x_old=a,x=b,Fx_old=Fa,Fx=Fb;assert(Fa*Fb<=0);
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
    if(iterations==max_iterations){
        LOG::cout<<"Bisection_Secant_Root failed with max iterations"<<std::endl;
        return 0;}
    return x;
}
//#####################################################################
// Function Bisection_Secant_Root_For_Thin_Shells
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Bisection_Secant_Root_For_Thin_Shells(NONLINEAR_FUNCTION<T(T)>& F,T a,T b)
{
    iterations=0;T Fa=F(a),Fb=F(b),x_old=a,x=b,Fx_old=Fa,Fx=Fb;assert(Fa!=FLT_MAX&&Fb!=FLT_MAX&&Fa*Fb<=0);
    while(b-a>tolerance && iterations++<max_iterations){
        if(abs(Fx-Fx_old)<tolerance){ // bisection method
            T m=(T).5*(a+b),Fm=F(m);
            if(Fm==FLT_MAX){return m;}// thin shells hack
            else if(Fa*Fm<=0){b=m;Fb=Fm;}else{a=m;Fa=Fm;}
            x_old=a;x=b;Fx_old=Fa;Fx=Fb;} // update secant points
        else{ // secant method
            T x_temp=x;x-=Fx*(x-x_old)/(Fx-Fx_old);
            if(a<x && x<b){x_old=x_temp;Fx_old=Fx;Fx=F(x);T m=x,Fm=Fx;
                if(Fm==FLT_MAX){return m;} // thin shells hack
                else if(Fa*Fm<=0){b=m;Fb=Fm;}else{a=m;Fa=Fm;}} // update bisection points
            else{ // throw out secant root - do bisection instead
                T m=(T).5*(a+b),Fm=F(m);
                if(Fm==FLT_MAX){return m;} // thin shells hack
                else if(Fa*Fm<=0){b=m;Fb=Fm;}else{a=m;Fa=Fm;}
                x_old=a;x=b;Fx_old=Fa;Fx=Fb;}}} // update secant points
    if(iterations==max_iterations){
        LOG::cout<<"Bisection_Secant_Root failed with max iterations"<<std::endl;
        return 0;}
    return x;
}
//#####################################################################
// Function Bisection_Newton_Root
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Bisection_Newton_Root(NONLINEAR_FUNCTION<T(T)>& F,T a,T b)
{
    iterations=0;T Fa=F(a),x=a,Fx=Fa,Fprime=F.Prime(x);assert(Fa*F(b)<=0);
    while(b-a>tolerance && iterations++<max_iterations){
        if(abs(Fprime)<tolerance){ // bisection method
            T m=(T).5*(a+b),Fm=F(m);if(Fa*Fm<=0){b=m;}else{a=m;Fa=Fm;}
            x=a;Fx=Fa;Fprime=F.Prime(x);} // update newton points
        else{ // newton method
            x-=Fx/Fprime;
            if(a<x && x<b){Fx=F(x);Fprime=F.Prime(x);T m=x,Fm=Fx;if(Fa*Fm<=0){b=m;}else{a=m;Fa=Fm;}} // update bisection points
            else{ // throw out secant root - do bisection instead
                T m=(T).5*(a+b),Fm=F(m);if(Fa*Fm<=0){b=m;}else{a=m;Fa=Fm;}
                x=a;Fx=Fa;Fprime=F.Prime(x);}}} // update newton points
    if(iterations==max_iterations){
        LOG::cout<<"Newton_Secant_Root failed with max iterations"<<std::endl;
        return 0;}
    return x;
}
//#####################################################################
// Function Golden_Minimum
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Golden_Minimum(NONLINEAR_FUNCTION<T(T)>& F,T a,T b)
{
    T tau=(T).5*(sqrt((T)5)-1);iterations=0;T x1=a+(1-tau)*(b-a),x2=a+tau*(b-a),Fx1=F(x1),Fx2=F(x2);
    while(b-a>tolerance && iterations++<max_iterations){
        if(Fx1>Fx2){a=x1;x1=x2;x2=a+tau*(b-a);Fx1=Fx2;Fx2=F(x2);}
        else{b=x2;x2=x1;x1=a+(1-tau)*(b-a);Fx2=Fx1;Fx1=F(x1);}}
    if(iterations==max_iterations){
        LOG::cout<<"Golden_Minimum failed with max iterations"<<std::endl;
        return 0;}
    return (T).5*(a+b);
}
//#####################################################################
// Function Parabolic_Minimum
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Parabolic_Minimum(NONLINEAR_FUNCTION<T(T)>& F,T x0,T x1,T x2)
{
    iterations=0;T d0=x2-x1,d1=x0-x2,d2=x1-x0,Fx0=F(x0),Fx1=F(x1),Fx2=F(x2);
    while(abs(d0)>tolerance && abs(d1)>tolerance && abs(d2)>tolerance && iterations++<max_iterations){
        T d0F=d0*Fx0,d1F=d1*Fx1,d2F=d2*Fx2,denominator=d0F+d1F+d2F;
        if(abs(denominator)<tolerance) return x2;
        T x=(d0F*(x2+x1)+d1F*(x0+x2)+d2F*(x1+x0))/(2*denominator);
        x0=x1;x1=x2;x2=x;d0=x2-x1;d1=x0-x2;d2=x1-x0;Fx0=Fx1;Fx1=Fx2;Fx2=F(x2);}
    if(iterations==max_iterations){
        LOG::cout<<"Parabolic_Minimum failed with max iterations"<<std::endl;
        return 0;}
    return x2;
}
//#####################################################################
// Function Golden_Parabolic_Minimum
//#####################################################################
template<class T> T ITERATIVE_SOLVER<T>::
Golden_Parabolic_Minimum(NONLINEAR_FUNCTION<T(T)>& F,T a,T b)
{
    iterations=0;T m=(T).5*(a+b),Fa=F(a),Fb=F(b),Fm=F(m),x0=a,x1=b,x2=m,Fx0=Fa,Fx1=Fb,Fx2=Fm;
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
    if(iterations==max_iterations){
        LOG::cout<<"Golden_Parabolic_Minimum failed with max iterations"<<std::endl;
        return 0;}
    return x2;
}
namespace PhysBAM{
template class ITERATIVE_SOLVER<float>;
template class ITERATIVE_SOLVER<double>;
}
