//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Polynomials/QUADRATIC.h>
#include <cmath>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> QUADRATIC<T>::
QUADRATIC()
    :a(1),b(0),c(0),roots(1),root1(0),root2(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> QUADRATIC<T>::
QUADRATIC(const T a_input,const T b_input,const T c_input)
    :a(a_input),b(b_input),c(c_input),roots(0),root1(0),root2(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> QUADRATIC<T>::
~QUADRATIC()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void QUADRATIC<T>::
Compute(const T x,T* ddf,T* df,T* f) const
{
    if(f) *f=(a*x+b)*x+c;
    if(df) *df=2*a*x+b;
    if(ddf) *ddf=2*a;
}
//#####################################################################
// Function Coefficients_From_Interpolation
//#####################################################################
template<class T> void QUADRATIC<T>::
Coefficients_From_Interpolation(T x0,T y0,T x1,T y1,T x2,T y2)
{
    T dx1=x0-x2,dx2=x1-x2,dy0=y0-y2,dy1=y1-y2,den=1/(dx1*dx2*(dx1-dx2));
    T m1=dy0*dx2,m2=dx1*dy1,p2=m1-m2,p1=dx1*m2-dx2*m1;
    a=den*p2;
    b=den*(p1-2*p2*x2);
    c=den*x2*(x2*p2-p1)+y2;
}
//#####################################################################
// Function Compute_Roots
//#####################################################################
template<class T> void QUADRATIC<T>::
Compute_Roots()
{
    if(a==0){
        if(b==0){
            if(c==0){roots=-1;return;} // function is identically zero - a=b=c=0 - always a root!
            else{roots=0;return;}} // when a=b=0 and c != 0, there are no roots
        else{roots=1;root1=-c/b;return;}} // when a=0 and b != 0, there is one root
    else{ // a != 0
        T d=Discriminant();
        if(d<0){roots=0;return;} // no real roots
        else if(d==0){roots=1;root1=-b/(2*a);return;} // one root
        else{ // d > 0 - two real roots
            T radical;if(b>0) radical=-b-sqrt(d);else radical=-b+sqrt(d);
            roots=2;root1=radical/(2*a);root2=2*c/radical;if(root1>root2) exchange(root1,root2);return;}}
}
//#####################################################################
// Function Compute_Roots_In_Interval
//#####################################################################
template<class T> void QUADRATIC<T>::
Compute_Roots_In_Interval(const T xmin,const T xmax)
{
    Compute_Roots();
    if(roots==1){if(root1<xmin || root1>xmax) roots=0;}
    else if(roots==2){
        if(root2<xmin || root2>xmax) roots--;
        if(root1<xmin || root1>xmax){
            root1=root2;
            roots--;}}
}
namespace PhysBAM{
template class QUADRATIC<float>;
template class QUADRATIC<double>;
}
