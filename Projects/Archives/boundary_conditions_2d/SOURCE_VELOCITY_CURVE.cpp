#include "HEADER.h"
using namespace PhysBAM;

template<class T>
T Source_Velocity_Curve(bool d0,bool d1,T x)
{
    if(d1 && d0) return x*(1-x)*(2*x+5);
    if(!d1 && d0) return x*((16*x-39)*x+30)/6;
    if(d1 && !d0) return (1-x)*((32*x+5)*x+5)/6;
    return ((6-4*x)*x*x+5)/6;
}

template float Source_Velocity_Curve(bool d0,bool d1,float x);
template double Source_Velocity_Curve(bool d0,bool d1,double x);
