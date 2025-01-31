//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ROBUST_FUNCTIONS
//#####################################################################
#ifndef __ROBUST_FUNCTIONS__
#define __ROBUST_FUNCTIONS__

#include <Core/Math_Tools/sqr.h>
#include <cmath>
namespace PhysBAM{

using ::std::abs;
using ::std::atan2;
using ::std::sin;

template<class T>
inline T sinc(const T x) // sin(x)/x
{if(abs(x)<1e-8) return 1;return sin(x)/x;}

template<class T>
inline T one_minus_cos_x_over_x_squared(const T x) // (1-cos(x))/x^2
{return (T).5*sqr(sinc((T).5*x));}

template<class T>
inline T one_minus_cos_x_over_x(const T x) // (1-cos(x))/x
{return one_minus_cos_x_over_x_squared(x)*x;}

template<class T>
inline T atan2_y_x_over_y(const T y,const T x) // atan2(y,x)/y
{if(abs(y)<1e-8) return 1;return atan2(y,x)/y;}

template<class T>
inline T log1p_x_over_x(const T x) // log(x+1)/x
{if(fabs(x)<1e-18) return 1;return log1p(x)/x;}

template<class T>
inline T diff_log_over_diff(const T x,const T y) // (log(x)-log(y))/(x-y)
{return log1p_x_over_x((y-x)/x)/x;}

template<class T>
inline T exp_x_minus_1_over_x(const T x) // (exp(x)-1)/x
{if(!x) return 1;return expm1(x)/x;}

template<class T>
inline T diff_exp_over_diff(const T x,const T y) // (exp(x)-exp(y))/(x-y)
{return exp(x)*exp_x_minus_1_over_x(y-x);}
}
#endif
