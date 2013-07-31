//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINE_SEARCH
//##################################################################### 
#ifndef __LINE_SEARCH__
#define __LINE_SEARCH__
#include <cfloat>
namespace PhysBAM{
template<class F> class NONLINEAR_FUNCTION;

template<class T>
struct LINE_SEARCH
{
    struct BRACKET
    {
        T a,m,b;
        T Fa,Fm,Fb;
    };

    struct WOLFE_HELPER
    {
        T a,fa,dfa;
    };

    static void Update_Interval(BRACKET& s,T d,T Fd); // a < b < c, update with new point d, where a < d < c.
    static int Compute_Quadratic_Minimum(const BRACKET& s,T& x,T tolerance); // 0=terminate, 1=fail, 2=good
    static bool Line_Search_Quadratic_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance,T quadratic_tolerance);
    static bool Line_Search_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance);
    static bool Line_Search_Wolfe_Conditions(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,T c1,T c2,T x_max=FLT_MAX);
    static bool Line_Search_Wolfe_Conditions_Zoom(NONLINEAR_FUNCTION<T(T)>& F,WOLFE_HELPER lo,WOLFE_HELPER hi,WOLFE_HELPER x0,T& x,T c1,T c2);
    static T Best_Value(const BRACKET& s);
};
}
#endif
