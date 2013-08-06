//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOCAL_LINE_SEARCH
//##################################################################### 
#ifndef __LOCAL_LINE_SEARCH__
#define __LOCAL_LINE_SEARCH__
namespace PhysBAM{
template<class F> class NONLINEAR_FUNCTION;

template<class T>
struct LOCAL_LINE_SEARCH
{
    struct BRACKET
    {
        T a,m,b;
        T Fa,Fm,Fb;
    };

    static void Update_Interval(BRACKET& s,T d,T Fd); // a < b < c, update with new point d, where a < d < c.
    static int Compute_Quadratic_Minimum(const BRACKET& s,T& x,T tolerance); // 0=terminate, 1=fail, 2=good
    static bool Line_Search_Quadratic_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance,T quadratic_tolerance);
    static bool Line_Search_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance);
    static bool Line_Search_For_Wolfe_Conditions(NONLINEAR_FUNCTION<T(T)>& F,T gx,T a,T c1,T c2,T& x,int max_iterations,T a_tolerance);
    static T Zoom_Interval(NONLINEAR_FUNCTION<T(T)>& F,T gx,T gx_ai,T F_0,T a_lo,T a_hi,T c1,T c2,int max_iterations);
    static T Interpolate_Quadratic_Alpha(T Fa1,T Fa2,T gx_a1,T a1,T a2);
    static T Interpolate_Cubic_Alpha(T Fa1,T Fa2,T gx_a1,T gx_a2,T a1,T a2); 
    static T Best_Value(const BRACKET& s);

};
}
#endif
