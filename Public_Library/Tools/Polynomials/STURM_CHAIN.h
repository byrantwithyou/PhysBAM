//#####################################################################
// Copyright 2018, Frank Madrid, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STURM_CHAIN
//##################################################################### 
#ifndef __STURM_CHAIN__
#define __STURM_CHAIN__
#include <Tools/Polynomials/FIXED_POLYNOMIAL.h>

namespace PhysBAM{

template<class T,int d>
class STURM_CHAIN_HELPER
{
public:
    FIXED_POLYNOMIAL<T,d> p;
    STURM_CHAIN_HELPER<T,d-1> s;

    STURM_CHAIN_HELPER(const FIXED_POLYNOMIAL<T,d>& a,const FIXED_POLYNOMIAL<T,d-1>& b)
        :p(a),s(b,-(a%b))
    {
    }

    int Sign_Changes(T x,bool positive) const
    {
        int b=p(x)>=0;
        return (positive^b)+Sign_Changes(x,b);
    }
};

template<class T>
class STURM_CHAIN_HELPER<T,1>
{
public:
    FIXED_POLYNOMIAL<T,1> p;
    FIXED_POLYNOMIAL<T,0> r;

    STURM_CHAIN_HELPER(const FIXED_POLYNOMIAL<T,1>& a,const FIXED_POLYNOMIAL<T,0>& b)
        :p(a),r(b)
    {
    }   

    int Sign_Changes(T x,bool positive) const
    {
        int b=p(x)>=0;
        int c=r(x)>=0;
        return (positive^b)+(b^c);
    }
};

template<class T,int d>
class STURM_CHAIN
{
public:
    STURM_CHAIN_HELPER<T,d> s;

    STURM_CHAIN(const FIXED_POLYNOMIAL<T,d>& p)
        :s(p,p.Derivative())
    {
    }

    int Sign_Changes(T x) const
    {
        return s.s.Sign_Changes(x,s.p(x)>=0);
    }
};

}
#endif
