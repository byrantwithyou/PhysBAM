//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_POISSON_TEST
//#####################################################################
#ifndef __ANALYTIC_POISSON_TEST__
#define __ANALYTIC_POISSON_TEST__
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_POISSON_TEST
{
    typedef typename TV::SCALAR T;

    ANALYTIC_LEVELSET<TV>* analytic_levelset;
    ARRAY<ANALYTIC_SCALAR<TV>*> analytic_solution;
    ARRAY<T> mu;

    ANALYTIC_POISSON_TEST()
        :analytic_levelset(0)
    {}
    
    ~ANALYTIC_POISSON_TEST()
    {
        delete analytic_levelset;
        analytic_solution.Delete_Pointers_And_Clean_Memory();
    }

    T u_jump(const TV& X,int color0,int color1)
    {return analytic_solution(color1)->f(X,0)-analytic_solution(color0)->f(X,0);}

    T j_surface(const TV& X,int color0,int color1)
    {
        TV du0=analytic_solution(color0)->dX(X,0),du1=analytic_solution(color1)->dX(X,0);
        return (mu(color1)*du1-mu(color0)*du0).Dot(analytic_levelset->N(X,0,color1));
    }

    void Test(const RANGE<TV>& domain) const;
};
}
#endif
