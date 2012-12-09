//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_POISSON_TEST
//#####################################################################
#ifndef __ANALYTIC_POISSON_TEST__
#define __ANALYTIC_POISSON_TEST__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_POISSON_SOLUTION.h>
#include <PhysBAM_Geometry/Finite_Elements/BOUNDARY_CONDITIONS_SCALAR_COLOR.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_POISSON_TEST:public BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>
{
    typedef typename TV::SCALAR T;

    ANALYTIC_LEVELSET<TV>* analytic_levelset;
    ARRAY<ANALYTIC_POISSON_SOLUTION<TV>*> analytic_solution;
    ARRAY<T> mu;

    ANALYTIC_POISSON_TEST()
        :analytic_levelset(0)
    {}
    
    virtual ~ANALYTIC_POISSON_TEST()
    {
        delete analytic_levelset;
        analytic_solution.Delete_Pointers_And_Clean_Memory();
    }

    T u_jump(const TV& X,int color0,int color1)
    {return analytic_solution(color1)->u(X)-analytic_solution(color0)->u(X);}

    T j_surface(const TV& X,int color0,int color1)
    {
        TV du0=analytic_solution(color0)->du(X),du1=analytic_solution(color1)->du(X);
        return (mu(color0)*du0-mu(color1)*du1).Dot(analytic_levelset->N(X,0,color1));
    }

    T n_surface(const TV& X,int color0,int color1)
    {
        TV du=analytic_solution(color1)->du(X);
        TV n=analytic_levelset->N(X,0,color1);
        return (mu(color1)*du).Dot(n);
    }

    T d_surface(const TV& X,int color0,int color1)
    {return analytic_solution(color1)->u(X);}

    void Test(const RANGE<TV>& domain) const;
};
}
#endif
