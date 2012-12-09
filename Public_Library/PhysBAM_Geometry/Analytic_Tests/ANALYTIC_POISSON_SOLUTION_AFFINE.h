//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_POISSON_SOLUTION_AFFINE
//#####################################################################
#ifndef __ANALYTIC_POISSON_SOLUTION_AFFINE__
#define __ANALYTIC_POISSON_SOLUTION_AFFINE__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_POISSON_SOLUTION.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_POISSON_SOLUTION_AFFINE:public ANALYTIC_POISSON_SOLUTION<TV>
{
    typedef typename TV::SCALAR T;
    TV a;
    T b;

    ANALYTIC_POISSON_SOLUTION_AFFINE(TV a,T b): a(a),b(b) {}
    virtual ~ANALYTIC_POISSON_SOLUTION_AFFINE(){}

    virtual T u(const TV& X) const {return X.Dot(a)+b;}
    virtual TV du(const TV& X) const {return a;}
    virtual T Laplacian(const TV& X) const {return 0;}
};
}
#endif
