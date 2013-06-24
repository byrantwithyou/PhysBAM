//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_POISSON_SOLUTION_QUADRATIC
//#####################################################################
#ifndef __ANALYTIC_POISSON_SOLUTION_QUADRATIC__
#define __ANALYTIC_POISSON_SOLUTION_QUADRATIC__
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Analytic_Tests/ANALYTIC_POISSON_SOLUTION.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_POISSON_SOLUTION_QUADRATIC:public ANALYTIC_POISSON_SOLUTION<TV>
{
    typedef typename TV::SCALAR T;
    SYMMETRIC_MATRIX<T,TV::m> M;
    TV a;
    T b;

    ANALYTIC_POISSON_SOLUTION_QUADRATIC(const MATRIX<T,TV::m>& M,TV a,T b): M(M.Symmetric_Part()),a(a),b(b) {}
    virtual ~ANALYTIC_POISSON_SOLUTION_QUADRATIC(){}

    virtual T u(const TV& X) const {return (M*X+a).Dot(X)+b;}
    virtual TV du(const TV& X) const {return M*X*2+a;}
    virtual T Laplacian(const TV& X) const {return M.Trace()*2;}
};
}
#endif
