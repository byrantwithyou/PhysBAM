//#####################################################################
// Copyright 2006-2008, Gina Ma, Mauricio Flores
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MINRES
//#####################################################################
#ifndef __MINRES__
#define __MINRES__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
namespace PhysBAM{

template<class T>
class MINRES:public KRYLOV_SOLVER<T>
{
    typedef KRYLOV_SOLVER<T> BASE;
public:
    using BASE::restart_iterations;using BASE::residual_magnitude_squared;using BASE::iterations_used;using BASE::print_diagnostics;using BASE::print_residuals;
    using BASE::nullspace_measure;using BASE::nullspace_tolerance;using BASE::Solve;using BASE::Ensure_Size;

    MINRES()
    {}

    virtual ~MINRES();

//#####################################################################
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
        ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,const T tolerance,const int min_iterations,const int max_iterations);
    void Print_Diagnostics(int iterations);
//#####################################################################
};
}
#endif
