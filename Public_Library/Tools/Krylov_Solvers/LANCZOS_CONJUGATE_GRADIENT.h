//#####################################################################
// Copyright 2006-2008, Gina Ma, Mauricio Flores
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LANCZOS_CONJUGATE_GRADIENT
//#####################################################################
#ifndef __LANCZOS_CONJUGATE_GRADIENT__
#define __LANCZOS_CONJUGATE_GRADIENT__

#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
namespace PhysBAM{

template<class T>
class LANCZOS_CONJUGATE_GRADIENT:public KRYLOV_SOLVER<T>
{
    typedef KRYLOV_SOLVER<T> BASE;
public:
    using BASE::restart_iterations;using BASE::residual_magnitude_squared;using BASE::iterations_used;
    using BASE::print_diagnostics;using BASE::print_residuals;using BASE::relative_tolerance;
    using BASE::nullspace_measure;using BASE::nullspace_tolerance;using BASE::Solve;using BASE::Ensure_Size;

    LANCZOS_CONJUGATE_GRADIENT()
    {}

    virtual ~LANCZOS_CONJUGATE_GRADIENT();

//#####################################################################
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
        ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,T tolerance,const int min_iterations,const int max_iterations);
    void Print_Diagnostics(int iterations);
//#####################################################################
};
}
#endif
