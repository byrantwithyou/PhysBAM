//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMQMR
//#####################################################################
#ifndef __SYMMQMR__
#define __SYMMQMR__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
namespace PhysBAM{

// Freund, R. and Nachtigal, N., "A New Krylov-Subspace Method For Symmetric Indefinite Linear Systems"
template<class T>
class SYMMQMR:public KRYLOV_SOLVER<T>
{
    typedef KRYLOV_SOLVER<T> BASE;
public:
    using BASE::restart_iterations;using BASE::residual_magnitude_squared;using BASE::iterations_used;using BASE::print_diagnostics;using BASE::print_residuals;
    using BASE::nullspace_measure;using BASE::nullspace_tolerance;using BASE::Solve;using BASE::Ensure_Size;

    SYMMQMR()
    {}

    virtual ~SYMMQMR();

//#####################################################################
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
        ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,const T tolerance,const int min_iterations,const int max_iterations);
//#####################################################################
};
}
#endif
