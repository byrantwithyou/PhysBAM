//#####################################################################
// Copyright 2008, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KRYLOV_SOLVER
//#####################################################################
#ifndef __KRYLOV_SOLVER__
#define __KRYLOV_SOLVER__
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{
template<class T_MATRIX> struct IS_MATRIX;
enum KRYLOV_SOLVER_TYPE {krylov_solver_cg,krylov_solver_cr,krylov_solver_symmqmr};

template<class T>
class KRYLOV_SOLVER
{
public:
    bool print_diagnostics,print_residuals,relative_tolerance;
    T nullspace_tolerance; // don't attempt to invert eigenvalues approximately less than nullspace_tolerance*max_eigenvalue
    int* iterations_used;
    T residual_magnitude_squared,nullspace_measure; // extra convergence information
    int restart_iterations;

    KRYLOV_SOLVER();
    virtual ~KRYLOV_SOLVER();

//#####################################################################
    static void Ensure_Size(ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,const KRYLOV_VECTOR_BASE<T>& v,int size);
    virtual bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
        ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,T tolerance,const int min_iterations,const int max_iterations)=0;
//#####################################################################
};
}
#endif
