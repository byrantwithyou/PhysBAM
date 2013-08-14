//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOCAL_TR_CONJUGATE_GRADIENT
//#####################################################################
#ifndef __LOCAL_TR_CONJUGATE_GRADIENT__
#define __LOCAL_TR_CONJUGATE_GRADIENT__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
//#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
namespace PhysBAM{

// see Golub and Van Loan, 10.2.6, p. 529 for details

// cost per iteration = 1 matrix multiply/project/precondition, 2 inner products, 1 convergence norm, 3 saxpy's
// approximate total flops = 2v + 10n

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
        ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,T tolerance,const int min_iterations,const int max_iterations,T trust_region)=0;
//#####################################################################
};

template<class T>
class LOCAL_TR_CONJUGATE_GRADIENT:public KRYLOV_SOLVER<T>
{
    typedef KRYLOV_SOLVER<T> BASE;
public:
    using BASE::restart_iterations;using BASE::residual_magnitude_squared;using BASE::iterations_used;
    using BASE::print_diagnostics;using BASE::print_residuals;using BASE::relative_tolerance;
    using BASE::nullspace_measure;using BASE::nullspace_tolerance;using BASE::Solve;using BASE::Ensure_Size;
    bool finish_before_indefiniteness;

    LOCAL_TR_CONJUGATE_GRADIENT()
        :finish_before_indefiniteness(false)
    {}

    virtual ~LOCAL_TR_CONJUGATE_GRADIENT();

//#####################################################################
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
        ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,T tolerance,const int min_iterations,const int max_iterations,T trust_region);
//#####################################################################
};
}
#endif
