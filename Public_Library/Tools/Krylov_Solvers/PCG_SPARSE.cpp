//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Frederic Gibou, Geoffrey Irving, Michael Lentine, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <Tools/Krylov_Solvers/PCG_SPARSE_SYSTEM.h>
#include <Tools/Krylov_Solvers/SYMMQMR.h>
using namespace PhysBAM;
//#####################################################################
// Function Solve  
//#####################################################################
template<class T> void PCG_SPARSE<T>::
Solve(SPARSE_MATRIX_FLAT_MXN<T>& A_matrix,ARRAY<T>& x,ARRAY<T>& b,ARRAY<KRYLOV_VECTOR_BASE<T>*>& vectors,
    const T tolerance,const bool recompute_preconditioner)
{
    int desired_iterations=A_matrix.n-enforce_compatibility;if(maximum_iterations) desired_iterations=maximum_iterations;

    CONJUGATE_GRADIENT<T> cg;
    CONJUGATE_RESIDUAL<T> cr;
    SYMMQMR<T> symmqmr;
    KRYLOV_SOLVER<T>* solver=0;
    if(evolution_solver_type==krylov_solver_cg){solver=&cg;}
    else if(evolution_solver_type==krylov_solver_cr){solver=&cr;}
    else if(evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;}
    else PHYSBAM_FATAL_ERROR("Invalid krylov solver");
    solver->print_diagnostics=show_results;solver->print_residuals=show_residual;
    solver->restart_iterations=cg_restart_iterations;
    PCG_SPARSE_SYSTEM<T> system(*this,A_matrix);

    if(incomplete_cholesky){
        system.use_preconditioner=true;
        system.preconditioner_commutes_with_projection=false;
        if(recompute_preconditioner || !A_matrix.C)
            A_matrix.Construct_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
                preconditioner_zero_tolerance,preconditioner_zero_replacement);}

    KRYLOV_VECTOR_WRAPPER<T,ARRAY<T>&> kx(x),kb(b);
    PHYSBAM_ASSERT(&kx.v==&x && &kb.v==&b);
    solver->Solve(system,kx,kb,vectors,tolerance,0,desired_iterations);
}
//#####################################################################
namespace PhysBAM{
template class PCG_SPARSE<float>;
template class PCG_SPARSE<double>;
}
