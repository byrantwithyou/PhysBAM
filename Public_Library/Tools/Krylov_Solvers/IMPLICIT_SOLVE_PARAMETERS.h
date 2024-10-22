//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_SOLVE_PARAMETERS
//#####################################################################
#ifndef __IMPLICIT_SOLVE_PARAMETERS__
#define __IMPLICIT_SOLVE_PARAMETERS__

#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
namespace PhysBAM{

template<class TV>
class IMPLICIT_SOLVE_PARAMETERS
{
    typedef typename TV::SCALAR T;
public:
    int cg_iterations;
    int cg_restart_iterations;
    T cg_tolerance;
    bool spectral_analysis;
    int lanczos_iterations;
    int cg_projection_iterations;
    bool throw_exception_on_backward_euler_failure;
    KRYLOV_SOLVER_TYPE evolution_solver_type;
    int project_nullspace_frequency;
    bool use_half_fully_implicit;
    bool test_system;
    bool print_matrix;

    IMPLICIT_SOLVE_PARAMETERS();
    IMPLICIT_SOLVE_PARAMETERS(const IMPLICIT_SOLVE_PARAMETERS&) = delete;
    void operator=(const IMPLICIT_SOLVE_PARAMETERS&) = delete;
    virtual ~IMPLICIT_SOLVE_PARAMETERS();
//#####################################################################
};
}
#endif
