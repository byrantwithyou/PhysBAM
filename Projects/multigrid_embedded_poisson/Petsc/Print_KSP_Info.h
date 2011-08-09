//#####################################################################
// Copyright 2010, Jeffrey Hellrung, Calvin Wang
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_PRINT_KSP_INFO_HPP
#define PHYSBAM_PUBLIC_LIBRARY_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_PRINT_KSP_INFO_HPP

#include <iostream>

#include <petsc.h>
#include <petscksp.h>

#include <Jeffrey_Utilities/Petsc/CALL_AND_CHKERRQ.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace Petsc
{

inline PetscErrorCode
Print_KSP_Info(const KSP petsc_ksp)
{
    KSPConvergedReason petsc_converged_reason;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPGetConvergedReason(petsc_ksp, &petsc_converged_reason) );
    std::cout << "  KSPGetConvergedReason = "
              << KSPConvergedReasons[petsc_converged_reason]
              << " (" << petsc_converged_reason << ')'
              << std::endl;

    PetscInt petsc_iterations;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPGetIterationNumber(petsc_ksp, &petsc_iterations) );
    std::cout << "  # of iterations = " << petsc_iterations << std::endl;

    PetscReal petsc_residual_norm;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPGetResidualNorm(petsc_ksp, &petsc_residual_norm) );
    std::cout << "  (approximate preconditioned) residual norm = " << petsc_residual_norm << std::endl;

    PetscReal petsc_emax, petsc_emin;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPComputeExtremeSingularValues(petsc_ksp, &petsc_emax, &petsc_emin) );
    std::cout << "  (estimated) extreme singular values = " << petsc_emax << ", " << petsc_emin << std::endl;
    std::cout << "  (estimated) condition # = " << petsc_emax / petsc_emin << std::endl;

    return 0;
}

} // namespace Petsc

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_PRINT_KSP_INFO_HPP
