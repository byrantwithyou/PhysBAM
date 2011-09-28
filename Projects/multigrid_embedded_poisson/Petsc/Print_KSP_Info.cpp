//#####################################################################
// Copyright 2010, Jeffrey Hellrung, Calvin Wang
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <ostream>

#include <petsc.h>
#include <petscksp.h>

#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/Petsc/CALL_AND_CHKERRQ.h>

#include <Jeffrey_Utilities/Petsc/Print_KSP_Info.h>

namespace PhysBAM
{

namespace Petsc
{

PetscErrorCode
Print_KSP_Info(
    const KSP petsc_ksp,
    std::ostream& lout /*= PhysBAM::nout*/)
{
    KSPConvergedReason petsc_converged_reason;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPGetConvergedReason(petsc_ksp, &petsc_converged_reason) );
    lout << "  KSPGetConvergedReason = "
         << KSPConvergedReasons[petsc_converged_reason]
         << " (" << petsc_converged_reason << ')'
         << std::endl;

    PetscInt petsc_iterations;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPGetIterationNumber(petsc_ksp, &petsc_iterations) );
    lout << "  # of iterations = " << petsc_iterations << std::endl;

    PetscReal petsc_residual_norm;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPGetResidualNorm(petsc_ksp, &petsc_residual_norm) );
    lout << "  (approximate preconditioned) residual norm = " << petsc_residual_norm << std::endl;

    PetscReal petsc_emax, petsc_emin;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPComputeExtremeSingularValues(petsc_ksp, &petsc_emax, &petsc_emin) );
    lout << "  (estimated) extreme singular values = " << petsc_emax << ", " << petsc_emin << std::endl;
    lout << "  (estimated) condition # = " << petsc_emax / petsc_emin << std::endl;

    return 0;
}

} // namespace Petsc

} // namespace PhysBAM
