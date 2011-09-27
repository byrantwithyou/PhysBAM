//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP

#include <iosfwd>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <Jeffrey_Utilities/GENERIC_SYSTEM_REFERENCE.h>
#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/SOLVER_PARAMS.h>

namespace PhysBAM
{

template< class T >
void
Solve_SPD_System_With_ICC_PCG(
    const unsigned int n_thread,
    const SOLVER_PARAMS params,
    const bool has_constant_vectors_in_null_space,
    const GENERIC_SYSTEM_REFERENCE<T> system,
    const ARRAY_VIEW<T> rhs,
    ARRAY_VIEW<T> x,
    std::ostream& lout = PhysBAM::nout);

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
