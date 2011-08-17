//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#include <cassert>

#include <algorithm>
#include <limits>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/NULL_FUNCTION.h>
#include <Jeffrey_Utilities/GENERIC_SYSTEM_REFERENCE.h>
#include <Jeffrey_Utilities/Krylov/DOT_PRODUCT_INNER_PRODUCT.h>
#include <Jeffrey_Utilities/Krylov/JACOBI_PRECONDITIONER.h>
#include <Jeffrey_Utilities/Krylov/KRYLOV_SYSTEM_COMPOSER.h>
#include <Jeffrey_Utilities/Krylov/KRYLOV_VECTOR_WRAPPER_ARRAY_VIEW.h>
#include <Jeffrey_Utilities/Krylov/MAXIMUM_MAGNITUDE_CONVERGENCE_NORM.h>
#include <Jeffrey_Utilities/SOLVER_PARAMS.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>

#include <Jeffrey_Utilities/Krylov/Solve_SPD_System_With_ICC_PCG.h>

namespace PhysBAM
{

template< class T >
void
Solve_SPD_System_With_ICC_PCG(
    const SOLVER_PARAMS params,
    const GENERIC_SYSTEM_REFERENCE<T> system,
    const ARRAY_VIEW<T> rhs,
    const bool has_constant_vectors_in_null_space,
    ARRAY_VIEW<T> x)
{
    // TODO: If params.precondition, then use ICC preconditioner.
    // TODO: If non-trivial null space, use project_nullspace.
    static_cast<void>(has_constant_vectors_in_null_space);
    const int n_index = rhs.Size();
    assert(n_index == rhs.Size());
    assert(n_index == x.Size());
    CONJUGATE_GRADIENT<T> cg;
    cg.print_diagnostics = params.print_diagnostics;
    cg.print_residuals = params.print_residuals;
    ARRAY<T> q(n_index); // init'ed to 0
    ARRAY<T> s(n_index); // init'ed to 0
    ARRAY<T> r(n_index); // init'ed to 0
    ARRAY<T> k(n_index); // init'ed to 0
    ARRAY<T> z(n_index); // init'ed to 0
    system.Apply(x, q);
    q -= rhs;
    const float initial_residual = static_cast< float >(ARRAYS_COMPUTATIONS::Maximum_Magnitude(q));
    const float tolerance = std::max(params.relative_tolerance * initial_residual, params.absolute_tolerance);
    const int max_iterations =
        params.max_iterations == std::numeric_limits< unsigned int >::max()
     && n_index <= std::numeric_limits< int >::max() / 4
      ? 4 * n_index
      : static_cast< int >(std::min(
            params.max_iterations,
            static_cast< unsigned int >(std::numeric_limits< int >::max())
        ));
    typedef KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> > KRYLOV_VECTOR_TYPE;
    cg.Solve(
        Make_Krylov_System_Composer< T, KRYLOV_VECTOR_TYPE >(
            PHYSBAM_BOUND_FAST_MEM_FN( system, &GENERIC_SYSTEM_REFERENCE<T>::Apply ),
            DOT_PRODUCT_INNER_PRODUCT(),
            MAXIMUM_MAGNITUDE_CONVERGENCE_NORM<T>(),
            NULL_FUNCTION(), // project
            NULL_FUNCTION(), // set_boundary_conditions
            NULL_FUNCTION(), // project_nullspace
            Make_Safe_Jacobi_Preconditioner(
                PHYSBAM_BOUND_FAST_MEM_FN( system, &GENERIC_SYSTEM_REFERENCE<T>::Diag )
            ),
            true, // use_preconditioner
            true  // preconditioner_commutes_with_projection
        ),
        const_cast< KRYLOV_VECTOR_TYPE& >(KRYLOV_VECTOR_TYPE(x)),
        KRYLOV_VECTOR_TYPE(rhs),
        const_cast< KRYLOV_VECTOR_TYPE& >(KRYLOV_VECTOR_TYPE(As_Array_View(q))),
        const_cast< KRYLOV_VECTOR_TYPE& >(KRYLOV_VECTOR_TYPE(As_Array_View(s))),
        const_cast< KRYLOV_VECTOR_TYPE& >(KRYLOV_VECTOR_TYPE(As_Array_View(r))),
        const_cast< KRYLOV_VECTOR_TYPE& >(KRYLOV_VECTOR_TYPE(As_Array_View(k))),
        const_cast< KRYLOV_VECTOR_TYPE& >(KRYLOV_VECTOR_TYPE(As_Array_View(z))),
        tolerance,
        0, // min_iterations
        max_iterations
    );
}

#define EXPLICIT_INSTANTIATION( T ) \
template void Solve_SPD_System_With_ICC_PCG<T>( \
    const SOLVER_PARAMS params, \
    const GENERIC_SYSTEM_REFERENCE<T> system, \
    const ARRAY_VIEW<T> rhs, \
    const bool has_constant_vectors_in_null_space, \
    ARRAY_VIEW<T> x);
EXPLICIT_INSTANTIATION( float )
EXPLICIT_INSTANTIATION( double )
#undef EXPLICIT_INSTANTIATION

} // namespace PhysBAM
