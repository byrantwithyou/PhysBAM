//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_KRYLOV_SYSTEM_COMPOSER_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_KRYLOV_SYSTEM_COMPOSER_HPP

#include <boost/cast.hpp>

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Jeffrey_Utilities/Krylov/Unwrap_Krylov_Vector.h>

namespace PhysBAM
{

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
struct KRYLOV_SYSTEM_COMPOSER
    : KRYLOV_SYSTEM_BASE<T>
{
    KRYLOV_SYSTEM_COMPOSER(
        const T_MULTIPLY& multiply,
        const T_INNER_PRODUCT& inner_product,
        const T_CONVERGENCE_NORM& convergence_norm,
        const T_PROJECT& project,
        const T_SET_BOUNDARY_CONDITIONS& set_boundary_conditions,
        const T_PROJECT_NULLSPACE& project_nullspace,
        const T_APPLY_PRECONDITIONER& apply_preconditioner,
        const bool use_preconditioner,
        const bool preconditioner_commutes_with_projection);

    void Multiply(
        const KRYLOV_VECTOR_BASE<T>& x,
        KRYLOV_VECTOR_BASE<T>& result) const PHYSBAM_OVERRIDE;

    double Inner_Product(
        const KRYLOV_VECTOR_BASE<T>& x,
        const KRYLOV_VECTOR_BASE<T>& y) const PHYSBAM_OVERRIDE;

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;

    void Project(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;

    // implements Sx=Px+x0, where x0 is the desired component of x in the nullspace of P
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;

    // removes component of x in nullspace of A (used to project residual for stopping conditions)
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;

private:
    void Apply_Preconditioner(
        const KRYLOV_VECTOR_BASE<T>& r,
        KRYLOV_VECTOR_BASE<T>& z) const PHYSBAM_OVERRIDE;

    const T_MULTIPLY m_multiply;
    const T_INNER_PRODUCT m_inner_product;
    const T_CONVERGENCE_NORM m_convergence_norm;
    const T_PROJECT m_project;
    const T_SET_BOUNDARY_CONDITIONS m_set_boundary_conditions;
    const T_PROJECT_NULLSPACE m_project_nullspace;
    const T_APPLY_PRECONDITIONER m_apply_preconditioner;
};

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>
Make_Krylov_System_Composer(
    const T_MULTIPLY& multiply,
    const T_INNER_PRODUCT& inner_product,
    const T_CONVERGENCE_NORM& convergence_norm,
    const T_PROJECT& project,
    const T_SET_BOUNDARY_CONDITIONS& set_boundary_conditions,
    const T_PROJECT_NULLSPACE& project_nullspace,
    const T_APPLY_PRECONDITIONER& apply_preconditioner,
    const bool use_preconditioner,
    const bool preconditioner_commutes_with_projection)
{
    return KRYLOV_SYSTEM_COMPOSER<
        T,
        T_VECTOR,
        T_MULTIPLY,
        T_INNER_PRODUCT,
        T_CONVERGENCE_NORM,
        T_PROJECT,
        T_SET_BOUNDARY_CONDITIONS,
        T_PROJECT_NULLSPACE,
        T_APPLY_PRECONDITIONER
    >(
        multiply,
        inner_product,
        convergence_norm,
        project,
        set_boundary_conditions,
        project_nullspace,
        apply_preconditioner,
        use_preconditioner,
        preconditioner_commutes_with_projection
    );
}

//#####################################################################
//#####################################################################

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline
KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>::
KRYLOV_SYSTEM_COMPOSER(
    const T_MULTIPLY& multiply,
    const T_INNER_PRODUCT& inner_product,
    const T_CONVERGENCE_NORM& convergence_norm,
    const T_PROJECT& project,
    const T_SET_BOUNDARY_CONDITIONS& set_boundary_conditions,
    const T_PROJECT_NULLSPACE& project_nullspace,
    const T_APPLY_PRECONDITIONER& apply_preconditioner,
    const bool use_preconditioner,
    const bool preconditioner_commutes_with_projection)
    : KRYLOV_SYSTEM_BASE<T>(
          use_preconditioner,
          preconditioner_commutes_with_projection
      ),
      m_multiply(multiply),
      m_inner_product(inner_product),
      m_convergence_norm(convergence_norm),
      m_project(project),
      m_set_boundary_conditions(set_boundary_conditions),
      m_project_nullspace(project_nullspace),
      m_apply_preconditioner(apply_preconditioner)
{ }

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline void
KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>::
Multiply(
    const KRYLOV_VECTOR_BASE<T>& x,
    KRYLOV_VECTOR_BASE<T>& result) const
{
    m_multiply(
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< const T_VECTOR* >(&x)),
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< T_VECTOR* >(&result))
    );
}

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline double
KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>::
Inner_Product(
    const KRYLOV_VECTOR_BASE<T>& x,
    const KRYLOV_VECTOR_BASE<T>& y) const
{
    return m_inner_product(
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< const T_VECTOR* >(&x)),
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< const T_VECTOR* >(&y))
    );
}

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline T
KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    return m_convergence_norm(
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< const T_VECTOR* >(&x))
    );
}

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline void
KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    m_project(
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< T_VECTOR* >(&x))
    );
}

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline void
KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
    m_set_boundary_conditions(
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< T_VECTOR* >(&x))
    );
}

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline void
KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    m_project_nullspace(
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< T_VECTOR* >(&x))
    );
}

template<
    class T,
    class T_VECTOR,
    class T_MULTIPLY,
    class T_INNER_PRODUCT,
    class T_CONVERGENCE_NORM,
    class T_PROJECT,
    class T_SET_BOUNDARY_CONDITIONS,
    class T_PROJECT_NULLSPACE,
    class T_APPLY_PRECONDITIONER
>
inline void
KRYLOV_SYSTEM_COMPOSER<
    T,
    T_VECTOR,
    T_MULTIPLY,
    T_INNER_PRODUCT,
    T_CONVERGENCE_NORM,
    T_PROJECT,
    T_SET_BOUNDARY_CONDITIONS,
    T_PROJECT_NULLSPACE,
    T_APPLY_PRECONDITIONER
>::
Apply_Preconditioner(
    const KRYLOV_VECTOR_BASE<T>& r,
    KRYLOV_VECTOR_BASE<T>& z) const
{
    m_apply_preconditioner(
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< const T_VECTOR* >(&r)),
        Unwrap_Krylov_Vector(*boost::polymorphic_downcast< T_VECTOR* >(&z))
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_KRYLOV_SYSTEM_COMPOSER_HPP
