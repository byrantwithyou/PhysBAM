//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_SOLVER_PARAMS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_SOLVER_PARAMS_HPP

#include <algorithm>
#include <limits>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

struct SOLVER_PARAMS
{
    enum SOLVER_ID
    {
        SOLVER_ID_NULL = 0,
        SOLVER_ID_PHYSBAM_CG,
        SOLVER_ID_PHYSBAM_MINRES,
        SOLVER_ID_PETSC_CG,
        SOLVER_ID_PETSC_MINRES,
        SOLVER_ID_MG,
        SOLVER_ID_MGPCG
    } solver_id;

    unsigned int max_iterations;
    float relative_tolerance;
    float absolute_tolerance;
    bool print_residuals;
    bool precondition;

    struct MULTIGRID_PARAMS
    {
        unsigned int n_level;
        unsigned int mu;

        enum SMOOTHER_ID
        {
            SMOOTHER_ID_NULL = 0,
            SMOOTHER_ID_GAUSS_SEIDEL,
            SMOOTHER_ID_WEIGHTED_JACOBI
        } smoother_id;

        struct NU_PARAMS
        {
            unsigned int n;
            unsigned int n_pre_full_embedded;
            unsigned int n_full;
            unsigned int n_post_full_embedded;

            NU_PARAMS()
                : n(0),
                  n_pre_full_embedded(0),
                  n_full(0),
                  n_post_full_embedded(0)
            { }

            NU_PARAMS Construct_Symmetric() const
            {
                NU_PARAMS result = *this;
                std::swap(result.n_pre_full_embedded, result.n_post_full_embedded);
                return result;
            }

            bool Symmetric(const NU_PARAMS& other) const
            {
                return n == other.n &&
                       n_pre_full_embedded == other.n_post_full_embedded &&
                       n_full == other.n_full &&
                       n_post_full_embedded == other.n_pre_full_embedded;
            }

            friend bool Symmetric(const NU_PARAMS& nu1, const NU_PARAMS& nu2)
            { return nu1.Symmetric(nu2); }
        };

        NU_PARAMS nu_pre_restrict;
        NU_PARAMS nu_post_prolong;

        MULTIGRID_PARAMS()
            : n_level(0),
              mu(0),
              smoother_id(SMOOTHER_ID_NULL)
        { }
    } multigrid;

    SOLVER_PARAMS()
        : solver_id(SOLVER_ID_NULL),
          max_iterations(std::numeric_limits< unsigned int >::max()),
          relative_tolerance(std::numeric_limits< float >::epsilon()),
          absolute_tolerance(std::numeric_limits< float >::min()),
          print_residuals(false),
          precondition(true)
    { }
};

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_SOLVER_PARAMS_HPP
