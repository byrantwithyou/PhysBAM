//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_SOLVER_PARAMS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_SOLVER_PARAMS_HPP

#include <limits>

namespace PhysBAM
{

struct SOLVER_PARAMS
{
    bool precondition;
    unsigned int max_iterations;
    float relative_tolerance;
    float absolute_tolerance;
    bool print_diagnostics;
    bool print_residuals;

    SOLVER_PARAMS()
        : precondition(true),
          max_iterations(std::numeric_limits< unsigned int >::max()),
          relative_tolerance(std::numeric_limits< float >::epsilon()),
          absolute_tolerance(std::numeric_limits< float >::min()),
          print_diagnostics(false),
          print_residuals(false)
    { }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_SOLVER_PARAMS_HPP
