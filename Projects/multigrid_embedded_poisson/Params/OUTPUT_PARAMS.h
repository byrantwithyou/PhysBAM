//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_OUTPUT_PARAMS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_OUTPUT_PARAMS_HPP

#include <string>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

struct OUTPUT_PARAMS
{
    std::string embedded_surface_filename;
    std::string cell_local_embedded_surface_filename_format;
    std::string infty_norm_error_filename;
};

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_OUTPUT_PARAMS_HPP
