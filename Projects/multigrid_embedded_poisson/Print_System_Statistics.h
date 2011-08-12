//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PRINT_SYSTEM_STATISTICS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PRINT_SYSTEM_STATISTICS_HPP

#include <cmath>

#include <algorithm>
#include <iosfwd>
#include <limits>

#include <Jeffrey_Utilities/BASIC_TIMER.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_SYSTEM >
void Print_System_Statistics(
    const T_SYSTEM& system,
    const int n_index,
    std::ostream& lout)
{
    typedef typename T_SYSTEM::SCALAR_TYPE SCALAR_TYPE;

    BASIC_TIMER timer;

    lout << "Computing system statistics...";
    lout.flush();
    timer.Restart();
    SCALAR_TYPE max_diag = 0;
    SCALAR_TYPE min_diag = std::numeric_limits<SCALAR_TYPE>::infinity();
    SCALAR_TYPE max_abs_stencil_sum = 0;
    for(int index = 1; index <= n_index; ++index) {
        const SCALAR_TYPE diag = system.Diag(index);
        const SCALAR_TYPE stencil_sum = system.Stencil_Sum(index);
        if(diag != 0) {
            max_diag = std::max(max_diag, diag);
            min_diag = std::min(min_diag, diag);
        }
        max_abs_stencil_sum = std::max(max_abs_stencil_sum, std::abs(stencil_sum));
    }
    lout << timer.Elapsed() << " s" << std::endl;
    lout << "  max diag = " << max_diag << '\n'
         << "  min diag = " << min_diag << '\n'
         << "    ratio = " << max_diag / min_diag << '\n'
         << "  max abs stencil sum = " << max_abs_stencil_sum
         << std::endl;
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PRINT_SYSTEM_STATISTICS_HPP
