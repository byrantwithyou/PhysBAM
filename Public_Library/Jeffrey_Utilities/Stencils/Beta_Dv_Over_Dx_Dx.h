//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_BETA_DV_OVER_DX_DX_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_BETA_DV_OVER_DX_DX_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D >
inline VECTOR<T,D>
Beta_Dv_Over_Dx_Dx(const T beta, VECTOR<T,D> dx)
{
    const T beta_dv = beta * dx.Product();
    for(int d = 1; d <= D; ++d)
        dx[d] = beta_dv / ((1<<(D-1)) * (dx[d] * dx[d]));
    return dx;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_BETA_DV_OVER_DX_DX_HPP
