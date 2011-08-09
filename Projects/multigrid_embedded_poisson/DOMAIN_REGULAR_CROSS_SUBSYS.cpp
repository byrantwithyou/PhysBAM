//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "DOMAIN_REGULAR_SUBSYS_BASE.ipp"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

#define EXPLICIT_INSTANTIATION( T, D ) \
template class DOMAIN_REGULAR_CROSS_SUBSYS<T,D>;
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
