//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "DOMAIN_EMBEDDING_SUBSYS_BASE.ipp"
#include "DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

#define EXPLICIT_INSTANTIATION( T, D ) \
template class DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>;
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM
