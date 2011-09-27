//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_EMBEDDING_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_EMBEDDING_HPP

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_FWD.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< int D >
void Coarsen_Embedding(
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound,
    const HASHTABLE<int,int>& fine_embedding_index_of_linear_index,
    HASHTABLE<int,int>& coarse_embedding_index_of_linear_index,
    ARRAY<int>& coarse_linear_index_of_embedding_index);

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_EMBEDDING_HPP
