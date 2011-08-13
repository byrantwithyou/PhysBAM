//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_DIRICHLET_CONSTRAINT_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_DIRICHLET_CONSTRAINT_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
struct POST_EMBEDDING_INIT_DIRICHLET_CONSTRAINT_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POST_EMBEDDING_INIT_DIRICHLET_CONSTRAINT_VISITOR,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( DIRICHLET_CONSTRAINT_SYSTEM<T,D> )) &, constraint_system ))
        (( typename ARRAY<T>&, constraint_rhs ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        const int n_constraint = constraint_system.cell_linear_index_of_stencil_index.Size();
        constraint_system.Init_Stencil_Index_Of_Cell_Linear_Index();
        constraint_system.stencils.Exact_Resize(n_constraint);
        constraint_rhs.Exact_Resize(n_constraint);
    }
};

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_DIRICHLET_CONSTRAINT_VISITOR_HPP
