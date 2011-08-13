//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_INTERFACE_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_INTERFACE_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D, class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_INTERFACE_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POST_EMBEDDING_INIT_INTERFACE_VISITOR,
        (( /******/ int const, n_material ))
        (( typename T_EMBEDDING_SUBSYS&, embedding_subsys ))
        (( typename ARRAY<T>&, system_rhs ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( INTERFACE_CONSTRAINT_SYSTEM<T,D> )) &, constraint_system ))
        (( typename ARRAY<T>&, constraint_rhs ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        const int n_virtual = embedding_subsys.index_of_stencil_index.Size();
        const int n_embedding = 2 * n_virtual;

        embedding_subsys.index_of_stencil_index.Preallocate(n_embedding);
        for(int i = 1; i <= n_virtual; ++i)
            embedding_subsys.index_of_stencil_index.Append(n_material + i);
        embedding_subsys.Init_Stencil_Index_Of_Index();
        embedding_subsys.stencils.Exact_Resize(n_embedding, false); // uninit'ed

        system_rhs.Exact_Resize(n_material + n_virtual); // init'ed to 0

        const int n_constraint = constraint_system.cell_index_of_stencil_index.Size();
        constraint_system.Init_Stencil_Index_Of_Cell_Index();
        constraint_system.stencils.Exact_Resize(n_constraint);
        constraint_system.stencils_containing_index.Initialize_New_Table(n_embedding);
        constraint_rhs.Exact_Resize(n_constraint);
    }
};

template< class T, int D, class T_EMBEDDING_SUBSYS >
inline POST_EMBEDDING_INIT_INTERFACE_VISITOR< T, D, T_EMBEDDING_SUBSYS >
Make_Post_Embedding_Init_Interface_Visitor(
    const int n_material,
    T_EMBEDDING_SUBSYS& embedding_subsys,
    ARRAY<T>& system_rhs,
    INTERFACE_CONSTRAINT_SYSTEM<T,D>& constraint_system,
    ARRAY<T>& constraint_rhs)
{
    return POST_EMBEDDING_INIT_INTERFACE_VISITOR< T, D, T_EMBEDDING_SUBSYS >(
        n_material,
        embedding_subsys, system_rhs,
        constraint_system, constraint_rhs
    );
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_INTERFACE_VISITOR_HPP
