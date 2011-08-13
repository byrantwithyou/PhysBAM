//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_DOMAIN_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_DOMAIN_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_DOMAIN_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POST_EMBEDDING_INIT_DOMAIN_VISITOR,
        (( typename T_EMBEDDING_SUBSYS&, embedding_subsys ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        const int n_embedding = embedding_subsys.linear_index_of_stencil_index.Size();
        embedding_subsys.Init_Stencil_Index_Of_Linear_Index();
        embedding_subsys.stencils.Exact_Resize(n_embedding, false); // uninit'ed
        embedding_subsys.Zero_Stencils();
    }
};

template< class T_EMBEDDING_SUBSYS >
inline POST_EMBEDDING_INIT_DOMAIN_VISITOR< T_EMBEDDING_SUBSYS >
Make_Post_Embedding_Init_Domain_Visitor(T_EMBEDDING_SUBSYS& embedding_subsys)
{ return POST_EMBEDDING_INIT_DOMAIN_VISITOR< T_EMBEDDING_SUBSYS >(embedding_subsys); }

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_POST_EMBEDDING_INIT_DOMAIN_VISITOR_HPP
