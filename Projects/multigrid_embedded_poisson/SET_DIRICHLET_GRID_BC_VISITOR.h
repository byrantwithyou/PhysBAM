//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_SET_DIRICHLET_GRID_BC_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_SET_DIRICHLET_GRID_BC_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template<
    class T_SYSTEM,
    class T_P_OF_INDEX,
    class T_RHS
>
struct SET_DIRICHLET_GRID_BC_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SET_DIRICHLET_GRID_BC_VISITOR,
        (( typename T_SYSTEM&, system ))
        (( typename T_P_OF_INDEX const, p_of_index ))
        (( typename T_RHS const, rhs ))
    )
public:
    typedef void result_type;
    template< class T_INDEX >
    void operator()(const T_INDEX& index) const
    { system.Set_Dirichlet_Grid_BC(index, p_of_index(index), rhs); }
};

template< class T_SYSTEM, class T_P_OF_INDEX, class T_RHS >
inline SET_DIRICHLET_GRID_BC_VISITOR< T_SYSTEM, T_P_OF_INDEX, T_RHS >
Make_Set_Dirichlet_Grid_BC_Visitor(
    T_SYSTEM& system,
    const T_P_OF_INDEX& p_of_index,
    const T_RHS& rhs)
{ return SET_DIRICHLET_GRID_BC_VISITOR< T_SYSTEM, T_P_OF_INDEX, T_RHS >(system, p_of_index, rhs); }

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_SET_DIRICHLET_GRID_BC_VISITOR_HPP
