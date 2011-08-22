//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_FULL_STENCILS_ON_GRID_BOUNDARY_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_FULL_STENCILS_ON_GRID_BOUNDARY_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/RESULT_OF.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_SYSTEM, class T_BETA_OF_CELL_INDEX >
struct INIT_FULL_STENCILS_ON_GRID_BOUNDARY_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INIT_FULL_STENCILS_ON_GRID_BOUNDARY_VISITOR,
        (( typename T_SYSTEM&, system ))
        (( typename T_BETA_OF_CELL_INDEX const, beta_of_cell_index ))
    )
public:
    typedef void result_type;
    template< int D >
    void operator()(const VECTOR<int,D> outside_cell_multi_index) const
    {
        typedef typename RESULT_OF< const T_BETA_OF_CELL_INDEX ( VECTOR<int,D> ) >::type BETA_TYPE;
        BETA_TYPE beta = beta_of_cell_index(outside_cell_multi_index);
        system.Init_Full_Stencils_On_Grid_Boundary(outside_cell_multi_index, beta);
    }
};

template< class T_SYSTEM, class T_BETA_OF_CELL_INDEX >
inline INIT_FULL_STENCILS_ON_GRID_BOUNDARY_VISITOR< T_SYSTEM, T_BETA_OF_CELL_INDEX >
Make_Init_Full_Stencils_On_Grid_Boundary_Visitor(
    T_SYSTEM& system,
    const T_BETA_OF_CELL_INDEX& beta_of_cell_index)
{
    return INIT_FULL_STENCILS_ON_GRID_BOUNDARY_VISITOR<
        T_SYSTEM, T_BETA_OF_CELL_INDEX
    >(system, beta_of_cell_index);
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_FULL_STENCILS_ON_GRID_BOUNDARY_VISITOR_HPP
