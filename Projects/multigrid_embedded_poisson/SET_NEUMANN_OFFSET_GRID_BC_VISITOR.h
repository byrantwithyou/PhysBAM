//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_OFFSET_GRID_BC_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_OFFSET_GRID_BC_VISITOR_HPP

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_SYSTEM,
    class T_BETA_OF_CELL_INDEX,
    class T_Q_OF_CELL_INDEX_AND_NORMAL,
    class T_RHS_OF_INDEX
>
struct SET_NEUMANN_OFFSET_GRID_BC_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SET_NEUMANN_OFFSET_GRID_BC_VISITOR,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, dx ))
        (( typename MULTI_INDEX_BOUND<D> const, cell_multi_index_bound ))
        (( typename T_SYSTEM&, system ))
        (( typename T_BETA_OF_CELL_INDEX const, beta_of_cell_index ))
        (( typename T_Q_OF_CELL_INDEX_AND_NORMAL const, q_of_cell_index_and_normal ))
        (( typename T_RHS_OF_INDEX const, rhs_of_index ))
    )
public:
    typedef void result_type;
    void operator()(const VECTOR<int,D> outside_cell_multi_index) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;
        const T beta = beta_of_cell_index(outside_cell_multi_index);
        system.Set_Pure_Neumann_Offset_Grid_BC(outside_cell_multi_index, beta);
        const T dv = dx.Product();
        const MULTI_INDEX_TYPE clamped_cell_multi_index = cell_multi_index_bound.Clamp(outside_cell_multi_index);
        assert((clamped_cell_multi_index - outside_cell_multi_index).Max_Abs() == 1);
        VECTOR<T,D> normal; // init'ed to 0
        for(int d = 1; d <= D; ++d) {
            normal[d] = static_cast<T>(outside_cell_multi_index[d] - clamped_cell_multi_index[d]);
            if(normal[d] == 0)
                continue;
            const T q = q_of_cell_index_and_normal(outside_cell_multi_index, normal);
            BOOST_FOREACH(
                const MULTI_INDEX_TYPE multi_index,
                Multi_Index_Box_Intersect(
                    MULTI_INDEX_CUBE<D,0,1>(outside_cell_multi_index),
                    MULTI_INDEX_CUBE<D,0,1>(clamped_cell_multi_index)
                )
            )
                rhs_of_index(multi_index) += dv * q / ((1 << (D-1)) * dx[d]);
            normal[d] = static_cast<T>(0);
        }
    }
};

template<
    class T, int D,
    class T_SYSTEM,
    class T_BETA_OF_CELL_INDEX,
    class T_Q_OF_CELL_INDEX_AND_NORMAL,
    class T_RHS_OF_INDEX
>
inline SET_NEUMANN_OFFSET_GRID_BC_VISITOR<
    T, D,
    T_SYSTEM,
    T_BETA_OF_CELL_INDEX,
    T_Q_OF_CELL_INDEX_AND_NORMAL,
    T_RHS_OF_INDEX
>
Make_Set_Neumann_Offset_Grid_BC_Visitor(
    const VECTOR<T,D>& dx,
    const MULTI_INDEX_BOUND<D>& cell_multi_index_bound,
    T_SYSTEM& system,
    const T_BETA_OF_CELL_INDEX& beta_of_cell_index,
    const T_Q_OF_CELL_INDEX_AND_NORMAL& q_of_cell_index_and_normal,
    const T_RHS_OF_INDEX& rhs_of_index)
{
    return SET_NEUMANN_OFFSET_GRID_BC_VISITOR<
        T, D,
        T_SYSTEM,
        T_BETA_OF_CELL_INDEX,
        T_Q_OF_CELL_INDEX_AND_NORMAL,
        T_RHS_OF_INDEX
    >(
        dx, cell_multi_index_bound,
        system,
        beta_of_cell_index,
        q_of_cell_index_and_normal,
        rhs_of_index
    );
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_OFFSET_GRID_BC_VISITOR_HPP
