//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_GRID_BC_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_GRID_BC_VISITOR_HPP

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_Q_OF_CELL_INDEX_AND_NORMAL,
    class T_F_OF_INDEX,
    class T_RHS_OF_INDEX
>
struct SET_NEUMANN_GRID_BC_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SET_NEUMANN_GRID_BC_VISITOR,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, dx ))
        (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
        (( typename T_Q_OF_CELL_INDEX_AND_NORMAL const, q_of_cell_index_and_normal ))
        (( typename T_F_OF_INDEX const, f_of_index ))
        (( typename T_RHS_OF_INDEX const, rhs_of_index ))
    )
public:
    typedef void result_type;
    void operator()(const VECTOR<int,D> outside_cell_multi_index) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;
        const MULTI_INDEX_BOUND<D> cell_multi_index_bound = multi_index_bound - 1;
        const T dv = dx.Product();

        BOOST_FOREACH(
            const MULTI_INDEX_TYPE multi_index,
            (MULTI_INDEX_CUBE<D,0,1>(outside_cell_multi_index))
        )
            if(multi_index_bound.Contains(multi_index))
                rhs_of_index(multi_index) += f_of_index(multi_index) * (dv / (1 << D));

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
    class T_Q_OF_CELL_INDEX_AND_NORMAL,
    class T_F_OF_INDEX,
    class T_RHS_OF_INDEX
>
inline SET_NEUMANN_GRID_BC_VISITOR<
    T, D,
    T_Q_OF_CELL_INDEX_AND_NORMAL,
    T_F_OF_INDEX,
    T_RHS_OF_INDEX
>
Make_Set_Neumann_Grid_BC_Visitor(
    const VECTOR<T,D>& dx,
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const T_Q_OF_CELL_INDEX_AND_NORMAL& q_of_cell_index_and_normal,
    const T_F_OF_INDEX& f_of_index,
    const T_RHS_OF_INDEX& rhs_of_index)
{
    return SET_NEUMANN_GRID_BC_VISITOR<
        T, D,
        T_Q_OF_CELL_INDEX_AND_NORMAL,
        T_F_OF_INDEX,
        T_RHS_OF_INDEX
    >(
        dx, multi_index_bound,
        q_of_cell_index_and_normal,
        f_of_index,
        rhs_of_index
    );
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_GRID_BC_VISITOR_HPP
