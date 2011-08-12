//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_OFFSET_GRID_BC_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_OFFSET_GRID_BC_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    int D,
    class T_Q_OF_FINE_INDEX_AND_NORMAL,
    class T_RHS_OF_INDEX
>
struct SET_NEUMANN_OFFSET_GRID_BC_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SET_NEUMANN_OFFSET_GRID_BC_VISITOR,
        (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
        (( typename T_Q_OF_FINE_INDEX_AND_NORMAL const, q_of_fine_index_and_normal ))
        (( typename T_RHS_OF_INDEX const, rhs_of_index ))
    )
public:
    typedef void result_type;
    void operator()(const VECTOR<int,D> index) const
    {
        VECTOR<int,D> fine_index = 2 * index - 1;
        VECTOR<int,D> normal; // init'ed to 0
        for(int d = 1; d <= D; ++d) {
            if(index[d] == 1) {
                --fine_index[d];
                normal[d] = -1;
                rhs_of_index(index) += q_of_fine_index_and_normal(fine_index, normal);
                normal[d] = 0;
                ++fine_index[d];
            }
            else if(index[d] == multi_index_bound.max_multi_index[d]) {
                ++fine_index[d];
                normal[d] = +1;
                rhs_of_index(index) += q_of_fine_index_and_normal(fine_index, normal);
                normal[d] = 0;
                --fine_index[d];
            }
        }
    }
};

template< int D, class T_Q_OF_FINE_INDEX_AND_NORMAL, class T_RHS_OF_INDEX >
inline SET_NEUMANN_OFFSET_GRID_BC_VISITOR< D, T_Q_OF_FINE_INDEX_AND_NORMAL, T_RHS_OF_INDEX >
Make_Set_Neumann_Offset_Grid_BC_Visitor(
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const T_Q_OF_FINE_INDEX_AND_NORMAL& q_of_fine_index_and_normal,
    const T_RHS_OF_INDEX& rhs_of_index)
{
    return SET_NEUMANN_OFFSET_GRID_BC_VISITOR< D, T_Q_OF_FINE_INDEX_AND_NORMAL, T_RHS_OF_INDEX >(
        multi_index_bound, q_of_fine_index_and_normal, rhs_of_index
    );
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SET_NEUMANN_OFFSET_GRID_BC_VISITOR_HPP
