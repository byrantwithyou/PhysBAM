//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SELECT_INDYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SELECT_INDYS_HPP

#include <cassert>
#include <cmath>

#include <algorithm>

#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/RESULT_OF.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace Detail_Aggregate_Constraints
{

template< class T >
struct HAS_GREATER_WEIGHT;

} // namespace Detail_Aggregate_Constraints

template<
    int D,
    class T_MULTI_INDEX_OF_INDEX,
    class T_CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_INDEX,
    class T_INDYABLE
>
void Select_Indys(
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound,
    T_MULTI_INDEX_OF_INDEX multi_index_of_index,
    const HASHTABLE<int,int>& constraint_index_of_cell_linear_index,
    T_CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_INDEX constraint_stencil_proxy_of_constraint_index,
    T_INDYABLE indyable,
    const float min_relative_weight,
    ARRAY_VIEW<int> indy_index_of_constraint_index,
    ARRAY<int>& index_of_indy_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    typedef typename RESULT_OF<
        T_CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_INDEX ( int )
    >::type CONSTRAINT_STENCIL_PROX_TYPE;
    BOOST_MPL_ASSERT((boost::is_convertible< typename CONSTRAINT_STENCIL_PROX_TYPE::INDEX_TYPE, int >));
    typedef typename CONSTRAINT_STENCIL_PROX_TYPE::SCALAR_TYPE SCALAR_TYPE;
    typedef INDEX_VALUE< int, SCALAR_TYPE > INDEX_VALUE_TYPE;

    assert(0 <= min_relative_weight && min_relative_weight <= 1);

    const int n_constraint = constraint_index_of_cell_linear_index.Size();
    assert(n_constraint == constraint_index_of_cell_linear_index.Size());
    assert(n_constraint == indy_index_of_constraint_index.Size());
    if(n_constraint == 0)
        return;

    static const int INVALID_INDY_INDEX = 0;
    indy_index_of_constraint_index.Fill(INVALID_INDY_INDEX);

    // Compute the weight of each indyable.
    HASHTABLE< int, SCALAR_TYPE > weight_of_index(3 * n_constraint);
    for(int constraint_index = 1; constraint_index <= n_constraint; ++constraint_index) {
        BOOST_FOREACH(
            const INDEX_VALUE_TYPE index_value,
            constraint_stencil_proxy_of_constraint_index(constraint_index)
        ) {
            if(index_value.value != 0 && indyable(index_value.index))
                weight_of_index.Get_Or_Insert(index_value.index) += std::abs(index_value.value);
        }
    }

    // Sort indyables by weight.
    const int n_indyable = weight_of_index.Size();
    ARRAY<int> index_of_indyable_index;
    weight_of_index.Get_Keys(index_of_indyable_index);
    std::sort(
        index_of_indyable_index.begin(), index_of_indyable_index.end(),
        Detail_Aggregate_Constraints::HAS_GREATER_WEIGHT< SCALAR_TYPE >(weight_of_index)
    );
    const SCALAR_TYPE max_weight = weight_of_index.Get(index_of_indyable_index(1));

    // Select indys and associate cells to an indy.
    int n_covered_constraint = 0;
    for(int indyable_index = 1; indyable_index <= n_indyable; ++indyable_index) {

        assert(n_covered_constraint <= n_constraint);
        const int index = index_of_indyable_index(indyable_index);
        const SCALAR_TYPE* const p_weight = weight_of_index.Get_Pointer(index);
        if(!p_weight)
            continue;
        const SCALAR_TYPE weight = *p_weight;
        assert(weight <= max_weight);
        if(
            weight == 0
         || (n_covered_constraint == n_constraint && weight <= min_relative_weight * max_weight)
        )
            break;
        index_of_indy_index.Append(index);

        const MULTI_INDEX_TYPE multi_index = multi_index_of_index(index);
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE cell_multi_index,
            Multi_Index_Box_Intersect(MULTI_INDEX_CUBE<D,-1,0>(multi_index), cell_multi_index_bound)
        ) {
            const int cell_linear_index = cell_multi_index_bound.Linear_Index(cell_multi_index);
            const int* const p_constraint_index = constraint_index_of_cell_linear_index.Get_Pointer(cell_linear_index);
            if(!p_constraint_index)
                continue;
            const int constraint_index = *p_constraint_index;
            BOOST_FOREACH(
                const INDEX_VALUE_TYPE index_value,
                constraint_stencil_proxy_of_constraint_index(constraint_index)
            )
                weight_of_index.Delete_If_Present(index_value.index);
        }

        BOOST_FOREACH(
            const MULTI_INDEX_TYPE cell_multi_index,
            Multi_Index_Box_Intersect(MULTI_INDEX_CUBE<D,-2,+1>(multi_index), cell_multi_index_bound)
        ) {
            const int cell_linear_index = cell_multi_index_bound.Linear_Index(cell_multi_index);
            const int* const p_constraint_index = constraint_index_of_cell_linear_index.Get_Pointer(cell_linear_index);
            if(!p_constraint_index)
                continue;
            const int constraint_index = *p_constraint_index;
            int& indy_index = indy_index_of_constraint_index(constraint_index);
            if(indy_index == INVALID_INDY_INDEX) {
                assert(n_covered_constraint < n_constraint);
                ++n_covered_constraint;
            }
            else {
                const MULTI_INDEX_TYPE multi_offset = cell_multi_index - multi_index;
                if(multi_offset.Min() < -1 || 0 < multi_offset.Max())
                    continue;
            }
            indy_index = index_of_indy_index.Size();
        }

    }

    assert(!indy_index_of_constraint_index.Contains(INVALID_INDY_INDEX));
}

namespace Detail_Aggregate_Constraints
{

template< class T >
struct HAS_GREATER_WEIGHT
{
    const HASHTABLE<int,T>& weight_of_index;

    HAS_GREATER_WEIGHT(const HASHTABLE<int,T>& weight_of_index_)
        : weight_of_index(weight_of_index_)
    { }

    typedef bool result_type;

    bool operator()(const int index0, const int index1) const
    {
        const T weight0 = weight_of_index.Get(index0);
        const T weight1 = weight_of_index.Get(index1);
        return weight0 > weight1 || (weight0 == weight1 && index0 < index1);
    }
};

} // namespace Detail_Aggregate_Constraints

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SELECT_INDYS_HPP
