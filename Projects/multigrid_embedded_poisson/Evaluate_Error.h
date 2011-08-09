//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EVALUATE_ERROR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EVALUATE_ERROR_HPP

#include <cassert>
#include <cmath>

#include <algorithm>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_X.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/Safe_Dv_Over_Dx.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_SIGN_OF_INDEX,
    class T_SIGN_OF_CELL_INDEX
>
void
Evaluate_Error(
    typename EXAMPLE_PARAMS<T,D>::DOMAIN_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    const int sign,
    T_SIGN_OF_INDEX sign_of_index,
    T_SIGN_OF_CELL_INDEX sign_of_cell_index,
    const ARRAY_VIEW<const T> u_approx,
    T& max_u_error, T& max_grad_u_error)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    assert(sign == -1 || sign == +1);

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = As_Vector<int>(main_params.grid.n_cell);
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

    const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
    const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);
    const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;

    assert(u_approx.Size() == multi_index_bound.Size());

    for(int linear_index = 1; linear_index <= multi_index_bound.Size(); ++linear_index) {
        if(sign_of_index(linear_index) != sign)
            continue;
        const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
        const VECTOR<T,D> x = Multi_Index_X(min_x, max_x, multi_index_bound, multi_index);
        const T u_continuous_x = problem.u(x);
        const T u_approx_x = u_approx(linear_index);
        max_u_error = std::max(max_u_error, std::abs(u_continuous_x - u_approx_x));
        const VECTOR<T,D> grad_u_continuous_x = problem.grad_u(x);
        VECTOR<T,D> grad_u_approx_x; // init'ed to 0
        int count = 0;
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE cell_multi_index,
            Multi_Index_Box_Intersect(MULTI_INDEX_CUBE<D,-1,0>(multi_index), cell_multi_index_bound)
        ) {
            if(sign_of_cell_index(cell_multi_index) == 0)
                continue;
            ++count;
            const MULTI_INDEX_TYPE multi_offset = 1 + 2 * (cell_multi_index - multi_index);
            MULTI_INDEX_TYPE other_multi_index = multi_index;
            for(int d = 1; d <= D; ++d) {
                other_multi_index[d] += multi_offset[d];
                assert(multi_index_bound.Contains(other_multi_index));
                const int other_linear_index = multi_index_bound.Linear_Index(other_multi_index);
                grad_u_approx_x[d] += multi_offset[d] * (u_approx(other_linear_index) - u_approx_x) / dx[d];
                other_multi_index[d] = multi_index[d];
            }
        }
        if(count != 0) {
            grad_u_approx_x /= static_cast<T>(count);
            max_grad_u_error = std::max(max_grad_u_error, (grad_u_continuous_x - grad_u_approx_x).Max_Abs());
        }
    }
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EVALUATE_ERROR_HPP
