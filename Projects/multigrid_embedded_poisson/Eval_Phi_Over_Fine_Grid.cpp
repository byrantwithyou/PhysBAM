//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Jeffrey_Utilities/Eval_Grid_Function.h>
#include <Jeffrey_Utilities/Level_Sets/Shift_Level_Set_Away_From_Vertices.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Params/MAIN_PARAMS.h"

#include "Eval_Phi_Over_Fine_Grid.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
void
Eval_Phi_Over_Fine_Grid(
    const MAIN_PARAMS<T,D>& main_params,
    ARRAY_VIEW<T> phi_of_fine_index,
    const int sign_of_zero /*= -1*/)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

    const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
    const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);

    Eval_Grid_Function_MT(
        main_params.general.n_thread,
        min_x, max_x, fine_multi_index_bound,
        main_params.level_set.phi,
        phi_of_fine_index
    );

    if(main_params.level_set.min_dist_to_vertex != 0)
        Shift_Level_Set_Away_From_Vertices_MT(
            main_params.general.n_thread,
            fine_multi_index_bound,
            phi_of_fine_index,
            main_params.level_set.min_dist_to_vertex,
            sign_of_zero
        );
}

#define EXPLICIT_INSTANTIATION( T, D ) \
template void \
Eval_Phi_Over_Fine_Grid<T,D>( \
    const MAIN_PARAMS<T,D>& main_params, \
    ARRAY_VIEW<T> phi_of_fine_index, \
    const int sign_of_zero);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
