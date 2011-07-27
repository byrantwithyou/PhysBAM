//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDING_POISSON_3D_V2_BUILD_DOMAIN_EMBEDDING_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDING_POISSON_3D_V2_BUILD_DOMAIN_EMBEDDING_SUBSYS_HPP

#include <cassert>

#include <iostream>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/Find_All_Signed_Cells_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Grid/Find_All_Signed_Vertices_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>

#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template<
    class T, int D,
    class T_SIGN_OF_CELL_INDEX,
    class T_POST_EMBEDDING_INIT_VISITOR,
    class T_EMBEDDING_CELL_VISITOR,
    class T_EMBEDDING_SUBSYS
>
int
Build_Domain_Embedding_Subsys(
    const MAIN_PARAMS<T,D>& main_params,
    T_SIGN_OF_CELL_INDEX sign_of_cell_index,
    T_POST_EMBEDDING_INIT_VISITOR post_embedding_init_visitor,
    T_EMBEDDING_CELL_VISITOR embedding_cell_visitor,
    T_EMBEDDING_SUBSYS& embedding_subsys,
    ARRAY<int>& embedding_cells)
{
    assert(embedding_cells.Size() == 0);

    BASIC_TIMER timer;

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;

    std::cout << "Initializing embedding cells...";
    std::cout.flush();
    timer.Restart();
    Find_All_Signed_Cells_Via_Cell_Sign_MT(
        main_params.general.n_thread,
        cell_multi_index_bound.Size(),
        0,
        sign_of_cell_index,
        embedding_cells
    );
    std::cout << timer.Elapsed() << " s" << std::endl;
    const int n_embedding_cells = embedding_cells.Size();
    std::cout << "  # of embedding cells = " << n_embedding_cells << std::endl;

    std::cout << "Initializing embedding indices in embedding system...";
    std::cout.flush();
    timer.Restart();
    Find_All_Signed_Vertices_Via_Cell_Sign_MT(
        main_params.general.n_thread,
        multi_index_bound,
        0,
        Make_Compose_Function(sign_of_cell_index, cell_multi_index_bound),
        embedding_subsys.linear_index_of_stencil_index
    );
    const int n_embedding = embedding_subsys.linear_index_of_stencil_index.Size();
    embedding_subsys.Init_Stencil_Index_Of_Linear_Index();
    embedding_subsys.stencils.Exact_Resize(n_embedding, false); // uninit'ed
    embedding_subsys.Zero_Stencils();
    std::cout << timer.Elapsed() << " s" << std::endl;
    std::cout << "  # of embedding vertices = " << n_embedding << std::endl;

    std::cout << "Initializing embedding system...";
    std::cout.flush();
    timer.Restart();
    post_embedding_init_visitor();
    for(int embedding_cell_index = 1; embedding_cell_index <= n_embedding_cells; ++embedding_cell_index) {
        const int cell_linear_index = embedding_cells(embedding_cell_index);
        assert(sign_of_cell_index(cell_linear_index) == 0);
        embedding_cell_visitor(cell_linear_index);
    }
    std::cout << timer.Elapsed() << " s" << std::endl;

    return 0;
}

template<
    class T, int D,
    class T_SIGN_OF_CELL_INDEX,
    class T_POST_EMBEDDING_INIT_VISITOR,
    class T_EMBEDDING_CELL_VISITOR,
    class T_EMBEDDING_SUBSYS
>
inline int
Build_Domain_Embedding_Subsys(
    const MAIN_PARAMS<T,D>& main_params,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    const T_POST_EMBEDDING_INIT_VISITOR& post_embedding_init_visitor,
    const T_EMBEDDING_CELL_VISITOR& embedding_cell_visitor,
    T_EMBEDDING_SUBSYS& embedding_subsys)
{
    ARRAY<int> embedding_cells;
    return Build_Domain_Embedding_Subsys(
        main_params,
        sign_of_cell_index,
        post_embedding_init_visitor,
        embedding_cell_visitor,
        embedding_subsys,
        embedding_cells
    );
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDING_POISSON_3D_V2_BUILD_DOMAIN_EMBEDDING_SUBSYS_HPP
