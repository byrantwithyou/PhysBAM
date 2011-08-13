//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDING_POISSON_3D_V2_BUILD_EMBEDDING_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDING_POISSON_3D_V2_BUILD_EMBEDDING_SUBSYS_HPP

#include <cassert>

#include <iosfwd>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/Find_All_Signed_Cells_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Grid/Find_All_Signed_Vertices_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    int D,
    class T_SIGN_OF_CELL_INDEX,
    class T_POST_EMBEDDING_INIT_VISITOR,
    class T_EMBEDDING_CELL_VISITOR
>
int
Build_Embedding_Subsys(
    const unsigned int n_thread,
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    T_SIGN_OF_CELL_INDEX sign_of_cell_index,
    T_POST_EMBEDDING_INIT_VISITOR post_embedding_init_visitor,
    T_EMBEDDING_CELL_VISITOR embedding_cell_visitor,
    ARRAY<int>& index_of_stencil_index,
    ARRAY<int>& embedding_cells,
    std::ostream& lout = PhysBAM::nout)
{
    assert(embedding_cells.Size() == 0);

    BASIC_TIMER timer;

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = multi_index_bound - 1;

    lout << "Initializing embedding cells...";
    lout.flush();
    timer.Restart();
    Find_All_Signed_Cells_Via_Cell_Sign_MT(
        n_thread,
        cell_multi_index_bound.Size(),
        0,
        sign_of_cell_index,
        embedding_cells
    );
    lout << timer.Elapsed() << " s" << std::endl;
    const int n_embedding_cells = embedding_cells.Size();
    lout << "  # of embedding cells = " << n_embedding_cells << std::endl;

    lout << "Initializing embedding indices in embedding system...";
    lout.flush();
    timer.Restart();
    Find_All_Signed_Vertices_Via_Cell_Sign_MT(
        n_thread,
        multi_index_bound,
        0,
        Make_Compose_Function(sign_of_cell_index, cell_multi_index_bound),
        index_of_stencil_index
    );
    lout << timer.Elapsed() << " s" << std::endl;
    const int n_embedding = index_of_stencil_index.Size();
    lout << "  # of embedding vertices = " << n_embedding << std::endl;

    lout << "Initializing embedding system...";
    lout.flush();
    timer.Restart();
    post_embedding_init_visitor();
    // TODO: MT
    BOOST_FOREACH( const int cell_linear_index, embedding_cells ) {
        assert(sign_of_cell_index(cell_linear_index) == 0);
        embedding_cell_visitor(cell_linear_index);
    }
    lout << timer.Elapsed() << " s" << std::endl;

    return 0;
}

template<
    int D,
    class T_SIGN_OF_CELL_INDEX,
    class T_POST_EMBEDDING_INIT_VISITOR,
    class T_EMBEDDING_CELL_VISITOR
>
inline int
Build_Embedding_Subsys(
    const unsigned int n_thread,
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    const T_POST_EMBEDDING_INIT_VISITOR& post_embedding_init_visitor,
    const T_EMBEDDING_CELL_VISITOR& embedding_cell_visitor,
    ARRAY<int>& index_of_stencil_index,
    std::ostream& lout = PhysBAM::nout)
{
    ARRAY<int> embedding_cells;
    return Build_Embedding_Subsys(
        n_thread,
        multi_index_bound,
        sign_of_cell_index,
        post_embedding_init_visitor,
        embedding_cell_visitor,
        index_of_stencil_index,
        embedding_cells,
        lout
    );
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDING_POISSON_3D_V2_BUILD_EMBEDDING_SUBSYS_HPP
