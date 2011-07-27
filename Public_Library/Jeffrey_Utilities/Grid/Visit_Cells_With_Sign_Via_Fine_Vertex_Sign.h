//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HPP

#include <cassert>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Grid/Cell_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

namespace PhysBAM
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX, class T_CELL_VISITOR >
inline void
Visit_Cells_With_Sign_Via_Fine_Vertex_Sign(
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound,
    const T_SIGN_OF_FINE_INDEX& sign_of_fine_index,
    T_CELL_VISITOR cell_visitor,
    const int sign_of_zero = -1)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    const int n_cell = cell_multi_index_bound.Size();
    for(int cell_linear_index = 1; cell_linear_index <= n_cell; ++cell_linear_index) {
        const MULTI_INDEX_TYPE cell_multi_index = cell_multi_index_bound.Multi_Index(cell_linear_index);
        const int cell_sign =
            Cell_Sign_Via_Fine_Vertex_Sign< FINE_FACTOR >(sign_of_fine_index, cell_multi_index, sign_of_zero);
        cell_visitor(cell_linear_index, cell_sign);
    }
}

//#####################################################################
//#####################################################################

namespace Detail_Visit_Cells_With_Sign_Via_Fine_Vertex_Sign
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX, class T_CELL_VISITOR >
struct VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER;

} // namespace Detail_Visit_Cells_With_Sign_Via_Fine_Vertex_Sign

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX, class T_CELL_VISITOR >
inline void
Visit_Cells_With_Sign_Via_Fine_Vertex_Sign_MT(
    const unsigned int n_thread,
    const MULTI_INDEX_BOUND<D>& cell_multi_index_bound,
    const T_SIGN_OF_FINE_INDEX& sign_of_fine_index,
    const T_CELL_VISITOR& cell_visitor,
    const int sign_of_zero = -1)
{
    typedef Detail_Visit_Cells_With_Sign_Via_Fine_Vertex_Sign::
        VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER<
            FINE_FACTOR, D, T_SIGN_OF_FINE_INDEX, T_CELL_VISITOR
        > VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER_;
    assert(n_thread >= 1);
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    For_Each_MT(
        n_thread,
        1, cell_multi_index_bound.Size(),
        VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER_(
            cell_multi_index_bound, sign_of_fine_index, cell_visitor, sign_of_zero
        )
    );
}

namespace Detail_Visit_Cells_With_Sign_Via_Fine_Vertex_Sign
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX, class T_CELL_VISITOR >
struct VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER,
        (( typename MULTI_INDEX_BOUND<D> const, cell_multi_index_bound ))
        (( typename T_SIGN_OF_FINE_INDEX const, sign_of_fine_index ))
        (( typename T_CELL_VISITOR const, cell_visitor ))
        (( /******/ int const, sign_of_zero ))
    )
public:
    typedef void result_type;
    void operator()(const int cell_linear_index) const
    {
        const VECTOR<int,D> cell_multi_index = cell_multi_index_bound.Multi_Index(cell_linear_index);
        const int cell_sign = Cell_Sign_Via_Fine_Vertex_Sign< FINE_FACTOR >(
            sign_of_fine_index, cell_multi_index, sign_of_zero
        );
        cell_visitor(cell_linear_index, cell_sign);
    }
};

} // namespace Detail_Visit_Cells_With_Sign_Via_Fine_Vertex_Sign

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_CELLS_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HPP
