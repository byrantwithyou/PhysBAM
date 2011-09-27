//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HPP

#include <cassert>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

namespace PhysBAM
{

template< class T_SIGN_OF_CELL_INDEX, class T_CELL_VISITOR >
inline void
Visit_Cells_With_Sign_Via_Cell_Sign(
    const int n_cell,
    T_SIGN_OF_CELL_INDEX sign_of_cell_index,
    T_CELL_VISITOR cell_visitor)
{
    for(int cell_linear_index = 1; cell_linear_index <= n_cell; ++cell_linear_index) {
        const int cell_sign = sign_of_cell_index(cell_linear_index);
        cell_visitor(cell_linear_index, cell_sign);
    }
}

//#####################################################################
//#####################################################################

namespace Detail_Visit_Cells_With_Sign_Via_Cell_Sign
{

template< class T_SIGN_OF_CELL_INDEX, class T_CELL_VISITOR >
struct VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HELPER;

} // namespace Detail_Visit_Cells_With_Sign_Via_Cell_Sign

template< class T_SIGN_OF_CELL_INDEX, class T_CELL_VISITOR >
inline void
Visit_Cells_With_Sign_Via_Cell_Sign_MT(
    const unsigned int n_thread,
    const int n_cell,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    const T_CELL_VISITOR& cell_visitor)
{
    typedef Detail_Visit_Cells_With_Sign_Via_Cell_Sign::
        VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HELPER<
            T_SIGN_OF_CELL_INDEX, T_CELL_VISITOR
        > VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HELPER_;
    assert(n_thread >= 1);
    For_Each_MT(
        n_thread,
        1, n_cell,
        VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HELPER_(
            sign_of_cell_index, cell_visitor
        )
    );
}

namespace Detail_Visit_Cells_With_Sign_Via_Cell_Sign
{

template< class T_SIGN_OF_CELL_INDEX, class T_CELL_VISITOR >
struct VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HELPER,
        (( typename T_SIGN_OF_CELL_INDEX const, sign_of_cell_index ))
        (( typename T_CELL_VISITOR const, cell_visitor ))
    )
public:
    typedef void result_type;
    void operator()(const int cell_linear_index) const
    {
        const int cell_sign = sign_of_cell_index(cell_linear_index);
        cell_visitor(cell_linear_index, cell_sign);
    }
};

} // namespace Detail_Visit_Cells_With_Sign_Via_Cell_Sign

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_CELLS_WITH_SIGN_VIA_CELL_SIGN_HPP
