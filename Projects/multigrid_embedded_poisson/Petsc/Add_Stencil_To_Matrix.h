//#####################################################################
// Copyright 2010, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_ADD_STENCIL_TO_MATRIX_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_ADD_STENCIL_TO_MATRIX_HPP

#include <cassert>

#include <limits>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <petsc.h>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

namespace Petsc
{

template< class T_STENCIL_PROXY >
inline int
Add_Stencil_To_Matrix(
    const ARRAY_VIEW<int> petsc_index_of_physbam_index,
    T_STENCIL_PROXY stencil_proxy,
    const int value_index_of_row,
    std::vector< PetscInt >& petsc_matrix_col_indices,
    std::vector< PetscScalar >& petsc_matrix_values)
{
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, int >));
    typedef typename T_STENCIL_PROXY::INDEX_VALUE_TYPE INDEX_VALUE_TYPE;
    int value_index = value_index_of_row;
    BOOST_FOREACH( const INDEX_VALUE_TYPE index_value, stencil_proxy ) {
        const int petsc_index = petsc_index_of_physbam_index(index_value.index);
        assert(petsc_index != std::numeric_limits<int>::max());
        // Assert that the column indices within this row will be sorted.
        assert(
            value_index == value_index_of_row ||
            petsc_matrix_col_indices[value_index - 1] < petsc_index
        );
        petsc_matrix_col_indices[value_index] = petsc_index;
        petsc_matrix_values[value_index] = static_cast< PetscScalar >(index_value.value);
        ++value_index;
    }
    return value_index;
}

} // namespace Petsc

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #endif // PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_ADD_STENCIL_TO_MATRIX_HPP
