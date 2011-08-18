//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_AS_SPARSE_MATRIX_FLAT_NXN_H
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_AS_SPARSE_MATRIX_FLAT_NXN_H

#include <algorithm>
#include <limits>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Algorithm/Fill.h>
#include <Jeffrey_Utilities/Algorithm/Find_All_If.h>
#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>

namespace PhysBAM
{

namespace Detail_As_Sparse_Matrix_Flat_NXN
{

template< class T_SYSTEM >
struct HAS_NONZERO_STENCIL;
struct SET_FLAT_INDEX_OF_ORIG_INDEX;

} // namespace Detail_As_Sparse_Matrix_Flat_NXN

template< class T, class T_SYSTEM >
void
As_Sparse_Matrix_Flat_NXN(
    const unsigned int n_thread,
    const T_SYSTEM& system,
    SPARSE_MATRIX_FLAT_NXN<T>& flat_system,
    ARRAY<int>& orig_index_of_flat_index,
    ARRAY<int>& flat_index_of_orig_index)
{
    typedef Detail_As_Sparse_Matrix_Flat_NXN::HAS_NONZERO_STENCIL< T_SYSTEM > HAS_NONZERO_STENCIL_;
    typedef Detail_As_Sparse_Matrix_Flat_NXN::SET_FLAT_INDEX_OF_ORIG_INDEX SET_FLAT_INDEX_OF_ORIG_INDEX_;

    const int n_orig_index = flat_index_of_orig_index.Size();
    orig_index_of_flat_index.Remove_All();
    Find_All_If_MT(
        n_thread,
        1, n_orig_index,
        HAS_NONZERO_STENCIL_(system),
        orig_index_of_flat_index
    );
    const int n_flat_index = orig_index_of_flat_index.Size();
    Fill_MT(
        n_thread,
        flat_index_of_orig_index,
        std::numeric_limits<int>::max()
    );
    For_Each_MT(
        n_thread,
        0, n_flat_index - 1,
        SET_FLAT_INDEX_OF_ORIG_INDEX_(
            orig_index_of_flat_index,
            flat_index_of_orig_index
        )
    );

    int matrix_nnz = 0;
    int max_stencil_nnz = 0;
    for(int flat_index = 1; flat_index <= n_flat_index; ++flat_index) {
        const int stencil_nnz = system.Stencil_N_Nonzero(flat_index);
        matrix_nnz += stencil_nnz;
        max_stencil_nnz = std::max(max_stencil_nnz, stencil_nnz);
    }
    flat_system.Reset();
    flat_system.offsets.Preallocate(n_flat_index);
    flat_system.A.Preallocate(matrix_nnz);
    {
        int value_index = 0;
        UNSTRUCTURED_STENCIL<int,T> system_stencil;
        system_stencil.values.Preallocate(max_stencil_nnz);
        UNSTRUCTURED_STENCIL_PROXY< UNSTRUCTURED_STENCIL<int,T> > system_stencil_proxy(system_stencil);
        for(int flat_index = 1; flat_index <= n_flat_index; ++flat_index) {
            const int orig_index = orig_index_of_flat_index(flat_index);
            system_stencil.values.Remove_All();
            system.Add_Stencil_To(orig_index, system_stencil_proxy);
            typedef INDEX_VALUE<int,T> INDEX_VALUE_TYPE;
            BOOST_FOREACH( const INDEX_VALUE_TYPE index_value, system_stencil_proxy ) {
                const int other_flat_index = flat_index_of_orig_index(index_value.index);
                assert(other_flat_index != std::numeric_limits<int>::max());
                flat_system.Append_Entry_To_Current_Row(other_flat_index, index_value.value);
            }
            flat_system.Finish_Row();
        }
    }
}

namespace Detail_As_Sparse_Matrix_Flat_NXN
{

template< class T_SYSTEM >
struct HAS_NONZERO_STENCIL
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        HAS_NONZERO_STENCIL,
        (( typename T_SYSTEM const &, system ))
    )
public:
    typedef bool result_type;
    bool operator()(const int index) const
    { return system.Stencil_N_Nonzero(index) != 0; }
};

struct SET_FLAT_INDEX_OF_ORIG_INDEX
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SET_FLAT_INDEX_OF_ORIG_INDEX,
        (( ARRAY<int> const &, orig_index_of_flat_index ))
        (( ARRAY<int> /***/ &, flat_index_of_orig_index ))
    )
public:
    typedef void result_type;
    void operator()(const int flat_index) const
    {
        const int orig_index = orig_index_of_flat_index(flat_index);
        flat_index_of_orig_index(orig_index) = flat_index;
    }
};

} // namespace Detail_As_Sparse_Matrix_Flat_NXN

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_AS_SPARSE_MATRIX_FLAT_NXN_H
