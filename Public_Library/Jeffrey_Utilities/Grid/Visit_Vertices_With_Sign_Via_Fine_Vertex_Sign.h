//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HPP

#include <cassert>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Grid/Vertex_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

namespace PhysBAM
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX, class T_VERTEX_VISITOR >
inline void
Visit_Vertices_With_Sign_Via_Fine_Vertex_Sign(
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    const T_SIGN_OF_FINE_INDEX& sign_of_fine_index,
    T_VERTEX_VISITOR vertex_visitor,
    const int sign_of_zero = -1)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    const int n = multi_index_bound.Size();
    for(int linear_index = 1; linear_index <= n; ++linear_index) {
        const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
        const int sign =
            Vertex_Sign_Via_Fine_Vertex_Sign< FINE_FACTOR >(
                multi_index_bound, sign_of_fine_index, multi_index, sign_of_zero
            );
        vertex_visitor(linear_index, sign);
    }
}

//#####################################################################
//#####################################################################

namespace Detail_Visit_Vertices_With_Sign_Via_Fine_Vertex_Sign
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX, class T_VERTEX_VISITOR >
struct VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER;

} // namespace Detail_Visit_Vertices_With_Sign_Via_Fine_Vertex_Sign

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX, class T_VERTEX_VISITOR >
inline void
Visit_Vertices_With_Sign_Via_Fine_Vertex_Sign_MT(
    const unsigned int n_thread,
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const T_SIGN_OF_FINE_INDEX& sign_of_fine_index,
    const T_VERTEX_VISITOR& vertex_visitor,
    const int sign_of_zero = -1)
{
    typedef Detail_Visit_Vertices_With_Sign_Via_Fine_Vertex_Sign::
        VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER<
            FINE_FACTOR, D, T_SIGN_OF_FINE_INDEX, T_VERTEX_VISITOR
        > VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER_;
    assert(n_thread >= 1);
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    For_Each_MT(
        n_thread,
        1, multi_index_bound.Size(),
        VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER_(
            multi_index_bound, sign_of_fine_index, vertex_visitor, sign_of_zero
        )
    );
}

namespace Detail_Visit_Vertices_With_Sign_Via_Fine_Vertex_Sign
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX, class T_VERTEX_VISITOR >
struct VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HELPER,
        (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
        (( typename T_SIGN_OF_FINE_INDEX const, sign_of_fine_index ))
        (( typename T_VERTEX_VISITOR const, vertex_visitor ))
        (( /******/ int const, sign_of_zero ))
    )
private:
    typedef void result_type;
    void operator()(const int linear_index) const
    {
        const VECTOR<int,D> multi_index = multi_index_bound.Multi_Index(linear_index);
        const int sign = Vertex_Sign_Via_Fine_Vertex_Sign< FINE_FACTOR >(
            multi_index_bound, sign_of_fine_index, multi_index, sign_of_zero
        );
        vertex_visitor(linear_index, sign);
    }
};

} // namespace Detail_Visit_Vertices_With_Sign_Via_Fine_Vertex_Sign

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_VERTICES_WITH_SIGN_VIA_FINE_VERTEX_SIGN_HPP
