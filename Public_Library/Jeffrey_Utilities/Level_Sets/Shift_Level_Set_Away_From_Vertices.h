//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HPP

#include <cassert>
#include <cmath>

#include <limits>

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

template< class T, class T_MULTI_INDEX_BOX >
inline void
Shift_Level_Set_Value_Away_From_Vertices(
    const T_MULTI_INDEX_BOX& multi_index_box,
    ARRAY_VIEW<T> phi_of_index,
    const float min_dist_to_vertex,
    const int linear_index,
    const int sign_of_zero = -1)
{
    static const int DIMENSION = T_MULTI_INDEX_BOX::DIMENSION;
    typedef VECTOR< int, DIMENSION > MULTI_INDEX_TYPE;

    static const T one_p_eps = 1 + std::numeric_limits<T>::epsilon();
    static const T one_m_eps = 1 - std::numeric_limits<T>::epsilon();

    assert(phi_of_index.Size() == multi_index_box.Size());
    assert(0 <= min_dist_to_vertex && min_dist_to_vertex <= 0.5f);
    assert(sign_of_zero == -1 || sign_of_zero == +1);

    const MULTI_INDEX_TYPE multi_index = multi_index_box.Multi_Index(linear_index);

    T& phi = phi_of_index(linear_index);
    if(phi == 0)
        phi = sign_of_zero * std::numeric_limits<T>::min();

    BOOST_FOREACH(
        const MULTI_INDEX_TYPE other_multi_index,
        Multi_Index_Box_Intersect(MULTI_INDEX_CUBE< DIMENSION, -1, +1 >(multi_index), multi_index_box)
    ) {
        const int other_linear_index = multi_index_box.Linear_Index(other_multi_index);
        const T other_phi = phi_of_index(other_linear_index);
        if(phi * other_phi >= 0)
            continue;
        assert(multi_index != other_multi_index);
        const MULTI_INDEX_TYPE multi_offset = other_multi_index - multi_index;
        const T scale = std::sqrt(static_cast<T>(multi_offset.L1_Norm()));
        assert(scale >= 1);
        const T x = scale * phi / (phi - other_phi);
        assert(0 <= x && x <= scale);
        if(x < min_dist_to_vertex)
            phi = -(one_p_eps * other_phi * min_dist_to_vertex) /
                   (one_m_eps * (scale - min_dist_to_vertex));
    }
}

//#####################################################################
//#####################################################################

template< class T, class T_MULTI_INDEX_BOX >
inline void
Shift_Level_Set_Away_From_Vertices(
    const T_MULTI_INDEX_BOX& multi_index_box,
    const ARRAY_VIEW<T>& phi_of_index,
    const float min_dist_to_vertex,
    const int sign_of_zero = -1)
{
    const int n = multi_index_box.Size();
    assert(phi_of_index.Size() == n);
    assert(0 <= min_dist_to_vertex && min_dist_to_vertex <= 0.5f);
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    for(int linear_index = 1; linear_index <= n; ++linear_index)
        Shift_Level_Set_Value_Away_From_Vertices(
            multi_index_box,
            phi_of_index,
            min_dist_to_vertex,
            linear_index,
            sign_of_zero
        );
}

//#####################################################################
//#####################################################################

namespace Detail_Shift_Level_Set_Away_From_Vertices
{

template< class T, class T_MULTI_INDEX_BOX >
struct SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HELPER;

} // namespace Detail_Shift_Level_Set_Away_From_Vertices

template< class T, class T_MULTI_INDEX_BOX >
inline void
Shift_Level_Set_Away_From_Vertices_MT(
    const unsigned int n_thread,
    const T_MULTI_INDEX_BOX& multi_index_box,
    ARRAY_VIEW<T> phi_of_index,
    const float min_dist_to_vertex,
    const int sign_of_zero = -1)
{
    // TODO: We actually need to stripe the grid domain for this to be 100%
    // correct...
    typedef Detail_Shift_Level_Set_Away_From_Vertices::
        SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HELPER<
            T, T_MULTI_INDEX_BOX
        > SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HELPER_;
    assert(0 <= min_dist_to_vertex && min_dist_to_vertex <= 0.5f);
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    For_Each_MT(
        n_thread,
        1, multi_index_box.Size(),
        SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HELPER_(
            multi_index_box,
            phi_of_index,
            min_dist_to_vertex,
            sign_of_zero
        )
    );
}

//#####################################################################
//#####################################################################

namespace Detail_Shift_Level_Set_Away_From_Vertices
{

template< class T, class T_MULTI_INDEX_BOX >
struct SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HELPER,
        (( typename T_MULTI_INDEX_BOX const, multi_index_box ))
        (( typename ARRAY_VIEW<T>&, phi_of_index ))
        (( /******/ float const, min_dist_to_vertex ))
        (( /******/ int const, sign_of_zero ))
    )
public:
    typedef void result_type;
    void operator()(const int linear_index) const
    {
        Shift_Level_Set_Value_Away_From_Vertices(
            multi_index_box,
            phi_of_index,
            min_dist_to_vertex,
            linear_index,
            sign_of_zero
        );
    }
};

} // namespace Detail_Shift_Level_Set_Away_From_Vertices

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_SHIFT_LEVEL_SET_AWAY_FROM_VERTICES_HPP
