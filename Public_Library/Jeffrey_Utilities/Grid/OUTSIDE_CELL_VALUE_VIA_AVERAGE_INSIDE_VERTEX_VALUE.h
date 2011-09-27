//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_OUTSIDE_CELL_VALUE_VIA_AVERAGE_INSIDE_VERTEX_VALUE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_OUTSIDE_CELL_VALUE_VIA_AVERAGE_INSIDE_VERTEX_VALUE_HPP

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>
#include <Jeffrey_Utilities/RESULT_OF.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>

namespace PhysBAM
{

template< class T_CELL_MULTI_INDEX_BOX, class T_VALUE_OF_INDEX >
struct OUTSIDE_CELL_VALUE_VIA_AVERAGE_INSIDE_VERTEX_VALUE
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        OUTSIDE_CELL_VALUE_VIA_AVERAGE_INSIDE_VERTEX_VALUE,
        (( typename T_CELL_MULTI_INDEX_BOX const, cell_multi_index_box ))
        (( typename T_VALUE_OF_INDEX const, value_of_index ))
    )
public:

    template<class> struct result;
    template< class T_THIS, class T_MULTI_INDEX >
    struct result< T_THIS ( T_MULTI_INDEX ) >
        : RESULT_OF< const T_VALUE_OF_INDEX (
              typename REMOVE_QUALIFIERS< T_MULTI_INDEX >::type
          ) >
    { };

    template< int D >
    typename RESULT_OF< const T_VALUE_OF_INDEX ( VECTOR<int,D> ) >::type
    operator()(const VECTOR<int,D> outside_cell_multi_index) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;
        typedef typename RESULT_OF< const T_VALUE_OF_INDEX ( MULTI_INDEX_TYPE ) >::type RESULT_TYPE;
        const MULTI_INDEX_TYPE clamped_cell_multi_index = cell_multi_index_box.Clamp(outside_cell_multi_index);
        assert((outside_cell_multi_index - clamped_cell_multi_index).Max_Abs() == 1);
        const int count = 1 << Count(outside_cell_multi_index - clamped_cell_multi_index, 0);
        RESULT_TYPE sum = RESULT_TYPE();
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE multi_index,
            Multi_Index_Box_Intersect(
                MULTI_INDEX_CUBE<D,0,1>(outside_cell_multi_index),
                MULTI_INDEX_CUBE<D,0,1>(clamped_cell_multi_index)
            )
        )
            sum += value_of_index(multi_index);
        return sum / count;
    }
};

template< class T_CELL_MULTI_INDEX_BOX, class T_VALUE_OF_INDEX >
inline OUTSIDE_CELL_VALUE_VIA_AVERAGE_INSIDE_VERTEX_VALUE<
    T_CELL_MULTI_INDEX_BOX, T_VALUE_OF_INDEX
>
Make_Outside_Cell_Value_Via_Average_Inside_Vertex_Value(
    const T_CELL_MULTI_INDEX_BOX& cell_multi_index_box,
    const T_VALUE_OF_INDEX& value_of_index)
{
    return OUTSIDE_CELL_VALUE_VIA_AVERAGE_INSIDE_VERTEX_VALUE<
        T_CELL_MULTI_INDEX_BOX, T_VALUE_OF_INDEX
    >(cell_multi_index_box, value_of_index);
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_OUTSIDE_CELL_VALUE_VIA_AVERAGE_INSIDE_VERTEX_VALUE_HPP
