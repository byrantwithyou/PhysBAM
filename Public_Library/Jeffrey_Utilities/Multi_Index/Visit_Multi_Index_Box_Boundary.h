//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_VISIT_MULTI_INDEX_BOX_BOUNDARY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_VISIT_MULTI_INDEX_BOX_BOUNDARY_HPP

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T_MULTI_INDEX_BOX, class T_VISITOR >
inline void
Visit_Multi_Index_Box_Boundary(
    const T_MULTI_INDEX_BOX& multi_index_box,
    T_VISITOR visitor)
{
    static const int DIMENSION = T_MULTI_INDEX_BOX::DIMENSION;
    typedef VECTOR< int, DIMENSION > MULTI_INDEX_TYPE;
    typedef VECTOR< int, DIMENSION - 1 > SUB_MULTI_INDEX_TYPE;

    const VECTOR< MULTI_INDEX_TYPE, 2 > minmax_multi_index(
        multi_index_box.Min_Multi_Index(),
        multi_index_box.Max_Multi_Index()
    );
    const MULTI_INDEX_TYPE& min_multi_index = minmax_multi_index[1];
    const MULTI_INDEX_TYPE& max_multi_index = minmax_multi_index[2];

    for(int d = 1; d <= DIMENSION; ++d) {
        const MULTI_INDEX_BOX< DIMENSION - 1 > sub_multi_index_box(
            min_multi_index.Remove_Index(d),
            max_multi_index.Remove_Index(d)
        );
        BOOST_FOREACH( const SUB_MULTI_INDEX_TYPE sub_multi_index, sub_multi_index_box ) {
            for(int e = 1; e < d; ++e)
                if(sub_multi_index[e] == min_multi_index[e] || sub_multi_index[e] == max_multi_index[e])
                    goto CONTINUE_FOREACH_SUB_MULTI_INDEX;
            for(int i = 1; i <= 2; ++i) {
                const MULTI_INDEX_TYPE multi_index = sub_multi_index.Insert(minmax_multi_index[i][d], d);
                visitor(multi_index);
            }
            CONTINUE_FOREACH_SUB_MULTI_INDEX:;
        }
    }
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_VISIT_MULTI_INDEX_BOX_BOUNDARY_HPP
