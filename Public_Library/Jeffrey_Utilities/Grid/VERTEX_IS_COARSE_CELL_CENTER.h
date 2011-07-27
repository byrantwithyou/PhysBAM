//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_IS_COARSE_CELL_CENTER_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_IS_COARSE_CELL_CENTER_HPP

#include <boost/mpl/assert.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

namespace PhysBAM
{

template< int COARSE_FACTOR, int D >
struct VERTEX_IS_COARSE_CELL_CENTER
{
    BOOST_MPL_ASSERT_RELATION( (COARSE_FACTOR & 1), ==, 0 );
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VERTEX_IS_COARSE_CELL_CENTER, (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
    )
public:
    typedef bool result_type;
    bool operator()(const int linear_index) const
    { return operator()(multi_index_bound.Multi_Index(linear_index)); }
    bool operator()(const VECTOR<int,D>& multi_index) const
    {
        for(int d = 1; d <= D; ++d)
            if((multi_index[d] - 1) % COARSE_FACTOR != COARSE_FACTOR / 2)
                return false;
        return true;
    }
};

template< int COARSE_FACTOR, int D >
inline VERTEX_IS_COARSE_CELL_CENTER< COARSE_FACTOR, D >
Make_Vertex_Is_Coarse_Cell_Center(const MULTI_INDEX_BOUND<D>& multi_index_bound)
{ return VERTEX_IS_COARSE_CELL_CENTER< COARSE_FACTOR, D >(multi_index_bound); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_IS_COARSE_CELL_CENTER_HPP
