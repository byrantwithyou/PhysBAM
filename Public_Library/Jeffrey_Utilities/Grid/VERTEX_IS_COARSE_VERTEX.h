//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_IS_COARSE_VERTEX_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_IS_COARSE_VERTEX_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

namespace PhysBAM
{

template< int COARSE_FACTOR, int D >
struct VERTEX_IS_COARSE_VERTEX
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VERTEX_IS_COARSE_VERTEX,
        (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
    )
public:
    typedef bool result_type;
    bool operator()(const int linear_index) const
    { return operator()(multi_index_bound.Multi_Index(linear_index)); }
    bool operator()(const VECTOR<int,D>& multi_index) const
    {
        for(int d = 1; d <= D; ++d)
            if((multi_index[d] - 1) % COARSE_FACTOR != 0)
                return false;
        return true;
    }
};

template< int COARSE_FACTOR, int D >
inline VERTEX_IS_COARSE_VERTEX< COARSE_FACTOR, D >
Make_Vertex_Is_Coarse_Vertex(const MULTI_INDEX_BOUND<D>& multi_index_bound)
{ return VERTEX_IS_COARSE_VERTEX< COARSE_FACTOR, D >(multi_index_bound); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_IS_COARSE_VERTEX_HPP
