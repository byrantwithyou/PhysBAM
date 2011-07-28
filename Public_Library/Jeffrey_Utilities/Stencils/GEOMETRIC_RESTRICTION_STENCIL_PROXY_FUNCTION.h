//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_FUNCTION_HPP

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Stencils/GEOMETRIC_RESTRICTION_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D >
struct GEOMETRIC_RESTRICTION_STENCIL_PROXY_FUNCTION
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    GEOMETRIC_RESTRICTION_STENCIL_PROXY_FUNCTION()
    { }
    explicit GEOMETRIC_RESTRICTION_STENCIL_PROXY_FUNCTION(const MULTI_INDEX_BOUND<D>& coarse_multi_index_bound)
        : m_coarse_multi_index_bound(coarse_multi_index_bound)
    { }

    typedef GEOMETRIC_RESTRICTION_STENCIL_PROXY<T,D> result_type;

    result_type operator()(const MULTI_INDEX_TYPE& coarse_base_multi_index) const
    { return result_type(coarse_base_multi_index); }
    result_type operator()(const int coarse_base_linear_index) const
    { return result_type(m_coarse_multi_index_bound.Multi_Index(coarse_base_linear_index)); }

private:
    const MULTI_INDEX_BOUND<D> m_coarse_multi_index_bound;
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_FUNCTION_HPP
