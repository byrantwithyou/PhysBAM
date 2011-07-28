//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_DIRICHLET_CONSTRAINT_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_DIRICHLET_CONSTRAINT_SYSTEM_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D >
struct DIRICHLET_CONSTRAINT_SYSTEM
{
    typedef T SCALAR_TYPE;
    static const int DIMENSION = D;

    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    typedef CUBE_STENCIL<T,D,0,1> STENCIL_TYPE;

    MULTI_INDEX_BOUND<D> multi_index_bound;
    HASHTABLE<int,int> stencil_index_of_cell_linear_index;
    ARRAY<int> cell_linear_index_of_stencil_index;
    ARRAY< STENCIL_TYPE > stencils;

    explicit DIRICHLET_CONSTRAINT_SYSTEM(const MULTI_INDEX_BOUND<D>& multi_index_bound_);

    MULTI_INDEX_BOUND<D> Cell_Multi_Index_Bound() const;

    void Init_Stencil_Index_Of_Cell_Linear_Index();

    template< class T_RHS_OF_STENCIL_INDEX >
    void Set_Dirichlet_Grid_BC(
        const int linear_index,
        const T p,
        const T_RHS_OF_STENCIL_INDEX& rhs_of_stencil_index);
    template< class T_RHS_OF_STENCIL_INDEX >
    void Set_Dirichlet_Grid_BC(
        const MULTI_INDEX_TYPE& multi_index,
        const T p,
        T_RHS_OF_STENCIL_INDEX rhs_of_stencil_index);

    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;
    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const;
    void Apply_Transpose(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;
    void Apply_Transpose(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const;
    T Apply(const int stencil_index, const ARRAY_VIEW<const T> x) const;
    T Apply(const int stencil_index, const ARRAY_VIEW<const T> x, const MULTI_INDEX_TYPE& strides) const;
    void Apply_Transpose(const int stencil_index, ARRAY_VIEW<T> y, const T x) const;
    void Apply_Transpose(const int stencil_index, ARRAY_VIEW<T> y, const T x, const MULTI_INDEX_TYPE& strides) const;

    typedef CUBE_STENCIL_PROXY< /***/ STENCIL_TYPE > /***/ MULTI_INDEX_STENCIL_PROXY_TYPE;
    typedef CUBE_STENCIL_PROXY< const STENCIL_TYPE > CONST_MULTI_INDEX_STENCIL_PROXY_TYPE;
    /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int stencil_index);
    CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int stencil_index) const;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
DIRICHLET_CONSTRAINT_SYSTEM(const MULTI_INDEX_BOUND<D>& multi_index_bound_)
    : multi_index_bound(multi_index_bound_)
{ }

template< class T, int D >
inline MULTI_INDEX_BOUND<D>
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Cell_Multi_Index_Bound() const
{ return multi_index_bound - 1; }

template< class T, int D >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Init_Stencil_Index_Of_Cell_Linear_Index()
{
    const int n_stencil = cell_linear_index_of_stencil_index.Size();
    stencil_index_of_cell_linear_index.Initialize_New_Table(n_stencil);
    for(int stencil_index = 1; stencil_index <= n_stencil; ++stencil_index) {
        const int cell_linear_index = cell_linear_index_of_stencil_index(stencil_index);
        stencil_index_of_cell_linear_index.Insert(cell_linear_index, stencil_index);
    }
}

template< class T, int D >
template< class T_RHS_OF_STENCIL_INDEX >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Set_Dirichlet_Grid_BC(
    const int linear_index,
    const T p,
    const T_RHS_OF_STENCIL_INDEX& rhs_of_stencil_index)
{
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    Set_Dirichlet_Grid_BC(multi_index, p, rhs_of_stencil_index);
}

template< class T, int D >
template< class T_RHS_OF_STENCIL_INDEX >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Set_Dirichlet_Grid_BC(
    const MULTI_INDEX_TYPE& multi_index,
    const T p,
    T_RHS_OF_STENCIL_INDEX rhs_of_stencil_index)
{
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = Cell_Multi_Index_Bound();
    BOOST_FOREACH(
        const MULTI_INDEX_TYPE cell_multi_index,
        Multi_Index_Box_Intersect(MULTI_INDEX_CUBE<D,-1,0>(multi_index), cell_multi_index_bound)
    ) {
        const int cell_linear_index = cell_multi_index_bound.Linear_Index(cell_multi_index);
        const int* const p_stencil_index = stencil_index_of_cell_linear_index.Get_Pointer(cell_linear_index);
        if(!p_stencil_index)
            continue;
        const int stencil_index = *p_stencil_index;
        const MULTI_INDEX_TYPE multi_offset = multi_index - cell_multi_index;
        STENCIL_TYPE& stencil = stencils(stencil_index);
        T& c = stencil(multi_offset);
        rhs_of_stencil_index(stencil_index) -= c * p;
        c = 0;
    }
}

template< class T, int D >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{ Apply(x, y, multi_index_bound.Strides()); }

template< class T, int D >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const
{
    assert(x.Size() == multi_index_bound.Size());
    assert(y.Size() == stencils.Size());
    assert(strides == multi_index_bound.Strides());
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index)
        y(stencil_index) += Apply(stencil_index, x, strides);
}

template< class T, int D >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Apply_Transpose(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{ Apply_Transpose(x, y, multi_index_bound.Strides()); }

template< class T, int D >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Apply_Transpose(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const
{
    assert(x.Size() == stencils.Size());
    assert(y.Size() == multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index)
        Apply_Transpose(stencil_index, y, x(stencil_index), strides);
}

template< class T, int D >
inline T
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Apply(const int stencil_index, const ARRAY_VIEW<const T> x) const
{ return Apply(stencil_index, x, multi_index_bound.Strides()); }

template< class T, int D >
inline T
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Apply(const int stencil_index, const ARRAY_VIEW<const T> x, const MULTI_INDEX_TYPE& strides) const
{
    assert(x.Size() == multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    const int cell_linear_index = cell_linear_index_of_stencil_index(stencil_index);
    const MULTI_INDEX_TYPE cell_multi_index = Cell_Multi_Index_Bound().Multi_Index(cell_linear_index);
    const int base_linear_index = multi_index_bound.Linear_Index(cell_multi_index);
    return stencils(stencil_index).Apply(base_linear_index, strides, x);
}

template< class T, int D >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Apply_Transpose(const int stencil_index, ARRAY_VIEW<T> y, const T x) const
{ Apply_Transpose(stencil_index, y, x, multi_index_bound.Strides()); }

template< class T, int D >
inline void
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Apply_Transpose(const int stencil_index, ARRAY_VIEW<T> y, const T x, const MULTI_INDEX_TYPE& strides) const
{
    assert(y.Size() == multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    const int cell_linear_index = cell_linear_index_of_stencil_index(stencil_index);
    const MULTI_INDEX_TYPE cell_multi_index = Cell_Multi_Index_Bound().Multi_Index(cell_linear_index);
    const int base_linear_index = multi_index_bound.Linear_Index(cell_multi_index);
    stencils(stencil_index).Apply_Transpose(base_linear_index, strides, y, x);
}

template< class T, int D >
inline typename DIRICHLET_CONSTRAINT_SYSTEM<T,D>::MULTI_INDEX_STENCIL_PROXY_TYPE
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Multi_Index_Stencil_Proxy(const int stencil_index)
{
    const int cell_linear_index = cell_linear_index_of_stencil_index(stencil_index);
    const MULTI_INDEX_TYPE cell_multi_index = Cell_Multi_Index_Bound().Multi_Index(cell_linear_index);
    return MULTI_INDEX_STENCIL_PROXY_TYPE(cell_multi_index, stencils(stencil_index));
}

template< class T, int D >
inline typename DIRICHLET_CONSTRAINT_SYSTEM<T,D>::CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
DIRICHLET_CONSTRAINT_SYSTEM<T,D>::
Multi_Index_Stencil_Proxy(const int stencil_index) const
{
    const int cell_linear_index = cell_linear_index_of_stencil_index(stencil_index);
    const MULTI_INDEX_TYPE cell_multi_index = Cell_Multi_Index_Bound().Multi_Index(cell_linear_index);
    return CONST_MULTI_INDEX_STENCIL_PROXY_TYPE(cell_multi_index, stencils(stencil_index));
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_DIRICHLET_CONSTRAINT_SYSTEM_HPP
