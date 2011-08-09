//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INTERFACE_CONSTRAINT_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INTERFACE_CONSTRAINT_SYSTEM_HPP

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/BOUNDED_LIST.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
struct INTERFACE_CONSTRAINT_SYSTEM
{
    typedef T SCALAR_TYPE;
    static const int DIMENSION = D;

    typedef UNSTRUCTURED_STENCIL< int, T, 2 * (1 << D) > STENCIL_TYPE;

    HASHTABLE<int,int> stencil_index_of_cell_index;
    ARRAY<int> cell_index_of_stencil_index;
    HASHTABLE< int, BOUNDED_LIST< int, (1 << D) > > stencils_containing_index;
    ARRAY< STENCIL_TYPE > stencils;

    void Init_Stencil_Index_Of_Cell_Index();
    void Init_Stencils_Containing_Index();

    template< class T_RHS_OF_STENCIL_INDEX >
    void Set_Dirichlet_Grid_BC(
        const int index,
        const T p,
        T_RHS_OF_STENCIL_INDEX rhs_of_stencil_index);

    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;
    void Apply_Transpose(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;
    T Apply(const int stencil_index, const ARRAY_VIEW<const T> x) const;
    void Apply_Transpose(const int stencil_index, ARRAY_VIEW<T> y, const T x) const;

    typedef UNSTRUCTURED_STENCIL_PROXY< /***/ STENCIL_TYPE > /***/ STENCIL_PROXY_TYPE;
    typedef UNSTRUCTURED_STENCIL_PROXY< const STENCIL_TYPE > CONST_STENCIL_PROXY_TYPE;
    /***/ STENCIL_PROXY_TYPE Stencil_Proxy(const int stencil_index);
    CONST_STENCIL_PROXY_TYPE Stencil_Proxy(const int stencil_index) const;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline void
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Init_Stencil_Index_Of_Cell_Index()
{
    const int n_stencil = cell_index_of_stencil_index.Size();
    stencil_index_of_cell_index.Initialize_New_Table(n_stencil);
    for(int stencil_index = 1; stencil_index <= n_stencil; ++stencil_index) {
        const int cell_index = cell_index_of_stencil_index(stencil_index);
        stencil_index_of_cell_index.Insert(cell_index, stencil_index);
    }
}

template< class T, int D >
inline void
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Init_Stencils_Containing_Index()
{
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index)
        BOOST_FOREACH( typename STENCIL_TYPE::INDEX_VALUE_TYPE const index_value, stencils(stencil_index) )
            stencils_containing_index.Get_Or_Insert(index_value.index).Append(stencil_index);
}

template< class T, int D >
template< class T_RHS_OF_STENCIL_INDEX >
inline void
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Set_Dirichlet_Grid_BC(
    const int index,
    const T p,
    T_RHS_OF_STENCIL_INDEX rhs_of_stencil_index)
{
    const BOUNDED_LIST< int, (1 << D) >* const p_stencil_indices = stencils_containing_index.Get_Pointer(index);
    if(!p_stencil_indices)
        return;
    BOOST_FOREACH( const int stencil_index, *p_stencil_indices ) {
        STENCIL_TYPE& stencil = stencils(stencil_index);
        T& c = stencil(index);
        rhs_of_stencil_index(stencil_index) -= c*p;
        c = 0;
    }
}

template< class T, int D >
inline void
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index)
        y(stencil_index) += Apply(stencil_index, x);
}

template< class T, int D >
inline void
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Apply_Transpose(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index)
        Apply_Transpose(stencil_index, y, x(stencil_index));
}

template< class T, int D >
inline T
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Apply(const int stencil_index, const ARRAY_VIEW<const T> x) const
{ return stencils(stencil_index).Apply(x); }

template< class T, int D >
inline void
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Apply_Transpose(const int stencil_index, ARRAY_VIEW<T> y, const T x) const
{ stencils(stencil_index).Apply_Transpose(y, x); }

template< class T, int D >
inline typename INTERFACE_CONSTRAINT_SYSTEM<T,D>::STENCIL_PROXY_TYPE
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Stencil_Proxy(const int stencil_index)
{ return STENCIL_PROXY_TYPE(stencils(stencil_index)); }

template< class T, int D >
inline typename INTERFACE_CONSTRAINT_SYSTEM<T,D>::CONST_STENCIL_PROXY_TYPE
INTERFACE_CONSTRAINT_SYSTEM<T,D>::
Stencil_Proxy(const int stencil_index) const
{ return CONST_STENCIL_PROXY_TYPE(stencils(stencil_index)); }

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INTERFACE_CONSTRAINT_SYSTEM_HPP
