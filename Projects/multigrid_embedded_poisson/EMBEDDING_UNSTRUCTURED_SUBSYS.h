//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EMBEDDING_UNSTRUCTURED_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EMBEDDING_UNSTRUCTURED_SUBSYS_HPP

#include <boost/concept/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL_PROXY.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T >
struct EMBEDDING_UNSTRUCTURED_SUBSYS
{
    typedef T SCALAR_TYPE;

    typedef UNSTRUCTURED_STENCIL<int,T> STENCIL_TYPE;

    HASHTABLE<int,int> stencil_index_of_index;
    ARRAY<int> index_of_stencil_index;
    ARRAY< STENCIL_TYPE > stencils;

    static bool Valid_Index(const int index);

    void Zero_Stencils();
    void Init_Stencil_Index_Of_Index();

    template< class T_STENCIL_PROXY >
    void Add_Stencil_To(
        const int index,
        T_STENCIL_PROXY stencil_proxy) const;

    T Diag(const int index) const;
    int Stencil_N_Nonzero(const int index) const;
    T Stencil_Sum(const int index) const;

    template< class T_RHS_OF_INDEX >
    void Set_Dirichlet_Grid_BC(
        const int index,
        const T p,
        T_RHS_OF_INDEX rhs_of_index);

    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;
    T Apply(const int index, const ARRAY_VIEW<const T> x) const;
    void Apply_Transpose(const int index, ARRAY_VIEW<T> y, const T x) const;

    typedef UNSTRUCTURED_STENCIL_PROXY< /***/ STENCIL_TYPE > /***/ STENCIL_PROXY_TYPE;
    typedef UNSTRUCTURED_STENCIL_PROXY< const STENCIL_TYPE > CONST_STENCIL_PROXY_TYPE;
    /***/ STENCIL_PROXY_TYPE Stencil_Proxy(const int index);
    CONST_STENCIL_PROXY_TYPE Stencil_Proxy(const int index) const;
    /***/ STENCIL_PROXY_TYPE Stencil_Proxy_Of_Stencil_Index(const int stencil_index);
    CONST_STENCIL_PROXY_TYPE Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const;
};

//#####################################################################
//#####################################################################

template< class T >
inline bool
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Valid_Index(const int /*index*/)
{ return true; }

template< class T >
inline void
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Zero_Stencils()
{
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index)
        stencils(stencil_index).values.Clean_Memory();
}

template< class T >
inline void
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Init_Stencil_Index_Of_Index()
{
    const int n_stencil = index_of_stencil_index.Size();
    stencil_index_of_index.Initialize_New_Table(n_stencil);
    for(int stencil_index = 1; stencil_index <= n_stencil; ++stencil_index) {
        const int index = index_of_stencil_index(stencil_index);
        stencil_index_of_index.Insert(index, stencil_index);
    }
}

template< class T >
template< class T_STENCIL_PROXY >
inline void
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Add_Stencil_To(
    const int index,
    T_STENCIL_PROXY stencil_proxy) const
{
    BOOST_CONCEPT_ASSERT((STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, int >));
    stencil_proxy += Stencil_Proxy(index);
}

template< class T >
inline T
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Diag(const int index) const
{
    const int* const p_stencil_index = stencil_index_of_index.Get_Pointer(index);
    return p_stencil_index ? stencils(*p_stencil_index)(index, static_cast<T>(0)) : static_cast<T>(0);
}

template< class T >
inline int
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Stencil_N_Nonzero(const int index) const
{
    const int* const p_stencil_index = stencil_index_of_index.Get_Pointer(index);
    return p_stencil_index ? stencils(*p_stencil_index).values.Size() : 0;
}

template< class T >
inline T
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Stencil_Sum(const int index) const
{
    const int* const p_stencil_index = stencil_index_of_index.Get_Pointer(index);
    return p_stencil_index ? stencils(*p_stencil_index).Sum() : static_cast<T>(0);
}

template< class T >
template< class T_RHS_OF_INDEX >
inline void
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Set_Dirichlet_Grid_BC(
    const int index,
    const T p,
    T_RHS_OF_INDEX rhs_of_index)
{
    const int* const p_stencil_index = stencil_index_of_index.Get_Pointer(index);
    if(!p_stencil_index)
        return;
    STENCIL_TYPE& stencil = stencils(*p_stencil_index);
    BOOST_FOREACH( typename STENCIL_TYPE::INDEX_VALUE_TYPE const index_value, stencil ) {
        if(index_value.value == 0)
            continue;
        const int other_index = index_value.index;
        if(other_index == index)
            continue;
        const int other_stencil_index = stencil_index_of_index.Get(other_index);
        STENCIL_TYPE& other_stencil = stencils(other_stencil_index);
        T& c = other_stencil(index);
        assert(c == index_value.value);
        rhs_of_index(other_index) -= c * p;
        c = 0;
    }
    stencil.values.Clean_Memory();
}

template< class T >
inline void
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index) {
        const int index = index_of_stencil_index(stencil_index);
        y(index) += stencils(stencil_index).Apply(x);
    }
}

template< class T >
inline T
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Apply(const int index, const ARRAY_VIEW<const T> x) const
{
    const int* const p_stencil_index = stencil_index_of_index.Get_Pointer(index);
    return p_stencil_index ? stencils(*p_stencil_index).Apply(x) : static_cast<T>(0);
}

template< class T >
inline void
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Apply_Transpose(const int index, ARRAY_VIEW<T> y, const T x) const
{
    const int* const p_stencil_index = stencil_index_of_index.Get_Pointer(index);
    if(p_stencil_index)
        stencils(*p_stencil_index).Apply_Transpose(y, x);
}

template< class T >
inline typename EMBEDDING_UNSTRUCTURED_SUBSYS<T>::STENCIL_PROXY_TYPE
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Stencil_Proxy(const int index)
{ return Stencil_Proxy_Of_Stencil_Index(stencil_index_of_index.Get(index)); }

template< class T >
inline typename EMBEDDING_UNSTRUCTURED_SUBSYS<T>::CONST_STENCIL_PROXY_TYPE
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Stencil_Proxy(const int index) const
{
    static const STENCIL_TYPE dummy_stencil = STENCIL_TYPE();
    const int* const p_stencil_index = stencil_index_of_index.Get_Pointer(index);
    return CONST_STENCIL_PROXY_TYPE(p_stencil_index ? stencils(*p_stencil_index) : dummy_stencil);
}

template< class T >
inline typename EMBEDDING_UNSTRUCTURED_SUBSYS<T>::STENCIL_PROXY_TYPE
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Stencil_Proxy_Of_Stencil_Index(const int stencil_index)
{ return STENCIL_PROXY_TYPE(stencils(stencil_index)); }

template< class T >
inline typename EMBEDDING_UNSTRUCTURED_SUBSYS<T>::CONST_STENCIL_PROXY_TYPE
EMBEDDING_UNSTRUCTURED_SUBSYS<T>::
Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const
{ return CONST_STENCIL_PROXY_TYPE(stencils(stencil_index)); }

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EMBEDDING_UNSTRUCTURED_SUBSYS_HPP
