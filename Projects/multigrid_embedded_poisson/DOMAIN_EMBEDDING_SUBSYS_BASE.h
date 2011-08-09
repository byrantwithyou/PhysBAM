//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_EMBEDDING_SUBSYS_BASE_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_EMBEDDING_SUBSYS_BASE_HPP

#include <cassert>

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Coarsen_Embedding.h"
#include "DOMAIN_SUBSYS_BASE.h"

//#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL.h>
//#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL_PROXY.h>
//#include <Jeffrey_Utilities/Stencils/GEOMETRIC_PROLONGATION_STENCIL_PROXY_FUNCTION.h>
//#include <Jeffrey_Utilities/Stencils/GEOMETRIC_RESTRICTION_STENCIL_PROXY.h>
//#include <Jeffrey_Utilities/Stencils/Multiply_Stencils.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_DERIVED, class T, int D, class T_STENCIL >
class DOMAIN_EMBEDDING_SUBSYS_BASE
    : public DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >
{
    typedef DOMAIN_SUBSYS_BASE< T_DERIVED, T, D > DOMAIN_SUBSYS_BASE_;
public:
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef T_STENCIL STENCIL_TYPE;

    HASHTABLE<int,int> stencil_index_of_linear_index;
    ARRAY<int> linear_index_of_stencil_index;
    ARRAY< STENCIL_TYPE > stencils;

    using DOMAIN_SUBSYS_BASE_::multi_index_bound;

    explicit DOMAIN_EMBEDDING_SUBSYS_BASE(const MULTI_INDEX_BOUND<D>& multi_index_bound_);

    static bool Valid_Index(const int linear_index);

    void Init_Stencil_Index_Of_Linear_Index();

    template< class T_FINE >
    void Coarsen(const T_FINE& fine);

    T Diag(const int linear_index) const;
    int Stencil_N_Nonzero(const int linear_index) const;
    T Stencil_Sum(const int linear_index) const;

    template< class T_RHS_OF_INDEX >
    void Set_Dirichlet_Grid_BC(
        const int linear_index,
        const T p,
        const T_RHS_OF_INDEX& rhs_of_index);

    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const;
    T Apply(const int linear_index, const ARRAY_VIEW<const T> x, const MULTI_INDEX_TYPE& strides) const;
    void Apply_Transpose(const int linear_index, ARRAY_VIEW<T> y, const T x, const MULTI_INDEX_TYPE& strides) const;

    using DOMAIN_SUBSYS_BASE_::Apply;
    using DOMAIN_SUBSYS_BASE_::Apply_Transpose;
    using DOMAIN_SUBSYS_BASE_::Set_Dirichlet_Grid_BC;

    // [Defined in T_DERIVED]
    // T Diag(const STENCIL_TYPE& stencil, const int linear_index) const;
    // int Stencil_N_Nonzero(const STENCIL_TYPE& stencil, const int linear_index) const;
    // T Stencil_Sum(const STENCIL_TYPE& stencil, const int linear_index) const;
    // void Set_Dirichlet_Grid_BC(
    //     const STENCIL_TYPE& stencil,
    //     const int linear_index,
    //     const T p,
    //     T_RHS_OF_INDEX rhs_of_index);
    // T Apply_Stencil(
    //     const STENCIL_TYPE& stencil,
    //     const int linear_index,
    //     const ARRAY_VIEW<const T> x,
    //     const MULTI_INDEX_TYPE& strides) const;
    // void Apply_Transpose_Stencil(
    //     const STENCIL_TYPE& stencil,
    //     const int linear_index,
    //     ARRAY_VIEW<T> y,
    //     const T x,
    //     const MULTI_INDEX_TYPE& strides) const;
    // typedef ... /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE
    // typedef ... CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
    // /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index);
    // CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index) const;
    // /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index);
    // CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const;
    // typedef ... /***/ MULTI_INDEX_STENCIL_PROXY_TYPE
    // typedef ... CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
    // /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index);
    // CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index) const;
    // /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index);
    // CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const;

protected:
    using DOMAIN_SUBSYS_BASE_::Derived;
};

//#####################################################################
//#####################################################################

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
DOMAIN_EMBEDDING_SUBSYS_BASE(const MULTI_INDEX_BOUND<D>& multi_index_bound_)
    : DOMAIN_SUBSYS_BASE_(multi_index_bound_)
{ }

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline bool
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Valid_Index(const int /*linear_index*/)
{ return true; }

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline void
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Init_Stencil_Index_Of_Linear_Index()
{
    const int n_stencil = linear_index_of_stencil_index.Size();
    stencil_index_of_linear_index.Initialize_New_Table(n_stencil);
    for(int stencil_index = 1; stencil_index <= n_stencil; ++stencil_index) {
        const int linear_index = linear_index_of_stencil_index(stencil_index);
        stencil_index_of_linear_index.Insert(linear_index, stencil_index);
    }
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
template< class T_FINE >
void
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Coarsen(const T_FINE& fine)
{
    assert(0 == stencil_index_of_linear_index.Size());
    assert(0 == linear_index_of_stencil_index.Size());
    assert(0 == stencils.Size());

    const int n_fine = fine.stencils.Size();
    assert(n_fine == fine.stencil_index_of_linear_index.Size());
    assert(n_fine == fine.linear_index_of_stencil_index.Size());
    assert(n_fine == fine.stencils.Size());

    multi_index_bound = (fine.multi_index_bound + 1)/2;

    Coarsen_Embedding(
        fine.multi_index_bound,
        fine.stencil_index_of_linear_index,
        stencil_index_of_linear_index,
        linear_index_of_stencil_index
    );

    const int n_coarse = linear_index_of_stencil_index.Size();
    stencils.Exact_Resize(n_coarse, false); // uninit'ed

    // TODO: This needs some work...
    throw 0;
#if 0
    static const int STENCIL_WIDTH = 7;
    typedef CUBE_STENCIL< T, D, -(STENCIL_WIDTH - 1)/2, +(STENCIL_WIDTH - 1)/2 > LOCAL_STENCIL_TYPE;
    typedef CUBE_STENCIL_PROXY<  LOCAL_STENCIL_TYPE > LOCAL_STENCIL_PROXY_TYPE;

    for(int stencil_index = 1; stencil_index <= n_coarse; ++stencil_index) {
        const int linear_index = linear_index_of_stencil_index(stencil_index);
        const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
        LOCAL_STENCIL_TYPE local_stencil = LOCAL_STENCIL_TYPE::Construct_Zero();
        LOCAL_STENCIL_PROXY_TYPE local_stencil_proxy(multi_index, local_stencil);
        Multiply_Stencils(
            GEOMETRIC_RESTRICTION_STENCIL_PROXY<T,D>(multi_index),
            fine.Multi_Index_Stencil_Proxy_Function(),
            GEOMETRIC_PROLONGATION_STENCIL_PROXY_FUNCTION<T,D>(fine.multi_index_bound),
            local_stencil_proxy
        );
        Derived().Multi_Index_Stencil_Proxy(stencils(stencil_index), linear_index) = local_stencil_proxy;
    }
#endif
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline T
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Diag(const int linear_index) const
{
    const int* const p_stencil_index = stencil_index_of_linear_index.Get_Pointer(linear_index);
    return p_stencil_index ?
           Derived().Diag(stencils(*p_stencil_index), linear_index) :
           static_cast<T>(0);
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline int
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Stencil_N_Nonzero(const int linear_index) const
{
    const int* const p_stencil_index = stencil_index_of_linear_index.Get_Pointer(linear_index);
    return p_stencil_index ?
           Derived().Stencil_N_Nonzero(stencils(*p_stencil_index), linear_index) :
           0;
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline T
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Stencil_Sum(const int linear_index) const
{
    const int* const p_stencil_index = stencil_index_of_linear_index.Get_Pointer(linear_index);
    return p_stencil_index ?
           Derived().Stencil_Sum(stencils(*p_stencil_index), linear_index) :
           static_cast<T>(0);
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
template< class T_RHS_OF_INDEX >
inline void
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Set_Dirichlet_Grid_BC(
    const int linear_index,
    const T p,
    const T_RHS_OF_INDEX& rhs_of_index)
{
    const int* const p_stencil_index = stencil_index_of_linear_index.Get_Pointer(linear_index);
    if(p_stencil_index)
        Derived().Set_Dirichlet_Grid_BC(stencils(*p_stencil_index), linear_index, p, rhs_of_index);
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline void
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const
{
    assert(x.Size() == multi_index_bound.Size());
    assert(y.Size() == multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index) {
        const int linear_index = linear_index_of_stencil_index(stencil_index);
        y(linear_index) += Derived().Apply_Stencil(stencils(stencil_index), linear_index, x, strides);
    }
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline T
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Apply(const int linear_index, const ARRAY_VIEW<const T> x, const MULTI_INDEX_TYPE& strides) const
{
    assert(x.Size() == multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    const int* const p_stencil_index = stencil_index_of_linear_index.Get_Pointer(linear_index);
    return p_stencil_index ?
           Derived().Apply_Stencil(stencils(*p_stencil_index), linear_index, x, strides) :
           static_cast<T>(0);
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline void
DOMAIN_EMBEDDING_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Apply_Transpose(const int linear_index, ARRAY_VIEW<T> y, const T x, const MULTI_INDEX_TYPE& strides) const
{
    assert(y.Size() == multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    const int* const p_stencil_index = stencil_index_of_linear_index.Get_Pointer(linear_index);
    if(p_stencil_index)
        Derived().Apply_Transpose_Stencil(stencils(*p_stencil_index), linear_index, y, x, strides);
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_EMBEDDING_SUBSYS_BASE_HPP
