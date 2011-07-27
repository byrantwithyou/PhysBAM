//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Stencils/INDEX_TRANSFORM_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL_PROXY.h>

#include "DOMAIN_EMBEDDING_SUBSYS_BASE.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D >
class DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS
    : public DOMAIN_EMBEDDING_SUBSYS_BASE<
          DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>,
          T, D,
          UNSTRUCTURED_STENCIL<int,T>
      >
{
    typedef DOMAIN_EMBEDDING_SUBSYS_BASE<
        DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>,
        T, D,
        UNSTRUCTURED_STENCIL<int,T>
    > DOMAIN_EMBEDDING_SUBSYS_BASE_;
public:
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef typename DOMAIN_EMBEDDING_SUBSYS_BASE_::STENCIL_TYPE STENCIL_TYPE;

    using DOMAIN_EMBEDDING_SUBSYS_BASE_::multi_index_bound;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::stencil_index_of_linear_index;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::linear_index_of_stencil_index;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::stencils;

    explicit DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS(const MULTI_INDEX_BOUND<D>& multi_index_bound_);

    void Zero_Stencils();

    using DOMAIN_EMBEDDING_SUBSYS_BASE_::Diag;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::Stencil_N_Nonzero;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::Stencil_Sum;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::Set_Dirichlet_Grid_BC;

    typedef UNSTRUCTURED_STENCIL_PROXY< /***/ STENCIL_TYPE > /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE;
    typedef UNSTRUCTURED_STENCIL_PROXY< const STENCIL_TYPE > CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE;
    /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index);
    CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index) const;
    /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index);
    CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const;

    typedef INDEX_TRANSFORM_STENCIL_PROXY<
        /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE,
        MULTI_INDEX_BOUND<D>
    > /***/ MULTI_INDEX_STENCIL_PROXY_TYPE;
    typedef INDEX_TRANSFORM_STENCIL_PROXY<
        CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE,
        MULTI_INDEX_BOUND<D>
    > CONST_MULTI_INDEX_STENCIL_PROXY_TYPE;
    /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index);
    CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index) const;
    /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index);
    CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const;

private:
    friend class DOMAIN_EMBEDDING_SUBSYS_BASE<
        DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>,
        T, D,
        UNSTRUCTURED_STENCIL<int,T>
    >;
    static T Diag(const STENCIL_TYPE& stencil, const int linear_index);
    static int Stencil_N_Nonzero(const STENCIL_TYPE& stencil, const int linear_index);
    static T Stencil_Sum(const STENCIL_TYPE& stencil, const int linear_index);
    template< class T_RHS_OF_INDEX >
    void Set_Dirichlet_Grid_BC(
        STENCIL_TYPE& stencil,
        const int linear_index,
        const T p,
        T_RHS_OF_INDEX rhs_of_index);
    static T Apply_Stencil(
        const STENCIL_TYPE& stencil,
        const int linear_index,
        const ARRAY_VIEW<const T> x,
        const MULTI_INDEX_TYPE& strides);
    static void Apply_Transpose_Stencil(
        const STENCIL_TYPE& stencil,
        const int linear_index,
        ARRAY_VIEW<T> y,
        const T x,
        const MULTI_INDEX_TYPE& strides);
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS(const MULTI_INDEX_BOUND<D>& multi_index_bound_)
    : DOMAIN_EMBEDDING_SUBSYS_BASE_(multi_index_bound_)
{ }

template< class T, int D >
inline void
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Zero_Stencils()
{
    for(int stencil_index = 1; stencil_index <= stencils.Size(); ++stencil_index)
        stencils(stencil_index).values.Clean_Memory();
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy(const int linear_index)
{ return Linear_Index_Stencil_Proxy_Of_Stencil_Index(stencil_index_of_linear_index.Get(linear_index)); }

template< class T, int D >
inline typename DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy(const int linear_index) const
{
    static const STENCIL_TYPE dummy_stencil = STENCIL_TYPE();
    const int* const p_stencil_index = stencil_index_of_linear_index.Get_Pointer(linear_index);
    return CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE(p_stencil_index ? stencils(*p_stencil_index) : dummy_stencil);
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index)
{ return LINEAR_INDEX_STENCIL_PROXY_TYPE(stencils(stencil_index)); }

template< class T, int D >
inline typename DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const
{ return CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE(stencils(stencil_index)); }

template< class T, int D >
inline typename DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy(const int linear_index)
{ return MULTI_INDEX_STENCIL_PROXY_TYPE(Linear_Index_Stencil_Proxy(linear_index), multi_index_bound); }

template< class T, int D >
inline typename DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy(const int linear_index) const
{ return CONST_MULTI_INDEX_STENCIL_PROXY_TYPE(Linear_Index_Stencil_Proxy(linear_index), multi_index_bound); }

template< class T, int D >
inline typename DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index)
{
    return MULTI_INDEX_STENCIL_PROXY_TYPE(
        Linear_Index_Stencil_Proxy_Of_Stencil_Index(stencil_index),
        multi_index_bound
    );
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const
{
    return CONST_MULTI_INDEX_STENCIL_PROXY_TYPE(
        Linear_Index_Stencil_Proxy_Of_Stencil_Index(stencil_index),
        multi_index_bound
    );
}

template< class T, int D >
inline T
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Diag(const STENCIL_TYPE& stencil, const int linear_index)
{ return stencil(linear_index, static_cast<T>(0)); }

template< class T, int D >
inline int
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Stencil_N_Nonzero(const STENCIL_TYPE& stencil, const int /*linear_index*/)
{ return stencil.values.Size(); }

template< class T, int D >
inline T
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Stencil_Sum(const STENCIL_TYPE& stencil, const int /*linear_index*/)
{ return stencil.Sum(); }

template< class T, int D >
template< class T_RHS_OF_INDEX >
inline void
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Set_Dirichlet_Grid_BC(
    STENCIL_TYPE& stencil,
    const int linear_index,
    const T p,
    T_RHS_OF_INDEX rhs_of_index)
{
    BOOST_FOREACH( typename STENCIL_TYPE::INDEX_VALUE_TYPE const index_value, stencil ) {
        if(index_value.value == 0)
            continue;
        const int other_linear_index = index_value.index;
        if(other_linear_index == linear_index)
            continue;
        const int other_stencil_index = stencil_index_of_linear_index.Get(other_linear_index);
        STENCIL_TYPE& other_stencil = stencils(other_stencil_index);
        T& c = other_stencil(linear_index);
        assert(c == index_value.value);
        rhs_of_index(other_linear_index) -= c * p;
        c = 0;
    }
    stencil.values.Clean_Memory();
}

template< class T, int D >
inline T
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Apply_Stencil(
    const STENCIL_TYPE& stencil,
    const int /*linear_index*/,
    const ARRAY_VIEW<const T> x,
    const MULTI_INDEX_TYPE& /*strides*/)
{ return stencil.Apply(x); }

template< class T, int D >
inline void
DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D>::
Apply_Transpose_Stencil(
    const STENCIL_TYPE& stencil,
    const int /*linear_index*/,
    ARRAY_VIEW<T> y,
    const T x,
    const MULTI_INDEX_TYPE& /*strides*/)
{ stencil.Apply_Transpose(y, x); }

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS_HPP
