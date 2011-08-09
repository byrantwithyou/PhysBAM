//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_EMBEDDING_CUBE_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_EMBEDDING_CUBE_SUBSYS_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/INDEX_TRANSFORM_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/SKIP_ZERO_VALUE_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>

#include "DOMAIN_EMBEDDING_SUBSYS_BASE.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
class DOMAIN_EMBEDDING_CUBE_SUBSYS
    : public DOMAIN_EMBEDDING_SUBSYS_BASE<
          DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>,
          T, D,
          CUBE_STENCIL<T,D,-1,+1>
      >
{
    typedef DOMAIN_EMBEDDING_SUBSYS_BASE<
        DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>,
        T, D,
        CUBE_STENCIL<T,D,-1,+1>
    > DOMAIN_EMBEDDING_SUBSYS_BASE_;
public:
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef typename DOMAIN_EMBEDDING_SUBSYS_BASE_::STENCIL_TYPE STENCIL_TYPE;

    using DOMAIN_EMBEDDING_SUBSYS_BASE_::multi_index_bound;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::stencil_index_of_linear_index;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::linear_index_of_stencil_index;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::stencils;

    explicit DOMAIN_EMBEDDING_CUBE_SUBSYS(const MULTI_INDEX_BOUND<D>& multi_index_bound_);

    void Zero_Stencils();

    using DOMAIN_EMBEDDING_SUBSYS_BASE_::Diag;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::Stencil_N_Nonzero;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::Stencil_Sum;
    using DOMAIN_EMBEDDING_SUBSYS_BASE_::Set_Dirichlet_Grid_BC;

    typedef CUBE_STENCIL_PROXY< /***/ STENCIL_TYPE > /***/ MULTI_INDEX_STENCIL_PROXY_TYPE;
    typedef CUBE_STENCIL_PROXY< const STENCIL_TYPE > CONST_MULTI_INDEX_STENCIL_PROXY_TYPE;
    /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index);
    CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index) const;
    /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index);
    CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const;

    typedef INDEX_TRANSFORM_STENCIL_PROXY<
        SKIP_ZERO_VALUE_STENCIL_PROXY< /***/ MULTI_INDEX_STENCIL_PROXY_TYPE >,
        MULTI_INDEX_BOUND<D>
    > /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE;
    typedef INDEX_TRANSFORM_STENCIL_PROXY<
        SKIP_ZERO_VALUE_STENCIL_PROXY< CONST_MULTI_INDEX_STENCIL_PROXY_TYPE >,
        MULTI_INDEX_BOUND<D>
    > CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE;
    /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index);
    CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index) const;
    /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index);
    CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const;

private:
    friend class DOMAIN_EMBEDDING_SUBSYS_BASE<
        DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>,
        T, D,
        CUBE_STENCIL<T,D,-1,+1>
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
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
DOMAIN_EMBEDDING_CUBE_SUBSYS(const MULTI_INDEX_BOUND<D>& multi_index_bound_)
    : DOMAIN_EMBEDDING_SUBSYS_BASE_(multi_index_bound_)
{ }

template< class T, int D >
inline void
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Zero_Stencils()
{ stencils.Fill(STENCIL_TYPE::Construct_Zero()); }

template< class T, int D >
inline typename DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy(const int linear_index)
{
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    const int stencil_index = stencil_index_of_linear_index.Get(linear_index);
    return MULTI_INDEX_STENCIL_PROXY_TYPE(multi_index, stencils(stencil_index));
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy(const int linear_index) const
{
    static const STENCIL_TYPE dummy_stencil = STENCIL_TYPE::Construct_Zero();
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    const int* const p_stencil_index = stencil_index_of_linear_index.Get_Pointer(linear_index);
    return CONST_MULTI_INDEX_STENCIL_PROXY_TYPE(
        multi_index,
        p_stencil_index ? stencils(*p_stencil_index) : dummy_stencil
    );
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index)
{
    const int linear_index = linear_index_of_stencil_index(stencil_index);
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    return MULTI_INDEX_STENCIL_PROXY_TYPE(multi_index, stencils(stencil_index));
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const
{
    const int linear_index = linear_index_of_stencil_index(stencil_index);
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    return CONST_MULTI_INDEX_STENCIL_PROXY_TYPE(multi_index, stencils(stencil_index));
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy(const int linear_index)
{
    return LINEAR_INDEX_STENCIL_PROXY_TYPE(
        Make_Skip_Zero_Value_Stencil_Proxy(Multi_Index_Stencil_Proxy(linear_index)),
        multi_index_bound
    );
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy(const int linear_index) const
{
    return CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE(
        Make_Skip_Zero_Value_Stencil_Proxy(Multi_Index_Stencil_Proxy(linear_index)),
        multi_index_bound
    );
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index)
{
    return LINEAR_INDEX_STENCIL_PROXY_TYPE(
        Make_Skip_Zero_Value_Stencil_Proxy(Multi_Index_Stencil_Proxy_Of_Stencil_Index(stencil_index)),
        multi_index_bound
    );
}

template< class T, int D >
inline typename DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy_Of_Stencil_Index(const int stencil_index) const
{
    return CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE(
        Make_Skip_Zero_Value_Stencil_Proxy(Multi_Index_Stencil_Proxy_Of_Stencil_Index(stencil_index)),
        multi_index_bound
    );
}

template< class T, int D >
inline T
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Diag(const STENCIL_TYPE& stencil, const int /*linear_index*/)
{ return stencil(MULTI_INDEX_TYPE()); }

template< class T, int D >
inline int
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Stencil_N_Nonzero(const STENCIL_TYPE& stencil, const int /*linear_index*/)
{ return stencil.N_Nonzero(); }

template< class T, int D >
inline T
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Stencil_Sum(const STENCIL_TYPE& stencil, const int /*linear_index*/)
{ return stencil.Sum(); }

template< class T, int D >
template< class T_RHS_OF_INDEX >
inline void
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Set_Dirichlet_Grid_BC(
    STENCIL_TYPE& stencil,
    const int linear_index,
    const T p,
    T_RHS_OF_INDEX rhs_of_index)
{
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    BOOST_FOREACH( const MULTI_INDEX_TYPE multi_offset, (STATIC_MULTI_INDEX_CUBE<D,-1,+1>()) ) {
        if(stencil(multi_offset) == 0)
            continue;
        const MULTI_INDEX_TYPE other_multi_index = multi_index + multi_offset;
        assert(multi_index_bound.Contains(other_multi_index));
        const int other_linear_index = multi_index_bound.Linear_Index(other_multi_index);
        if(other_linear_index == linear_index)
            continue;
        const int other_stencil_index = stencil_index_of_linear_index.Get(other_linear_index);
        STENCIL_TYPE& other_stencil = stencils(other_stencil_index);
        T& c = other_stencil(-multi_offset);
        assert(c == stencil(multi_offset));
        rhs_of_index(other_linear_index) -= c * p;
        c = 0;
    }
    stencil.Zero();
}

template< class T, int D >
inline T
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Apply_Stencil(
    const STENCIL_TYPE& stencil,
    const int linear_index,
    const ARRAY_VIEW<const T> x,
    const MULTI_INDEX_TYPE& strides)
{ return stencil.Apply(linear_index, strides, x); }

template< class T, int D >
inline void
DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D>::
Apply_Transpose_Stencil(
    const STENCIL_TYPE& stencil,
    const int linear_index,
    ARRAY_VIEW<T> y,
    const T x,
    const MULTI_INDEX_TYPE& strides)
{ stencil.Apply_Transpose(linear_index, strides, y, x); }

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_EMBEDDING_CUBE_SUBSYS_HPP
