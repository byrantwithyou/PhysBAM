//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS_HPP

#include <cassert>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CROSS_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/Beta_Dv_Over_Dx_Dx.h>
#include <Jeffrey_Utilities/Stencils/CROSS_CONSTBETA_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/CROSS_CONSTBETA_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/INDEX_TRANSFORM_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/SKIP_ZERO_VALUE_STENCIL_PROXY.h>

#include "DOMAIN_REGULAR_SUBSYS_BASE.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
class DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS
    : public DOMAIN_REGULAR_SUBSYS_BASE<
          DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>,
          T, D,
          CROSS_CONSTBETA_STENCIL<D>
      >
{
    typedef DOMAIN_REGULAR_SUBSYS_BASE<
        DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>,
        T, D,
        CROSS_CONSTBETA_STENCIL<D>
    > DOMAIN_REGULAR_SUBSYS_BASE_;
public:
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef typename DOMAIN_REGULAR_SUBSYS_BASE_::STENCIL_TYPE STENCIL_TYPE;

    T beta;
    VECTOR<T,D> beta_dv_over_dx_dx;

    using DOMAIN_REGULAR_SUBSYS_BASE_::multi_index_bound;
    using DOMAIN_REGULAR_SUBSYS_BASE_::dx;
    using DOMAIN_REGULAR_SUBSYS_BASE_::stencil_of_index;

    DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS(
        const MULTI_INDEX_BOUND<D>& multi_index_bound_,
        const VECTOR<T,D>& dx_,
        const T beta_);

    void Init_Beta_Dv_Over_Dx_Dx();

    T Beta_Of_Cell_Index(const int cell_linear_index) const;
    T Diag(const int linear_index) const;
    int Stencil_N_Nonzero(const int linear_index) const;
    T Stencil_Sum(const int linear_index) const;

    void Zero_Stencil(const int linear_index);
    template< class T_RHS_OF_INDEX >
    void Set_Dirichlet_Grid_BC(
        const int linear_index,
        const T p,
        T_RHS_OF_INDEX rhs_of_index);

    T Apply_Stencil(
        const STENCIL_TYPE& stencil,
        const int linear_index,
        const ARRAY_VIEW<const T> x,
        const MULTI_INDEX_TYPE& strides) const;
    void Apply_Transpose_Stencil(
        const STENCIL_TYPE& stencil,
        const int linear_index,
        ARRAY_VIEW<T> y,
        const T x,
        const MULTI_INDEX_TYPE& strides) const;

    typedef CROSS_CONSTBETA_STENCIL_PROXY< /***/ STENCIL_TYPE, T > /***/ MULTI_INDEX_STENCIL_PROXY_TYPE;
    typedef CROSS_CONSTBETA_STENCIL_PROXY< const STENCIL_TYPE, T > CONST_MULTI_INDEX_STENCIL_PROXY_TYPE;
    /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index);
    CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index) const;

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

    using DOMAIN_REGULAR_SUBSYS_BASE_::Cell_Multi_Index_Bound;
    using DOMAIN_REGULAR_SUBSYS_BASE_::Resize;
    using DOMAIN_REGULAR_SUBSYS_BASE_::Zero_Stencils;
    using DOMAIN_REGULAR_SUBSYS_BASE_::Zero_Stencils_MT;
    using DOMAIN_REGULAR_SUBSYS_BASE_::Set_Dirichlet_Grid_BC;

private:
    friend class DOMAIN_REGULAR_SUBSYS_BASE<
        DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>,
        T, D,
        CROSS_CONSTBETA_STENCIL<D>
    >;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS(
    const MULTI_INDEX_BOUND<D>& multi_index_bound_,
    const VECTOR<T,D>& dx_,
    const T beta_)
    : DOMAIN_REGULAR_SUBSYS_BASE_(multi_index_bound_, dx_),
      beta(beta_),
      beta_dv_over_dx_dx(Beta_Dv_Over_Dx_Dx(beta, dx))
{ }

template< class T, int D >
inline void
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Init_Beta_Dv_Over_Dx_Dx()
{ beta_dv_over_dx_dx = Beta_Dv_Over_Dx_Dx(beta, dx); }

template< class T, int D >
inline T
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Beta_Of_Cell_Index(const int /*cell_linear_index*/) const
{ return beta; }

template< class T, int D >
inline T
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Diag(const int linear_index) const
{ return stencil_of_index(linear_index).Center(beta_dv_over_dx_dx); }

template< class T, int D >
inline int
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Stencil_N_Nonzero(const int linear_index) const
{ return stencil_of_index(linear_index).N_Nonzero(); }

template< class T, int D >
inline T
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Stencil_Sum(const int linear_index) const
{ return stencil_of_index(linear_index).Sum(beta_dv_over_dx_dx); }

template< class T, int D >
inline void
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Zero_Stencil(const int linear_index)
{
    STENCIL_TYPE& stencil = stencil_of_index(linear_index);
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    for(MULTI_INDEX_CROSS_ITERATOR<D> it(multi_index); it.Valid(); ++it) {
        const MULTI_INDEX_TYPE other_multi_index = *it;
        const MULTI_INDEX_TYPE multi_offset = other_multi_index - multi_index;
        if(stencil(multi_offset) == 0)
            continue;
        assert(multi_index_bound.Contains(other_multi_index));
        const int other_linear_index = multi_index_bound.Linear_Index(other_multi_index);
        if(other_linear_index == linear_index)
            continue;
        assert(multi_offset != MULTI_INDEX_TYPE());
        STENCIL_TYPE& other_stencil = stencil_of_index(other_linear_index);
        assert(other_stencil(-multi_offset) == stencil(multi_offset));
        other_stencil(-multi_offset) = static_cast< unsigned char >(0);
    }
    stencil.Zero();
}

template< class T, int D >
template< class T_RHS_OF_INDEX >
inline void
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Set_Dirichlet_Grid_BC(
    const int linear_index,
    const T p,
    T_RHS_OF_INDEX rhs_of_index)
{
    STENCIL_TYPE& stencil = stencil_of_index(linear_index);
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    for(MULTI_INDEX_CROSS_ITERATOR<D> it(multi_index); it.Valid(); ++it) {
        const MULTI_INDEX_TYPE other_multi_index = *it;
        const MULTI_INDEX_TYPE multi_offset = other_multi_index - multi_index;
        if(stencil(multi_offset) == 0)
            continue;
        assert(multi_index_bound.Contains(other_multi_index));
        const int other_linear_index = multi_index_bound.Linear_Index(other_multi_index);
        if(other_linear_index == linear_index)
            continue;
        assert(multi_offset != MULTI_INDEX_TYPE());
        STENCIL_TYPE& other_stencil = stencil_of_index(other_linear_index);
        assert(other_stencil(-multi_offset) == stencil(multi_offset));
        rhs_of_index(other_linear_index) -= other_stencil(-multi_offset, beta_dv_over_dx_dx) * p;
        other_stencil(-multi_offset) = static_cast< unsigned char >(0);
    }
    stencil.Zero();
    stencil.Center() = static_cast< unsigned char >(1);
    rhs_of_index(linear_index) = stencil.Center(beta_dv_over_dx_dx) * p;
}

template< class T, int D >
inline T
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Apply_Stencil(
    const STENCIL_TYPE& stencil,
    const int linear_index,
    const ARRAY_VIEW<const T> x,
    const MULTI_INDEX_TYPE& strides) const
{ return stencil.Apply(beta_dv_over_dx_dx, linear_index, strides, x); }

template< class T, int D >
inline void
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Apply_Transpose_Stencil(
    const STENCIL_TYPE& stencil,
    const int linear_index,
    ARRAY_VIEW<T> y,
    const T x,
    const MULTI_INDEX_TYPE& strides) const
{ return stencil.Apply_Transpose(beta_dv_over_dx_dx, linear_index, strides, y, x); }

template< class T, int D >
inline typename DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy(const int linear_index)
{
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    return MULTI_INDEX_STENCIL_PROXY_TYPE(beta_dv_over_dx_dx, multi_index, stencil_of_index(linear_index));
}

template< class T, int D >
inline typename DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Multi_Index_Stencil_Proxy(const int linear_index) const
{
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    return CONST_MULTI_INDEX_STENCIL_PROXY_TYPE(beta_dv_over_dx_dx, multi_index, stencil_of_index(linear_index));
}

template< class T, int D >
inline typename DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy(const int linear_index)
{
    return LINEAR_INDEX_STENCIL_PROXY_TYPE(
        Make_Skip_Zero_Value_Stencil_Proxy(Multi_Index_Stencil_Proxy(linear_index)),
        multi_index_bound
    );
}

template< class T, int D >
inline typename DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS<T,D>::
Linear_Index_Stencil_Proxy(const int linear_index) const
{
    return CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE(
        Make_Skip_Zero_Value_Stencil_Proxy(Multi_Index_Stencil_Proxy(linear_index)),
        multi_index_bound
    );
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_REGULAR_CROSS_CONSTBETA_SUBSYS_HPP
