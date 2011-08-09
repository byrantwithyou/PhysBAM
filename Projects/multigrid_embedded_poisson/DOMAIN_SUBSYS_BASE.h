//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// TODO: Should probably clean up all the *_Index_Stencil_Proxy_Function's...
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_SUBSYS_BASE_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_SUBSYS_BASE_HPP

#include <boost/concept/assert.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/PROPAGATE_CONST.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_DERIVED, class T, int D >
struct DOMAIN_SUBSYS_BASE
{
    typedef T SCALAR_TYPE;
    static const int DIMENSION = D;
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        DOMAIN_SUBSYS_BASE,
        (( typename MULTI_INDEX_BOUND<D>, multi_index_bound ))
    )

    MULTI_INDEX_BOUND<D> Cell_Multi_Index_Bound() const;

    bool Valid_Index(const int linear_index) const;

    template< class T_STENCIL_PROXY >
    typename boost::enable_if<
        boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, int >
    >::type Add_Stencil_To(
        const int linear_index,
        T_STENCIL_PROXY stencil_proxy) const;
    template< class T_STENCIL_PROXY >
    typename boost::enable_if<
        boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, VECTOR<int,D> >
    >::type Add_Stencil_To(
        const int linear_index,
        T_STENCIL_PROXY stencil_proxy) const;

    template< class T_RHS_OF_INDEX >
    void Set_Dirichlet_Grid_BC(
        const MULTI_INDEX_TYPE& multi_index,
        const T p,
        const T_RHS_OF_INDEX& rhs_of_index);

    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;
    T Apply(const int linear_index, const ARRAY_VIEW<const T> x) const;
    void Apply_Transpose(const int linear_index, ARRAY_VIEW<T> y, const T x) const;

    // [Defined in T_DERIVED]
    // void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const;
    // T Apply(const int linear_index, const ARRAY_VIEW<const T> x, const MULTI_INDEX_TYPE strides) const;
    // void Apply_Transpose(const int linear_index, ARRAY_VIEW<T> y, const T x, const MULTI_INDEX_TYPE strides) const;
    // typedef ... /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE
    // typedef ... CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
    // /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index);
    // CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index) const;
    // typedef ... /***/ MULTI_INDEX_STENCIL_PROXY_TYPE
    // typedef ... CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
    // /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index);
    // CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index) const;

protected:
    typedef T_DERIVED DERIVED_TYPE;
    /***/ DERIVED_TYPE& Derived() /***/ { return static_cast< /***/ DERIVED_TYPE& >(*this); }
    const DERIVED_TYPE& Derived() const { return static_cast< const DERIVED_TYPE& >(*this); }
};

//#####################################################################
//#####################################################################

template< class T_DERIVED, class T, int D >
inline MULTI_INDEX_BOUND<D>
DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >::
Cell_Multi_Index_Bound() const
{ return multi_index_bound - 1; }

template< class T_DERIVED, class T, int D >
inline bool
DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >::
Valid_Index(const int linear_index) const
{ return linear_index <= multi_index_bound.Size(); }

template< class T_DERIVED, class T, int D >
template< class T_STENCIL_PROXY >
inline typename boost::enable_if<
    boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, int >
>::type
DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >::
Add_Stencil_To(
    const int linear_index,
    T_STENCIL_PROXY stencil_proxy) const
{
    BOOST_CONCEPT_ASSERT((STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    stencil_proxy += Derived().Linear_Index_Stencil_Proxy(linear_index);
}

template< class T_DERIVED, class T, int D >
template< class T_STENCIL_PROXY >
inline typename boost::enable_if<
    boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, VECTOR<int,D> >
>::type
DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >::
Add_Stencil_To(
    const int linear_index,
    T_STENCIL_PROXY stencil_proxy) const
{
    BOOST_CONCEPT_ASSERT((STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    stencil_proxy += Derived().Multi_Index_Stencil_Proxy(linear_index);
}

template< class T_DERIVED, class T, int D >
template< class T_RHS_OF_INDEX >
inline void
DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >::
Set_Dirichlet_Grid_BC(
    const MULTI_INDEX_TYPE& multi_index,
    const T p,
    const T_RHS_OF_INDEX& rhs_of_index)
{
    const int linear_index = multi_index_bound.Linear_Index(multi_index);
    Derived().Set_Dirichlet_Grid_BC(linear_index, p, rhs_of_index);
}

template< class T_DERIVED, class T, int D >
inline void
DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{ Derived().Apply(x, y, multi_index_bound.Strides()); }

template< class T_DERIVED, class T, int D >
inline T
DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >::
Apply(const int linear_index, const ARRAY_VIEW<const T> x) const
{ return Derived().Apply(linear_index, x, multi_index_bound.Strides()); }

template< class T_DERIVED, class T, int D >
inline void
DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >::
Apply_Transpose(const int linear_index, ARRAY_VIEW<T> y, const T x) const
{ Derived().Apply_Transpose(linear_index, y, x, multi_index_bound.Strides()); }

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_SUBSYS_BASE_HPP
