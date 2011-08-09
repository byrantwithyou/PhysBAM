//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_REGULAR_SUBSYS_BASE_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_REGULAR_SUBSYS_BASE_HPP

#include <cassert>

#include <Jeffrey_Utilities/Algorithm/Fill.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Coarsen_Sign_Of_Cell_Index.h"
#include "DOMAIN_SUBSYS_BASE.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_DERIVED, class T, int D, class T_STENCIL >
class DOMAIN_REGULAR_SUBSYS_BASE
    : public DOMAIN_SUBSYS_BASE< T_DERIVED, T, D >
{
    typedef DOMAIN_SUBSYS_BASE< T_DERIVED, T, D > DOMAIN_SUBSYS_BASE_;
public:
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef T_STENCIL STENCIL_TYPE;
    typedef signed char SIGN_TYPE;

    VECTOR<T,D> dx;
    ARRAY< SIGN_TYPE > sign_of_cell_index;
    ARRAY< STENCIL_TYPE > stencil_of_index;

    using DOMAIN_SUBSYS_BASE_::multi_index_bound;

    DOMAIN_REGULAR_SUBSYS_BASE(
        const MULTI_INDEX_BOUND<D>& multi_index_bound_,
        const VECTOR<T,D>& dx_);

private:
    class BETA_OF_CELL_INDEX_FUNCTION;
public:
    typedef BETA_OF_CELL_INDEX_FUNCTION BETA_OF_CELL_INDEX_FUNCTION_TYPE;
    BETA_OF_CELL_INDEX_FUNCTION_TYPE Beta_Of_Cell_Index_Function() const;

    void Resize(const MULTI_INDEX_BOUND<D>& new_multi_index_bound);
    void Zero_Stencils();
    void Zero_Stencils_MT(const unsigned int n_thread);

    template< class T_FINE >
    void Coarsen(const T_FINE& fine);

    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const;
    T Apply(const int linear_index, const ARRAY_VIEW<const T> x, const MULTI_INDEX_TYPE& strides) const;
    void Apply_Transpose(const int linear_index, ARRAY_VIEW<T> y, const T x, const MULTI_INDEX_TYPE& strides) const;

    using DOMAIN_SUBSYS_BASE_::Cell_Multi_Index_Bound;
    using DOMAIN_SUBSYS_BASE_::Apply;
    using DOMAIN_SUBSYS_BASE_::Apply_Transpose;

    // [Defined in T_DERIVED]
    // T Beta_Of_Cell_Index(const int cell_linear_index) const;
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
    // template< class T_FINE >
    // void Coarsen_Stencils(const T_FINE& fine);
    // typedef ... /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE
    // typedef ... CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
    // /***/ LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index);
    // CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE Linear_Index_Stencil_Proxy(const int linear_index) const;
    // typedef ... /***/ MULTI_INDEX_STENCIL_PROXY_TYPE
    // typedef ... CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
    // /***/ MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index);
    // CONST_MULTI_INDEX_STENCIL_PROXY_TYPE Multi_Index_Stencil_Proxy(const int linear_index) const;

protected:
    using DOMAIN_SUBSYS_BASE_::Derived;
};

//#####################################################################
//#####################################################################

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
DOMAIN_REGULAR_SUBSYS_BASE(
    const MULTI_INDEX_BOUND<D>& multi_index_bound_,
    const VECTOR<T,D>& dx_)
    : DOMAIN_SUBSYS_BASE_(multi_index_bound_),
      dx(dx_),
      sign_of_cell_index((multi_index_bound_ - 1).Size(), false), // uninit'ed
      stencil_of_index(multi_index_bound_.Size(), false) // uninit'ed
{ }

template< class T_DERIVED, class T, int D, class T_STENCIL >
class DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
BETA_OF_CELL_INDEX_FUNCTION
{
    friend class DOMAIN_REGULAR_SUBSYS_BASE;
    const T_DERIVED& m_this;
    explicit BETA_OF_CELL_INDEX_FUNCTION(const T_DERIVED& this_)
        : m_this(this_)
    { }
public:
    typedef T result_type;
    T operator()(const int cell_linear_index) const
    { return m_this.Beta_Of_Cell_Index(cell_linear_index); }
    T operator()(const MULTI_INDEX_TYPE& cell_multi_index) const
    { return operator()(m_this.Cell_Multi_Index_Bound().Linear_Index(cell_multi_index)); }
};

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline typename DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::BETA_OF_CELL_INDEX_FUNCTION_TYPE
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Beta_Of_Cell_Index_Function() const
{ return BETA_OF_CELL_INDEX_FUNCTION_TYPE(Derived()); }

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline void
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Resize(const MULTI_INDEX_BOUND<D>& new_multi_index_bound)
{
    multi_index_bound = new_multi_index_bound;
    sign_of_cell_index.Resize(Cell_Multi_Index_Bound().Size(), false, false); // uninit'ed
    stencil_of_index.Resize(multi_index_bound.Size(), false, false); // uninit'ed
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline void
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Zero_Stencils()
{ Fill(stencil_of_index, STENCIL_TYPE::Construct_Zero()); }

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline void
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Zero_Stencils_MT(const unsigned int n_thread)
{ Fill_MT(n_thread, stencil_of_index, STENCIL_TYPE::Construct_Zero()); }

template< class T_DERIVED, class T, int D, class T_STENCIL >
template< class T_FINE >
void
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Coarsen(const T_FINE& fine)
{
    dx = 2 * fine.dx;
    Resize((fine.multi_index_bound + 1) / 2);
    Coarsen_Sign_Of_Cell_Index(
        fine.Cell_Multi_Index_Bound(),
        fine.sign_of_cell_index,
        sign_of_cell_index
    );
    //Derived().Coarsen_Stencils(fine);
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline void
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y, const MULTI_INDEX_TYPE strides) const
{
    assert(x.Size() >= multi_index_bound.Size());
    assert(y.Size() >= multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    for(int linear_index = 1; linear_index <= multi_index_bound.Size(); ++linear_index)
        y(linear_index) += Apply(linear_index, x, strides);
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline T
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Apply(const int linear_index, const ARRAY_VIEW<const T> x, const MULTI_INDEX_TYPE& strides) const
{
    assert(x.Size() >= multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    return Derived().Apply_Stencil(stencil_of_index(linear_index), linear_index, x, strides);
}

template< class T_DERIVED, class T, int D, class T_STENCIL >
inline void
DOMAIN_REGULAR_SUBSYS_BASE< T_DERIVED, T, D, T_STENCIL >::
Apply_Transpose(const int linear_index, ARRAY_VIEW<T> y, const T x, const MULTI_INDEX_TYPE& strides) const
{
    assert(y.Size() >= multi_index_bound.Size());
    assert(strides == multi_index_bound.Strides());
    Derived().Apply_Transpose_Stencil(stencil_of_index(linear_index), linear_index, y, x, strides);
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_REGULAR_SUBSYS_BASE_HPP
