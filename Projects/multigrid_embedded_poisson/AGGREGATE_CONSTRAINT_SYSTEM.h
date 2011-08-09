//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_AGGREGATE_CONSTRAINT_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_AGGREGATE_CONSTRAINT_SYSTEM_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/BOUNDED_LIST.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
struct AGGREGATE_CONSTRAINT_SYSTEM
{
    typedef T SCALAR_TYPE;

    typedef UNSTRUCTURED_STENCIL<int,T> STENCIL_TYPE;

    HASHTABLE<int,int> indy_index_of_index;
    ARRAY<int> index_of_indy_index;
    HASHTABLE< int, BOUNDED_LIST< int, (1 << D) > > stencils_containing_index;
    ARRAY<T> value_of_indy_index;
    ARRAY< STENCIL_TYPE > indyless_stencils;

    void Init_Indy_Index_Of_Index();
    void Init_Stencils_Containing_Index();

    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;
    void Apply_Transpose(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;
    T Apply(const int indy_index, const ARRAY_VIEW<const T> x) const;
    void Apply_Transpose(const int indy_index, ARRAY_VIEW<T> y, const T x) const;

    void Apply_Z(ARRAY_VIEW<T> x) const;
    void Apply_Z_Transpose(ARRAY_VIEW<T> x) const;

    typedef UNSTRUCTURED_STENCIL_PROXY< /***/ STENCIL_TYPE > /***/ INDYLESS_STENCIL_PROXY_TYPE;
    typedef UNSTRUCTURED_STENCIL_PROXY< const STENCIL_TYPE > CONST_INDYLESS_STENCIL_PROXY_TYPE;
    /***/ INDYLESS_STENCIL_PROXY_TYPE Indyless_Stencil_Proxy(const int indy_index);
    CONST_INDYLESS_STENCIL_PROXY_TYPE Indyless_Stencil_Proxy(const int indy_index) const;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline void
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Init_Indy_Index_Of_Index()
{
    const int n_indy = index_of_indy_index.Size();
    indy_index_of_index.Initialize_New_Table(n_indy);
    for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
        const int index = index_of_indy_index(indy_index);
        indy_index_of_index.Insert(index, indy_index);
    }
}

template< class T, int D >
inline void
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Init_Stencils_Containing_Index()
{
    for(int indy_index = 1; indy_index <= indyless_stencils.Size(); ++indy_index)
        BOOST_FOREACH( typename STENCIL_TYPE::INDEX_VALUE_TYPE const index_value, indyless_stencils(indy_index) )
            stencils_containing_index.Get_Or_Insert(index_value.index).Append(indy_index);
}

template< class T, int D >
inline void
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{
    assert(y.Size() == indyless_stencils.Size());
    for(int indy_index = 1; indy_index <= indyless_stencils.Size(); ++indy_index)
        y(indy_index) += Apply(indy_index, x);
}

template< class T, int D >
inline void
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Apply_Transpose(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{
    assert(x.Size() == indyless_stencils.Size());
    for(int indy_index = 1; indy_index <= indyless_stencils.Size(); ++indy_index)
        Apply_Transpose(indy_index, y, x(indy_index));
}

template< class T, int D >
inline T
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Apply(const int indy_index, const ARRAY_VIEW<const T> x) const
{
    const int index = index_of_indy_index(indy_index);
    return value_of_indy_index(indy_index) * x(index) + indyless_stencils(indy_index).Apply(x);
}

template< class T, int D >
inline void
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Apply_Transpose(const int indy_index, ARRAY_VIEW<T> y, const T x) const
{
    const int index = index_of_indy_index(indy_index);
    y(index) += value_of_indy_index(indy_index) * x;
    indyless_stencils(indy_index).Apply_Tranpose(y, x);
}

template< class T, int D >
inline void
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Apply_Z(ARRAY_VIEW<T> x) const
{
    for(int indy_index = 1; indy_index <= indyless_stencils.Size(); ++indy_index) {
        const int index = index_of_indy_index(indy_index);
        const T indy_value = value_of_indy_index(indy_index);
        x(index) = -indyless_stencils(indy_index).Apply(x) / indy_value;
    }
}

template< class T, int D >
inline void
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Apply_Z_Transpose(ARRAY_VIEW<T> x) const
{
    for(int indy_index = 1; indy_index <= indyless_stencils.Size(); ++indy_index) {
        const int index = index_of_indy_index(indy_index);
        const T indy_value = value_of_indy_index(indy_index);
        indyless_stencils(indy_index).Apply_Transpose(x, -x(index) / indy_value);
        x(index) = 0;
    }
}

template< class T, int D >
inline typename AGGREGATE_CONSTRAINT_SYSTEM<T,D>::INDYLESS_STENCIL_PROXY_TYPE
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Indyless_Stencil_Proxy(const int indy_index)
{ return INDYLESS_STENCIL_PROXY_TYPE(indyless_stencils(indy_index)); }

template< class T, int D >
inline typename AGGREGATE_CONSTRAINT_SYSTEM<T,D>::CONST_INDYLESS_STENCIL_PROXY_TYPE
AGGREGATE_CONSTRAINT_SYSTEM<T,D>::
Indyless_Stencil_Proxy(const int indy_index) const
{ return CONST_INDYLESS_STENCIL_PROXY_TYPE(indyless_stencils(indy_index)); }

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_AGGREGATE_CONSTRAINT_SYSTEM_HPP
