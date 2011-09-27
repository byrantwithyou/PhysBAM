//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_KRYLOV_VECTOR_WRAPPER_ARRAY_VIEW_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_KRYLOV_VECTOR_WRAPPER_ARRAY_VIEW_HPP

#include <boost/cast.hpp>

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T >
struct KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> >
    : KRYLOV_VECTOR_BASE<T>
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        KRYLOV_VECTOR_WRAPPER,
        (( typename ARRAY_VIEW<T>, v ))
    )

    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& y) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& y) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator*=(const T c) PHYSBAM_OVERRIDE;
    void Copy(const T c, const KRYLOV_VECTOR_BASE<T>& y) PHYSBAM_OVERRIDE;
    void Copy(const T c, const KRYLOV_VECTOR_BASE<T>& y, const KRYLOV_VECTOR_BASE<T>& z) PHYSBAM_OVERRIDE;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
};

//#####################################################################
//#####################################################################

template< class T >
inline KRYLOV_VECTOR_BASE<T>&
KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> >::
operator+=(const KRYLOV_VECTOR_BASE<T>& y)
{
    v += boost::polymorphic_cast< const KRYLOV_VECTOR_WRAPPER* >(&y)->v;
    return *this;
}

template< class T >
inline KRYLOV_VECTOR_BASE<T>&
KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> >::
operator-=(const KRYLOV_VECTOR_BASE<T>& y)
{
    v -= boost::polymorphic_cast< const KRYLOV_VECTOR_WRAPPER* >(&y)->v;
    return *this;
}

template< class T >
inline KRYLOV_VECTOR_BASE<T>&
KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> >::
operator*=(const T c)
{
    v *= c;
    return *this;
}

template< class T >
inline void
KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> >::
Copy(const T c, const KRYLOV_VECTOR_BASE<T>& y)
{ ARRAY_VIEW<T>::Copy(c, boost::polymorphic_cast< const KRYLOV_VECTOR_WRAPPER* >(&y)->v, v); }

template< class T >
inline void
KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> >::
Copy(const T c, const KRYLOV_VECTOR_BASE<T>& y, const KRYLOV_VECTOR_BASE<T>& z)
{
    ARRAY_VIEW<T>::Copy(
        c,
        boost::polymorphic_cast< const KRYLOV_VECTOR_WRAPPER* >(&y)->v,
        boost::polymorphic_cast< const KRYLOV_VECTOR_WRAPPER* >(&z)->v,
        v
    );
}

template< class T >
inline int
KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> >::
Raw_Size() const
{ return Scalar_Size(v); }

template< class T >
inline T&
KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> >::
Raw_Get(int i)
{ return Scalar_At1(v, i); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_KRYLOV_VECTOR_WRAPPER_ARRAY_VIEW_HPP
