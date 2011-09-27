//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GENERIC_SYSTEM_REFERENCE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GENERIC_SYSTEM_REFERENCE_HPP

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL_PROXY.h>

namespace PhysBAM
{

template< class T >
struct GENERIC_SYSTEM_REFERENCE
{
    typedef T SCALAR_TYPE;

    template< class T_SYSTEM >
    explicit GENERIC_SYSTEM_REFERENCE(const T_SYSTEM& system);

    T Diag(const int index) const;
    int Stencil_N_Nonzero(const int index) const;
    T Stencil_Sum(const int index) const;

    typedef UNSTRUCTURED_STENCIL_PROXY<
        UNSTRUCTURED_STENCIL<int,T>
    > ADD_STENCIL_TO_STENCIL_PROXY_TYPE;
    void Add_Stencil_To(
        const int index,
        const ADD_STENCIL_TO_STENCIL_PROXY_TYPE& stencil_proxy) const;

    void Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const;

private:
    struct DISPATCHER;
    template< class T_SYSTEM > struct ACCESSOR;
    const void* const mp_system;
    const DISPATCHER& m_dispatcher;
};

//#####################################################################
//#####################################################################

template< class T >
template< class T_SYSTEM >
inline
GENERIC_SYSTEM_REFERENCE<T>::
GENERIC_SYSTEM_REFERENCE(const T_SYSTEM& system)
    : mp_system(&system),
      m_dispatcher(ACCESSOR< T_SYSTEM >::Dispatcher())
{ }

template< class T >
inline T
GENERIC_SYSTEM_REFERENCE<T>::
Diag(const int index) const
{ return (*m_dispatcher.Diag)(mp_system, index); }

template< class T >
inline int
GENERIC_SYSTEM_REFERENCE<T>::
Stencil_N_Nonzero(const int index) const
{ return (*m_dispatcher.Stencil_N_Nonzero)(mp_system, index); }

template< class T >
inline T
GENERIC_SYSTEM_REFERENCE<T>::
Stencil_Sum(const int index) const
{ return (*m_dispatcher.Stencil_Sum)(mp_system, index); }

template< class T >
inline void
GENERIC_SYSTEM_REFERENCE<T>::
Add_Stencil_To(const int index, const ADD_STENCIL_TO_STENCIL_PROXY_TYPE& stencil_proxy) const
{ (*m_dispatcher.Add_Stencil_To)(mp_system, index, stencil_proxy); }

template< class T >
inline void
GENERIC_SYSTEM_REFERENCE<T>::
Apply(const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y) const
{ (*m_dispatcher.Apply)(mp_system, x, y); }

template< class T >
struct GENERIC_SYSTEM_REFERENCE<T>::DISPATCHER
{
    T (*Diag)(const void* const, const int);
    int (*Stencil_N_Nonzero)(const void* const, const int);
    T (*Stencil_Sum)(const void* const, const int);
    void (*Add_Stencil_To)(const void* const, const int, const ADD_STENCIL_TO_STENCIL_PROXY_TYPE&);
    void (*Apply)(const void* const, const ARRAY_VIEW<const T>, ARRAY_VIEW<T>);
};

template< class T >
template< class T_SYSTEM >
struct GENERIC_SYSTEM_REFERENCE<T>::ACCESSOR
{
    static const DISPATCHER& Dispatcher()
    {
        static const DISPATCHER dispatcher = {
            &Diag,
            &Stencil_N_Nonzero,
            &Stencil_Sum,
            &Add_Stencil_To,
            &Apply
        };
        return dispatcher;
    }

    static T Diag(const void* const p_system, const int index)
    { return static_cast< const T_SYSTEM* >(p_system)->Diag(index); }
    static int Stencil_N_Nonzero(const void* const p_system, const int index)
    { return static_cast< const T_SYSTEM* >(p_system)->Stencil_N_Nonzero(index); }
    static T Stencil_Sum(const void* const p_system, const int index)
    { return static_cast< const T_SYSTEM* >(p_system)->Stencil_Sum(index); }
    static void Add_Stencil_To(const void* const p_system, const int index, const ADD_STENCIL_TO_STENCIL_PROXY_TYPE& stencil_proxy)
    { static_cast< const T_SYSTEM* >(p_system)->Add_Stencil_To(index, stencil_proxy); }
    static void Apply(const void* const p_system, const ARRAY_VIEW<const T> x, ARRAY_VIEW<T> y)
    { static_cast< const T_SYSTEM* >(p_system)->Apply(x, y); }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GENERIC_SYSTEM_REFERENCE_HPP
