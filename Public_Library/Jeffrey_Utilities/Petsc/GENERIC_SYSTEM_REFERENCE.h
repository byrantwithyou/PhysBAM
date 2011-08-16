//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_GENERIC_SYSTEM_REFERENCE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_GENERIC_SYSTEM_REFERENCE_HPP

#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL_PROXY.h>

namespace PhysBAM
{

namespace Petsc
{

template< class T >
struct GENERIC_SYSTEM_REFERENCE
{
    typedef T SCALAR_TYPE;

    template< class T_SYSTEM >
    explicit GENERIC_SYSTEM_REFERENCE(const T_SYSTEM& system);

    T Diag(const int index) const;
    int Stencil_N_Nonzero(const int index) const;

    typedef UNSTRUCTURED_STENCIL_PROXY<
        UNSTRUCTURED_STENCIL<int,T>
    > ADD_STENCIL_TO_STENCIL_PROXY_TYPE;
    void Add_Stencil_To(
        const int index,
        const ADD_STENCIL_TO_STENCIL_PROXY_TYPE& stencil_proxy) const;

private:
    struct DISPATCHER;
    template< class T_SYSTEM > struct ACCESSOR;
    const void* const p_system;
    const DISPATCHER* const p_dispatcher;
};

//#####################################################################
//#####################################################################

template< class T >
template< class T_SYSTEM >
inline
GENERIC_SYSTEM_REFERENCE<T>::
GENERIC_SYSTEM_REFERENCE(const T_SYSTEM& system)
    : p_system(&system),
      p_dispatcher(ACCESSOR< T_SYSTEM >::Dispatcher())
{ }

template< class T >
inline T
GENERIC_SYSTEM_REFERENCE<T>::
Diag(const int index) const
{ return (*p_dispatcher->Diag)(p_system, index); }

template< class T >
inline int
GENERIC_SYSTEM_REFERENCE<T>::
Stencil_N_Nonzero(const int index) const
{ return (*p_dispatcher->Stencil_N_Nonzero)(p_system, index); }

template< class T >
inline void
GENERIC_SYSTEM_REFERENCE<T>::
Add_Stencil_To(const int index, const ADD_STENCIL_TO_STENCIL_PROXY_TYPE& stencil_proxy) const
{ (*p_dispatcher->Add_Stencil_To)(p_system, index, stencil_proxy); }

template< class T >
struct GENERIC_SYSTEM_REFERENCE<T>::DISPATCHER
{
    T (*Diag)(const void* const, const int);
    int (*Stencil_N_Nonzero)(const void* const, const int);
    void (*Add_Stencil_To)(const void* const, const int, const ADD_STENCIL_TO_STENCIL_PROXY_TYPE&);
};

template< class T >
template< class T_SYSTEM >
struct GENERIC_SYSTEM_REFERENCE<T>::ACCESSOR
{
    static const DISPATCHER* Dispatcher()
    {
        static const DISPATCHER dispatcher = {
            &Diag,
            &Stencil_N_Nonzero,
            &Add_Stencil_To
        };
        return &dispatcher;
    }

    static T Diag(const void* const p_system, const int index)
    { return static_cast< const T_SYSTEM* >(p_system)->Diag(index); }
    static int Stencil_N_Nonzero(const void* const p_system, const int index)
    { return static_cast< const T_SYSTEM* >(p_system)->Stencil_N_Nonzero(index); }
    static void Add_Stencil_To(const void* const p_system, const int index, const ADD_STENCIL_TO_STENCIL_PROXY_TYPE& stencil_proxy)
    { static_cast< const T_SYSTEM* >(p_system)->Add_Stencil_To(index, stencil_proxy); }
};

} // namespace Petsc

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_GENERIC_SYSTEM_REFERENCE_HPP
