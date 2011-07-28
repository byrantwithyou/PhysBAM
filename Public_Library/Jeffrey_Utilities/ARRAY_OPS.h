//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ARRAY_OPS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ARRAY_OPS_HPP

#include <vector>

#include <Jeffrey_Utilities/PROPAGATE_CONST.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM
{

//#####################################################################
// struct ARRAY_VALUE< T_ARRAY >
//#####################################################################

template< class T_ARRAY >
struct ARRAY_VALUE
{ typedef typename T_ARRAY::value_type type; };

template< class T_ARRAY >
struct ARRAY_VALUE< const T_ARRAY >
    : ARRAY_VALUE< T_ARRAY >
{ };

//#####################################################################
// Front(T_ARRAY& a) -> ...
//#####################################################################

template< class T, class T_DERIVED, class T_ID >
inline T&
Front(ARRAY_BASE< T, T_DERIVED, T_ID >& a)
{ return static_cast< T_DERIVED& >(a)(1); }

template< class T, class T_DERIVED, class T_ID >
inline const T&
Front(const ARRAY_BASE< T, T_DERIVED, T_ID >& a)
{ return static_cast< const T_DERIVED& >(a)(1); }

template< class T, class T_ALLOCATOR >
inline typename std::vector< T, T_ALLOCATOR >::reference
Front(std::vector< T, T_ALLOCATOR >& v)
{ return v.front(); }

template< class T, class T_ALLOCATOR >
inline typename std::vector< T, T_ALLOCATOR >::const_reference
Front(const std::vector< T, T_ALLOCATOR >& v)
{ return v.front(); }

//#####################################################################
// Back(T_ARRAY& a) -> ...
//#####################################################################

template< class T, class T_DERIVED, class T_ID >
inline T&
Back(ARRAY_BASE< T, T_DERIVED, T_ID >& a)
{ return static_cast< T_DERIVED& >(a)(static_cast< const T_DERIVED& >(a).Size()); }

template< class T, class T_DERIVED, class T_ID >
inline const T&
Back(const ARRAY_BASE< T, T_DERIVED, T_ID >& a)
{ return static_cast< const T_DERIVED& >(a)(static_cast< const T_DERIVED& >(a).Size()); }

template< class T, class T_ALLOCATOR >
inline typename std::vector< T, T_ALLOCATOR >::reference
Back(std::vector< T, T_ALLOCATOR >& v)
{ return v.back(); }

template< class T, class T_ALLOCATOR >
inline typename std::vector< T, T_ALLOCATOR >::const_reference
Back(const std::vector< T, T_ALLOCATOR >& v)
{ return v.back(); }

//#####################################################################
// Size(const T_ARRAY& a) -> ...
//#####################################################################

template< class T, class T_DERIVED, class T_ID >
inline T_ID
Size(const ARRAY_BASE< T, T_DERIVED, T_ID >& a)
{ return static_cast< const T_DERIVED& >(a).Size(); }

template< class T, class T_ALLOCATOR >
inline typename std::vector< T, T_ALLOCATOR >::size_type
Size(const std::vector< T, T_ALLOCATOR >& v)
{ return v.size(); }

//#####################################################################
// Resize(T_ARRAY& a, ... new_size) -> void
//#####################################################################

template< class T, class T_DERIVED, class T_ID >
inline void
Resize(
    ARRAY_BASE< T, T_DERIVED, T_ID >& a, const T_ID new_size,
    const bool initialize_new_elements = true,
    const bool copy_existing_elements = true)
{ static_cast< T_DERIVED& >(a).Resize(new_size, initialize_new_elements, copy_existing_elements); }

template< class T, class T_ALLOCATOR >
inline void
Resize(
    std::vector< T, T_ALLOCATOR >& v,
    typename std::vector< T, T_ALLOCATOR >::size_type const new_size,
    const bool /*initialize_new_elements*/ = true,
    const bool /*copy_existing_elements*/ = true)
{ v.resize(new_size); }

//#####################################################################
// Exact_Resize(T_ARRAY& a, ... new_size) -> void
//#####################################################################

template< class T_ARRAY, class T_SIZE >
inline void
Exact_Resize(
    T_ARRAY& a, const T_SIZE new_size,
    const bool initialize_new_elements = true,
    const bool copy_existing_elements = true)
{ Resize(a, new_size, initialize_new_elements, copy_existing_elements); }

template< class T, class T_ID >
inline void
Exact_Resize(
    ARRAY<T,T_ID>& a, const T_ID new_size,
    const bool initialize_new_elements = true,
    const bool /*copy_existing_elements*/ = true)
{ a.Exact_Resize(new_size, initialize_new_elements); }

//#####################################################################
// As_Array_View(T_ARRAY& a) -> ARRAY_VIEW< ... >
// As_Const_Array_View(const T_ARRAY& a) -> ARRAY_VIEW< ... >
//#####################################################################

template< class T_ARRAY >
inline ARRAY_VIEW<
    typename PROPAGATE_CONST< T_ARRAY, typename T_ARRAY::ELEMENT >::type,
    typename T_ARRAY::INDEX
>
As_Array_View(T_ARRAY& a)
{
    return ARRAY_VIEW<
        typename PROPAGATE_CONST< T_ARRAY, typename T_ARRAY::ELEMENT >::type,
        typename T_ARRAY::INDEX
    >(a);
}

template< class T_ARRAY >
inline ARRAY_VIEW< typename T_ARRAY::ELEMENT const, typename T_ARRAY::INDEX >
As_Const_Array_View(const T_ARRAY& a)
{ return ARRAY_VIEW< typename T_ARRAY::ELEMENT const, typename T_ARRAY::INDEX >(a); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ARRAY_OPS_HPP
