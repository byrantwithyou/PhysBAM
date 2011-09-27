//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ARRAY_OPS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ARRAY_OPS_HPP

#include <cassert>

#include <algorithm>
#include <vector>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <Jeffrey_Utilities/PROPAGATE_CONST.h>

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
// struct ARRAY_SCALAR_VALUE< T_ARRAY >
//#####################################################################

template< class T_ARRAY >
struct ARRAY_SCALAR_VALUE
    : ARRAY_SCALAR_VALUE< typename ARRAY_VALUE< T_ARRAY >::type >
{ };

template< class T_ARRAY >
struct ARRAY_SCALAR_VALUE< const T_ARRAY >
    : ARRAY_SCALAR_VALUE< T_ARRAY >
{ };

template<>
struct ARRAY_SCALAR_VALUE< float >
{ typedef float type; };

template<>
struct ARRAY_SCALAR_VALUE< double >
{ typedef double type; };

//#####################################################################
// Size(const T_ARRAY& a) -> ...
//#####################################################################

template< class T, class T_DERIVED, class T_ID >
inline T_ID
Size(const ARRAY_BASE< T, T_DERIVED, T_ID >& a)
{ return static_cast< const T_DERIVED& >(a).Size(); }

template< class T, class T_DERIVED >
inline int
Size(const VECTOR_BASE< T, T_DERIVED >& v)
{ return static_cast< const T_DERIVED& >(v).Size(); }

template< class T, class T_ALLOCATOR >
inline typename std::vector< T, T_ALLOCATOR >::size_type
Size(const std::vector< T, T_ALLOCATOR >& v)
{ return v.size(); }

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

template< class T, class T_DERIVED >
inline T&
Front(VECTOR_BASE< T, T_DERIVED >& v)
{ return static_cast< T_DERIVED& >(v)(1); }

template< class T, class T_DERIVED >
inline const T&
Front(const VECTOR_BASE< T, T_DERIVED >& v)
{ return static_cast< const T_DERIVED& >(v)(1); }

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
{ return static_cast< T_DERIVED& >(a).Last(); }

template< class T, class T_DERIVED, class T_ID >
inline const T&
Back(const ARRAY_BASE< T, T_DERIVED, T_ID >& a)
{ return static_cast< const T_DERIVED& >(a).Last(); }

template< class T, class T_DERIVED >
inline T&
Back(VECTOR_BASE< T, T_DERIVED >& v)
{ return static_cast< T_DERIVED& >(v)(Size(v)); }

template< class T, class T_DERIVED >
inline const T&
Back(const VECTOR_BASE< T, T_DERIVED >& v)
{ return static_cast< const T_DERIVED& >(v)(Size(v)); }

template< class T, class T_ALLOCATOR >
inline typename std::vector< T, T_ALLOCATOR >::reference
Back(std::vector< T, T_ALLOCATOR >& v)
{ return v.back(); }

template< class T, class T_ALLOCATOR >
inline typename std::vector< T, T_ALLOCATOR >::const_reference
Back(const std::vector< T, T_ALLOCATOR >& v)
{ return v.back(); }

//#####################################################################
// At1(T_ARRAY& a, const ... index) -> ...
//#####################################################################

template< class T, class T_DERIVED, class T_ID >
inline T&
At1(ARRAY_BASE< T, T_DERIVED, T_ID >& a, const T_ID index)
{ return static_cast< T_DERIVED& >(a)(index); }

template< class T, class T_DERIVED, class T_ID >
inline const T&
At1(const ARRAY_BASE< T, T_DERIVED, T_ID >& a, const T_ID index)
{ return static_cast< const T_DERIVED& >(a)(index); }

template< class T, class T_DERIVED >
inline T&
At1(VECTOR_BASE< T, T_DERIVED >& v, const int index)
{ return static_cast< T_DERIVED& >(v)(index); }

template< class T, class T_DERIVED >
inline const T&
At1(const VECTOR_BASE< T, T_DERIVED >& v, const int index)
{ return static_cast< const T_DERIVED& >(v)(index); }

template< class T, class T_ALLOCATOR, class T_INDEX >
inline typename std::vector< T, T_ALLOCATOR >::reference
At1(std::vector< T, T_ALLOCATOR >& v, const T_INDEX index)
{
    typedef typename std::vector< T, T_ALLOCATOR >::size_type size_type;
    return v[static_cast< size_type >(index - 1)];
}

template< class T, class T_ALLOCATOR, class T_INDEX >
inline typename std::vector< T, T_ALLOCATOR >::const_reference
At1(const std::vector< T, T_ALLOCATOR >& v, const T_INDEX index)
{
    typedef typename std::vector< T, T_ALLOCATOR >::size_type size_type;
    return v[static_cast< size_type >(index - 1)];
}

//#####################################################################
// Fill(T_ARRAY& a, const U& x) -> void
//#####################################################################

template< class T, class T_DERIVED, class T_ID, class U >
inline void
Fill(ARRAY_BASE< T, T_DERIVED, T_ID >& a, const U& x)
{ ARRAYS_COMPUTATIONS::Fill(a, x); }

template< class T, class T_DERIVED, class U >
inline void
Fill(VECTOR_BASE< T, T_DERIVED >& v, const U& x)
{ static_cast< T_DERIVED& >(v).Fill(x); }

template< class T, class T_ALLOCATOR, class U >
inline void
Fill(std::vector< T, T_ALLOCATOR >& v, const U& x)
{ std::fill(v.begin(), v.end(), x); }

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
// Make_Array_View(int n, T* p) -> ARRAY_VIEW<T>
// Make_Const_Array_View(int n, const T* p) -> ARRAY_VIEW< const T >
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

template< class T >
inline ARRAY_VIEW<T>
Make_Array_View(const int n, T* const p)
{ return ARRAY_VIEW<T>(n, p); }

template< class T >
inline ARRAY_VIEW< const T >
Make_Const_Array_View(const int n, const T* const p)
{ return ARRAY_VIEW< const T >(n, p); }

//#####################################################################
// Scalar_Size(const T_ARRAY& a) -> ...
//#####################################################################

inline int
Scalar_Size(float)
{ return 1; }

inline int
Scalar_Size(double)
{ return 1; }

template< class T_ARRAY >
inline int
Scalar_Size(const T_ARRAY& a)
{ return Size(a) == 0 ? 0 : Size(a) * Scalar_Size(Front(a)); }

//#####################################################################
// Scalar_At1(T_ARRAY& a, const int index) -> ...
//#####################################################################

inline float&
Scalar_At1(float& x, const int index)
{
    static_cast<void>(index);
    assert(index == 1);
    return x;
}

inline const float&
Scalar_At1(const float& x, const int index)
{
    static_cast<void>(index);
    assert(index == 1);
    return x;
}

inline double&
Scalar_At1(double& x, const int index)
{
    static_cast<void>(index);
    assert(index == 1);
    return x;
}

inline const double&
Scalar_At1(const double& x, const int index)
{
    static_cast<void>(index);
    assert(index == 1);
    return x;
}

template< class T_ARRAY >
inline typename ARRAY_SCALAR_VALUE< T_ARRAY >::type &
Scalar_At1(T_ARRAY& a, const int index)
{
    const int size = Scalar_Size(Front(a));
    return Scalar_At1(At1(a, 1 + (index - 1) / size), 1 + (index - 1) % size);
}

template< class T_ARRAY >
inline typename ARRAY_SCALAR_VALUE< T_ARRAY >::type const &
Scalar_At1(const T_ARRAY& a, const int index)
{
    const int size = Scalar_Size(Front(a));
    return Scalar_At1(At1(a, 1 + (index - 1) / size), 1 + (index - 1) % size);
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ARRAY_OPS_HPP
