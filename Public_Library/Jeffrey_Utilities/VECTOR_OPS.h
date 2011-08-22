//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Seriously, why is this needed???
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VECTOR_OPS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VECTOR_OPS_HPP

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_EXPRESSION.h>

namespace PhysBAM
{

//#####################################################################
// As_Vector(const T (&array)[N]) -> VECTOR<T,N>
// As_Vector<T>(const U (&array)[N]) -> VECTOR<T,N>
//#####################################################################

template< class T >
inline VECTOR<T,1>
As_Vector(const T (&array)[1])
{ return VECTOR<T,1>(array[0]); }

template< class T >
inline VECTOR<T,2>
As_Vector(const T (&array)[2])
{ return VECTOR<T,2>(array[0], array[1]); }

template< class T >
inline VECTOR<T,3>
As_Vector(const T (&array)[3])
{ return VECTOR<T,3>(array[0], array[1], array[2]); }

template< class T, int N >
inline VECTOR<T,N>
As_Vector(const T (&array)[N])
{
    VECTOR<T,N> result;
    for(int i = 1; i <= N; ++i)
        result[i] = array[i-1];
    return result;
}

template< class T, class U >
inline VECTOR<T,1>
As_Vector(const U (&array)[1])
{ return VECTOR<T,1>(static_cast<T>(array[0])); }

template< class T, class U >
inline VECTOR<T,2>
As_Vector(const U (&array)[2])
{ return VECTOR<T,2>(static_cast<T>(array[0]), static_cast<T>(array[1])); }

template< class T, class U >
inline VECTOR<T,3>
As_Vector(const U (&array)[3])
{ return VECTOR<T,3>(static_cast<T>(array[0]), static_cast<T>(array[1]), static_cast<T>(array[2])); }

template< class T, class U, int N >
inline VECTOR<T,N>
As_Vector(const U (&array)[N])
{
    VECTOR<T,N> result;
    for(int i = 1; i <= N; ++i)
        result[i] = static_cast<T>(array[i-1]);
    return result;
}

//#####################################################################
// Count(const VECTOR<T,N>& v, const T& x) -> int
//#####################################################################

template< class T, int N >
inline int
Count(const VECTOR<T,N>& v, const T& x)
{
    int count = 0;
    for(int i = 1; i <= N; ++i)
        if(v[i] == x)
            ++count;
    return count;
}

//#####################################################################
// Replace(VECTOR<T,N>& v, const T& x, const T& y) -> void
//#####################################################################

template< class T, int N >
inline void
Replace(VECTOR<T,N>& v, const T& x, const T& y)
{
    for(int i = 1; i <= N; ++i)
        if(v[i] == x)
            v[i] = y;
}

//#####################################################################
// operator*(const int c, VECTOR<T,N> v) -> VECTOR<T,N>
// operator*(VECTOR<T,N> v, const int c) -> VECTOR<T,N>
// operator*(const VECTOR<int,N>& c, VECTOR<T,N> v) -> VECTOR<T,N>
// operator*(VECTOR<T,N> v, const VECTOR<int,N>& c) -> VECTOR<T,N>
// operator/(VECTOR<T,N> v, const VECTOR<int,N>& c) -> VECTOR<T,N>
// operator/(const VECTOR<int,N>& c, VECTOR<T,N> v) -> VECTOR<T,N>
//#####################################################################

template< class T, int N >
inline typename boost::disable_if_c<
    (N > 3) || boost::is_same< T, int >::value,
    VECTOR<T,N>
>::type
operator*(const int c, VECTOR<T,N> v)
{
    for(int i = 1; i <= N; ++i)
        v[i] *= c;
    return v;
}

template< class T_C, class T, int N >
inline typename boost::disable_if_c<
    (N <= 3)
 || !boost::is_same< T_C, int >::value
 || boost::is_same< T, int >::value,
    VECTOR_SCALE< int, VECTOR<T,N> >
>::type
operator*(const T_C c, const VECTOR<T,N>& v)
{ return VECTOR_SCALE< int, VECTOR<T,N> >(c, v); }

template< class T, int N >
inline typename boost::disable_if_c<
    (N > 3) || boost::is_same< T, int >::value,
    VECTOR<T,N>
>::type
operator*(VECTOR<T,N> v, const int c)
{
    for(int i = 1; i <= N; ++i)
        v[i] *= c;
    return v;
}

template< class T, int N, class T_C >
inline typename boost::disable_if_c<
    (N <= 3)
 || boost::is_same< T, int >::value
 || !boost::is_same< T_C, int >::value,
    VECTOR_SCALE< int, VECTOR<T,N> >
>::type
operator*(const VECTOR<T,N>& v, const T_C c)
{ return VECTOR_SCALE< int, VECTOR<T,N> >(c, v); }

template< class T, int N >
inline typename boost::disable_if<
    boost::is_same< T, int >,
    VECTOR<T,N>
>::type
operator*(const VECTOR<int,N> c, VECTOR<T,N> v)
{
    for(int i = 1; i <= N; ++i)
        v[i] *= c[i];
    return v;
}

template< class T, int N >
inline typename boost::disable_if<
    boost::is_same< T, int >,
    VECTOR<T,N>
>::type
operator*(VECTOR<T,N> v, const VECTOR<int,N> c)
{
    for(int i = 1; i <= N; ++i)
        v[i] *= c[i];
    return v;
}

template< class T, int N >
inline typename boost::disable_if<
    boost::is_same< T, int >,
    VECTOR<T,N>
>::type
operator/(VECTOR<T,N> v, const VECTOR<int,N> c)
{
    for(int i = 1; i <= N; ++i)
        v[i] /= c[i];
    return v;
}

template< class T, int N >
inline typename boost::disable_if<
    boost::is_same< T, int >,
    VECTOR<T,N>
>::type
operator/(const VECTOR<int,N> c, VECTOR<T,N> v)
{
    for(int i = 1; i <= N; ++i)
        v[i] = c[i] / v[i];
    return v;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VECTOR_OPS_HPP
