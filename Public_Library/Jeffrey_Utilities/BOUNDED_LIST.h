//#####################################################################
// Copyright 2009, Jeffrey Hellrung, Jan Hegemann.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// struct BOUNDED_LIST
//#####################################################################
// A BOUNDED_LIST is like a ARRAY with a compile-time bounded
// size.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_BOUNDED_LIST_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_BOUNDED_LIST_HPP

#include <boost/mpl/assert.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int N >
struct BOUNDED_LIST
{
    typedef T ELEMENT;
    typedef int INDEX;

    BOUNDED_LIST();

    template<class T_ARRAY>
    explicit BOUNDED_LIST(const T_ARRAY& array);

    int Size() const;
    T* Data();
    const T* Data() const;

    /***/ T& operator()(int i);
    const T& operator()(int i) const;
    /***/ T& Last();
    const T& Last() const;

    int Append(const T& x);
    int Append_Unique(const T& x);
    void Insert(const T& x, const int i);

    int Remove_Index(const int index);
    void Remove_All();

    template< class T_INDICES >
    INDIRECT_ARRAY< BOUNDED_LIST, const T_INDICES& > Subset(const T_INDICES& indices);
    template< class T_INDICES >
    INDIRECT_ARRAY< const BOUNDED_LIST, const T_INDICES& > Subset(const T_INDICES& indices) const;

    typedef /***/ T* /***/ iterator;
    typedef const T* const_iterator;
    /***/ iterator begin();
    const_iterator begin() const;
    /***/ iterator end();
    const_iterator end() const;

private:
    VECTOR<T,N> m_values;
    int m_size;
};

//#####################################################################
//#####################################################################

template< class T, int N >
inline
BOUNDED_LIST<T,N>::
BOUNDED_LIST()
    : m_size(0)
{ }

template< class T, int N >
template<class T_ARRAY>
inline
BOUNDED_LIST<T,N>::
BOUNDED_LIST(const T_ARRAY& array)
    : m_size(0)
{
    for(int i = 1; i <= array.Size(); ++i)
        Append(array(i));
}

template< class T, int N >
inline int
BOUNDED_LIST<T,N>::
Size() const
{ return m_size; }

template< class T, int N >
inline T*
BOUNDED_LIST<T,N>::
Data()
{ return m_values.array; }

template< class T, int N >
inline const T*
BOUNDED_LIST<T,N>::
Data() const
{ return m_values.array; }

template< class T, int N >
inline T&
BOUNDED_LIST<T,N>::
operator()(int i)
{
    assert(1 <= i);
    assert(i <= m_size);
    return m_values[i];
}

template< class T, int N >
inline const T&
BOUNDED_LIST<T,N>::
operator()(int i) const
{
    assert(1 <= i);
    assert(i <= m_size);
    return m_values[i];
}

template< class T, int N >
inline T&
BOUNDED_LIST<T,N>::
Last()
{
    BOOST_MPL_ASSERT_RELATION( N, >, 0 );
    return m_values[m_size];
}

template< class T, int N >
inline const T&
BOUNDED_LIST<T,N>::
Last() const
{
    BOOST_MPL_ASSERT_RELATION( N, >, 0 );
    return m_values[m_size];
}

template< class T, int N >
inline int
BOUNDED_LIST<T,N>::
Append(const T& x)
{
    m_values[++m_size] = x;
    return m_size;
}

template< class T, int N >
inline int
BOUNDED_LIST<T,N>::
Append_Unique(const T& x)
{
    for(int i = 1; i <= m_size; ++i)
        if(m_values[i] == x)
            return m_size;
    return Append(x);
}

template< class T, int N >
inline void
BOUNDED_LIST<T,N>::
Insert(const T& x, const int i)
{
    for(int j = ++m_size; j > i; --j)
        m_values[j] = m_values[j-1];
    m_values[i] = x;
}

template< class T, int N >
inline int
BOUNDED_LIST<T,N>::
Remove_Index(const int index)
{
    for(int i = index; i < m_size; ++i)
        m_values[i] = m_values[i + 1];
    m_values[m_size] = T();
    return --m_size;
}

template< class T, int N >
inline void
BOUNDED_LIST<T,N>::
Remove_All()
{
    m_values.Fill(T());
    m_size = 0;
}

template< class T, int N >
template< class T_INDICES >
inline INDIRECT_ARRAY< BOUNDED_LIST<T,N>, const T_INDICES& >
BOUNDED_LIST<T,N>::
Subset(const T_INDICES& indices)
{ return INDIRECT_ARRAY< BOUNDED_LIST, const T_INDICES& >(*this, indices); }

template< class T, int N >
template< class T_INDICES >
inline INDIRECT_ARRAY< const BOUNDED_LIST<T,N>, const T_INDICES& >
BOUNDED_LIST<T,N>::
Subset(const T_INDICES& indices) const
{ return INDIRECT_ARRAY< const BOUNDED_LIST, const T_INDICES& >(*this, indices); }

template< class T, int N >
inline typename BOUNDED_LIST<T,N>::iterator
BOUNDED_LIST<T,N>::
begin()
{ return &m_values[1]; }

template< class T, int N >
inline typename BOUNDED_LIST<T,N>::const_iterator
BOUNDED_LIST<T,N>::
begin() const
{ return &m_values[1]; }

template< class T, int N >
inline typename BOUNDED_LIST<T,N>::iterator
BOUNDED_LIST<T,N>::
end()
{ return &m_values[1] + m_size; }

template< class T, int N >
inline typename BOUNDED_LIST<T,N>::const_iterator
BOUNDED_LIST<T,N>::
end() const
{ return &m_values[1] + m_size; }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_BOUNDED_LIST_HPP
