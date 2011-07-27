//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CUBE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CUBE_HPP

#include <boost/mpl/assert.hpp>
#include <boost/mpl/bitxor.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE_ITERATOR.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
class MULTI_INDEX_CUBE_ITERATOR;

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ /* = -MIN_OFFSET_ */ >
struct MULTI_INDEX_CUBE
{
    BOOST_MPL_ASSERT_RELATION( MIN_OFFSET_, <=, MAX_OFFSET_ );

    static const int DIMENSION = D;
    static const int MIN_OFFSET = MIN_OFFSET_;
    static const int MAX_OFFSET = MAX_OFFSET_;
    static const int WIDTH = 1 + MAX_OFFSET - MIN_OFFSET;
    static const int SIZE = STATIC_POW_C< WIDTH, DIMENSION >::value;
    typedef STATIC_MULTI_INDEX_CUBE< D, MIN_OFFSET, MAX_OFFSET > STATIC_MULTI_OFFSET_CUBE_TYPE;

    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    MULTI_INDEX_TYPE base_multi_index;

    MULTI_INDEX_CUBE();
    explicit MULTI_INDEX_CUBE(const MULTI_INDEX_TYPE& base_multi_index_);

    MULTI_INDEX_TYPE Min_Multi_Index() const;
    MULTI_INDEX_TYPE Max_Multi_Index() const;
    static int Size();

    bool Contains(const MULTI_INDEX_TYPE& multi_index) const;

    static MULTI_INDEX_TYPE Strides();

    MULTI_INDEX_TYPE Multi_Offset(const MULTI_INDEX_TYPE& multi_index) const;

    int Linear_Index(const MULTI_INDEX_TYPE& multi_index) const;
    MULTI_INDEX_TYPE Multi_Index(const int linear_index) const;

    template<class> struct result;
    template< class T_THIS, class T >
    struct result< T_THIS ( T ) >
    {
        BOOST_MPL_ASSERT((boost::mpl::bitxor_<
            boost::is_convertible< T, MULTI_INDEX_TYPE >,
            boost::is_convertible< T, int >
        >));
        typedef typename boost::mpl::if_<
            boost::is_convertible< T, MULTI_INDEX_TYPE >,
            int,
            MULTI_INDEX_TYPE
        >::type type;
    };
    int operator()(const MULTI_INDEX_TYPE& multi_index) const;
    MULTI_INDEX_TYPE operator()(const int linear_index) const;

    typedef MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET, MAX_OFFSET > iterator;
    typedef iterator const_iterator;
    iterator begin() const;
    iterator end() const;
};

template<
    int D,
    int MIN_OFFSET1, int MAX_OFFSET1,
    int MIN_OFFSET2, int MAX_OFFSET2
>
inline typename boost::enable_if_c<
    MAX_OFFSET1 - MIN_OFFSET1 == MAX_OFFSET2 - MIN_OFFSET2,
    bool
>::type
operator==(
    const MULTI_INDEX_CUBE< D, MIN_OFFSET1, MAX_OFFSET1 >& cube1,
    const MULTI_INDEX_CUBE< D, MIN_OFFSET2, MAX_OFFSET2 >& cube2)
{
    for(int d = 1; d <= D; ++d)
        if(cube1.base_multi_index[d] - cube2.base_multi_index[d] != MIN_OFFSET2 - MIN_OFFSET1)
            return false;
    return true;
}

template<
    int D,
    int MIN_OFFSET1, int MAX_OFFSET1,
    int MIN_OFFSET2, int MAX_OFFSET2
>
inline typename boost::enable_if_c<
    MAX_OFFSET1 - MIN_OFFSET1 != MAX_OFFSET2 - MIN_OFFSET2,
    bool
>::type
operator==(
    const MULTI_INDEX_CUBE< D, MIN_OFFSET1, MAX_OFFSET1 >& /*cube1*/,
    const MULTI_INDEX_CUBE< D, MIN_OFFSET2, MAX_OFFSET2 >& /*cube2*/)
{ return false; }

template<
    int D,
    int MIN_OFFSET1, int MAX_OFFSET1,
    int MIN_OFFSET2, int MAX_OFFSET2
>
inline bool
operator!=(
    const MULTI_INDEX_CUBE< D, MIN_OFFSET1, MAX_OFFSET1 >& cube1,
    const MULTI_INDEX_CUBE< D, MIN_OFFSET2, MAX_OFFSET2 >& cube2)
{ return !operator==(cube1, cube2); }

//#####################################################################
//#####################################################################

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE()
    : base_multi_index() // init'ed to 0
{ }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE(const MULTI_INDEX_TYPE& base_multi_index_)
    : base_multi_index(base_multi_index_)
{ }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_TYPE
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
Min_Multi_Index() const
{ return base_multi_index + MIN_OFFSET_; }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_TYPE
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
Max_Multi_Index() const
{ return base_multi_index + MAX_OFFSET_; }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline int
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
Size()
{ return SIZE; }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline bool
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
Contains(const MULTI_INDEX_TYPE& multi_index) const
{ return STATIC_MULTI_OFFSET_CUBE_TYPE::Contains(multi_index - base_multi_index); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_TYPE
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
Strides()
{ return STATIC_MULTI_OFFSET_CUBE_TYPE::Strides(); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_TYPE
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
Multi_Offset(const MULTI_INDEX_TYPE& multi_index) const
{ return multi_index - base_multi_index; }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline int
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
Linear_Index(const MULTI_INDEX_TYPE& multi_index) const
{ return STATIC_MULTI_OFFSET_CUBE_TYPE::Linear_Index(multi_index - base_multi_index); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_TYPE
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
Multi_Index(const int linear_index) const
{ return base_multi_index + STATIC_MULTI_OFFSET_CUBE_TYPE::Multi_Index(linear_index); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline int
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
operator()(const MULTI_INDEX_TYPE& multi_index) const
{ return Linear_Index(multi_index); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_TYPE
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
operator()(const int linear_index) const
{ return Multi_Index(linear_index); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::iterator
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
begin() const
{ return iterator(*this, BEGIN_TAG()); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::iterator
MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ >::
end() const
{ return iterator(*this, END_TAG()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CUBE_HPP
