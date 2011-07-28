//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_STATIC_MULTI_INDEX_CUBE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_STATIC_MULTI_INDEX_CUBE_HPP

#include <cassert>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/bitxor.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/HAS_ISC_VALUE.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE_ITERATOR.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
class STATIC_MULTI_INDEX_CUBE_ITERATOR;

template< int D, int MIN_INDEX_, int MAX_INDEX_ /* = -MIN_INDEX_ */ >
struct STATIC_MULTI_INDEX_CUBE
{
    BOOST_MPL_ASSERT_RELATION( MIN_INDEX_, <=, MAX_INDEX_ );

    static const int DIMENSION = D;
    static const int MIN_INDEX = MIN_INDEX_;
    static const int MAX_INDEX = MAX_INDEX_;
    static const int WIDTH = 1 + MAX_INDEX - MIN_INDEX;
    static const int SIZE = STATIC_POW_C< WIDTH, DIMENSION >::value;

    template< class T_MULTI_INDEX > struct CONTAINS;
    template< class T_MULTI_INDEX > struct LINEAR_INDEX;
    template< int LINEAR_INDEX_ > struct MULTI_INDEX_C;
    template< class T_LINEAR_INDEX > struct MULTI_INDEX;

    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    static MULTI_INDEX_TYPE Min_Multi_Index();
    static MULTI_INDEX_TYPE Max_Multi_Index();
    static int Size();

    static bool Contains(const MULTI_INDEX_TYPE multi_index);
    template< class T_MULTI_INDEX >
    static typename boost::enable_if<
        boost::mpl::is_sequence< T_MULTI_INDEX >,
        bool
    >::type Contains(T_MULTI_INDEX);

    static MULTI_INDEX_TYPE Strides();

    static int Linear_Index(const MULTI_INDEX_TYPE& multi_index);
    template< class T_MULTI_INDEX >
    static typename boost::enable_if<
        boost::mpl::is_sequence< T_MULTI_INDEX >,
        int
    >::type Linear_Index(T_MULTI_INDEX);
    static MULTI_INDEX_TYPE Multi_Index(const int linear_index);
    template< class T_LINEAR_INDEX >
    static typename boost::enable_if<
        HAS_ISC_VALUE< T_LINEAR_INDEX >,
        MULTI_INDEX_TYPE
    >::type Multi_Index(T_LINEAR_INDEX);

    template<class> struct result;
    template< class T_THIS, class T >
    struct result< T_THIS ( T ) >
    {
        BOOST_MPL_ASSERT((boost::mpl::or_<
            boost::is_convertible< T, MULTI_INDEX_TYPE >,
            boost::mpl::is_sequence< typename REMOVE_QUALIFIERS<T>::type >,
            boost::is_convertible< T, int >,
            HAS_ISC_VALUE< typename REMOVE_QUALIFIERS<T>::type >
        >));
        typedef typename boost::mpl::if_<
            boost::mpl::or_<
                boost::is_convertible< T, MULTI_INDEX_TYPE >,
                boost::mpl::is_sequence< typename REMOVE_QUALIFIERS<T>::type >
            >,
            int,
            MULTI_INDEX_TYPE
        >::type type;
    };
    int operator()(const MULTI_INDEX_TYPE& multi_index) const;
    template< class T_MULTI_INDEX >
    typename boost::enable_if<
        boost::mpl::is_sequence< T_MULTI_INDEX >,
        int
    >::type operator()(T_MULTI_INDEX) const;
    MULTI_INDEX_TYPE operator()(const int linear_index) const;
    template< class T_LINEAR_INDEX >
    typename boost::enable_if<
        HAS_ISC_VALUE< T_LINEAR_INDEX >,
        MULTI_INDEX_TYPE
    >::type operator()(T_LINEAR_INDEX) const;

    typedef STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX, MAX_INDEX > iterator;
    typedef iterator const_iterator;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

namespace Detail_STATIC_MULTI_INDEX_CUBE
{

template< int D, int MIN_INDEX, int MAX_INDEX, class T_MULTI_INDEX >
struct CONTAINS_DISPATCH
{
    static const int LOCAL_INDEX = boost::mpl::at_c< T_MULTI_INDEX, D-1 >::type::value;
    static const bool value =
        MIN_INDEX <= LOCAL_INDEX && LOCAL_INDEX <= MAX_INDEX
     && CONTAINS_DISPATCH< D-1, MIN_INDEX, MAX_INDEX, T_MULTI_INDEX >::value;
};

template< int MIN_INDEX, int MAX_INDEX, class T_MULTI_INDEX >
struct CONTAINS_DISPATCH< 0, MIN_INDEX, MAX_INDEX, T_MULTI_INDEX >
{ static const bool value = true; };

} // namespace Detail_STATIC_MULTI_INDEX_CUBE

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< class T_MULTI_INDEX >
struct STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
CONTAINS
{
    BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_MULTI_INDEX >));
    typedef CONTAINS type;
    static const bool value = Detail_STATIC_MULTI_INDEX_CUBE::
        CONTAINS_DISPATCH< D, MIN_INDEX_, MAX_INDEX_, T_MULTI_INDEX >::value;
};

namespace Detail_STATIC_MULTI_INDEX_CUBE
{

template< int D, int MIN_INDEX, int MAX_INDEX, class T_MULTI_INDEX >
struct LINEAR_INDEX_DISPATCH
{
    static const int WIDTH = 1 + MAX_INDEX - MIN_INDEX;
    static const int LOCAL_INDEX = boost::mpl::at_c< T_MULTI_INDEX, D-1 >::type::value;
    static const int value =
        1 + (LOCAL_INDEX - MIN_INDEX)
      + WIDTH * (LINEAR_INDEX_DISPATCH< D-1, MIN_INDEX, MAX_INDEX, T_MULTI_INDEX >::value - 1);
};

template< int MIN_INDEX, int MAX_INDEX, class T_MULTI_INDEX >
struct LINEAR_INDEX_DISPATCH< 0, MIN_INDEX, MAX_INDEX, T_MULTI_INDEX >
{ static const int value = 1; };

} // namespace Detail_STATIC_MULTI_INDEX_CUBE

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< class T_MULTI_INDEX >
struct STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
LINEAR_INDEX
{
    BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_MULTI_INDEX >));
    BOOST_MPL_ASSERT((CONTAINS< T_MULTI_INDEX >));
    typedef LINEAR_INDEX type;
    static const int value = Detail_STATIC_MULTI_INDEX_CUBE::
        LINEAR_INDEX_DISPATCH< D, MIN_INDEX_, MAX_INDEX_, T_MULTI_INDEX >::value;
};

namespace Detail_STATIC_MULTI_INDEX_CUBE
{

template< int D, int MIN_INDEX, int MAX_INDEX, int LINEAR_INDEX >
struct MULTI_INDEX_DISPATCH
{
    static const int WIDTH = 1 + (MAX_INDEX - MIN_INDEX);
    typedef typename boost::mpl::push_back<
        typename MULTI_INDEX_DISPATCH< D-1, MIN_INDEX, MAX_INDEX, 1 + (LINEAR_INDEX - 1) / WIDTH >::type,
        boost::mpl::int_< MIN_INDEX + (LINEAR_INDEX - 1) % WIDTH >
    >::type type;
};

template< int MIN_INDEX, int MAX_INDEX, int LINEAR_INDEX >
struct MULTI_INDEX_DISPATCH< 0, MIN_INDEX, MAX_INDEX, LINEAR_INDEX >
{ typedef boost::mpl::vector0<> type; };

} // namespace Detail_STATIC_MULTI_INDEX_CUBE

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< int LINEAR_INDEX_ >
struct STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
MULTI_INDEX_C
{
    typedef typename Detail_STATIC_MULTI_INDEX_CUBE::
        MULTI_INDEX_DISPATCH< D, MIN_INDEX, MAX_INDEX, LINEAR_INDEX_ >::type type;
};

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< class T_LINEAR_INDEX >
struct STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
MULTI_INDEX
    : MULTI_INDEX_C< T_LINEAR_INDEX::value >
{ };

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::MULTI_INDEX_TYPE
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Min_Multi_Index()
{ return MIN_INDEX_ * MULTI_INDEX_TYPE::All_Ones_Vector(); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::MULTI_INDEX_TYPE
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Max_Multi_Index()
{ return MAX_INDEX_ * MULTI_INDEX_TYPE::All_Ones_Vector(); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline int
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Size()
{ return SIZE; }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline bool
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Contains(const MULTI_INDEX_TYPE multi_index)
{
    for(int d = 1; d <= D; ++d)
        if(multi_index[d] < MIN_INDEX || MAX_INDEX < multi_index[d])
            return false;
    return true;
}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< class T_MULTI_INDEX >
inline typename boost::enable_if<
    boost::mpl::is_sequence< T_MULTI_INDEX >,
    bool
>::type
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Contains(T_MULTI_INDEX)
{ return CONTAINS< T_MULTI_INDEX >::value; }

namespace Detail_STATIC_MULTI_INDEX_CUBE
{

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
struct DISPATCH;

template< int MIN_INDEX_, int MAX_INDEX_ >
struct DISPATCH< 1, MIN_INDEX_, MAX_INDEX_ >
{
    static const int WIDTH = 1 + MAX_INDEX_ - MIN_INDEX_;
    static VECTOR<int,1> Strides()
    { return VECTOR<int,1>(1); }
    static int Linear_Index(const VECTOR<int,1>& multi_index)
    { return multi_index[1] + (1 - MIN_INDEX_); }
    static VECTOR<int,1> Multi_Index(const int linear_index)
    { return VECTOR<int,1>(linear_index + (MIN_INDEX_ - 1)); }
};

template< int MIN_INDEX_, int MAX_INDEX_ >
struct DISPATCH< 2, MIN_INDEX_, MAX_INDEX_ >
{
    static const int WIDTH = 1 + MAX_INDEX_ - MIN_INDEX_;
    static VECTOR<int,2> Strides()
    { return VECTOR<int,2>(WIDTH, 1); }
    static int Linear_Index(const VECTOR<int,2>& multi_index)
    {
        return WIDTH * multi_index[1] + multi_index[2]
             + (1 - (1 + WIDTH) * MIN_INDEX_);
    }
    static VECTOR<int,2> Multi_Index(const int linear_index)
    {
        return VECTOR<int,2>(
            MIN_INDEX_ + (linear_index - 1) / WIDTH,
            MIN_INDEX_ + (linear_index - 1) % WIDTH
        );
    }
};

template< int MIN_INDEX_, int MAX_INDEX_ >
struct DISPATCH< 3, MIN_INDEX_, MAX_INDEX_ >
{
    static const int WIDTH = 1 + MAX_INDEX_ - MIN_INDEX_;
    static VECTOR<int,3> Strides()
    { return VECTOR<int,3>(WIDTH * WIDTH, WIDTH, 1); }
    static int Linear_Index(const VECTOR<int,3>& multi_index)
    {
        return (WIDTH * WIDTH) * multi_index[1] + WIDTH * multi_index[2] + multi_index[3]
             + (1 - (1 + WIDTH + WIDTH * WIDTH) * MIN_INDEX_);
    }
    static VECTOR<int,3> Multi_Index(const int linear_index)
    {
        return VECTOR<int,3>(
            MIN_INDEX_ + (linear_index - 1) / (WIDTH * WIDTH),
            MIN_INDEX_ + (linear_index - 1) / WIDTH % WIDTH,
            MIN_INDEX_ + (linear_index - 1) % WIDTH
        );
    }
};

} // namespace Detail_STATIC_MULTI_INDEX_CUBE

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::MULTI_INDEX_TYPE
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Strides()
{ return Detail_STATIC_MULTI_INDEX_CUBE::DISPATCH< D, MIN_INDEX, MAX_INDEX >::Strides(); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline int
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Linear_Index(const MULTI_INDEX_TYPE& multi_index)
{
#ifndef NDEBUG
    for(int d = 1; d <= D; ++d)
        assert(MIN_INDEX <= multi_index[d] && multi_index[d] <= MAX_INDEX);
#endif // #ifndef NDEBUG
    return Detail_STATIC_MULTI_INDEX_CUBE::DISPATCH< D, MIN_INDEX, MAX_INDEX >::Linear_Index(multi_index);
}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< class T_MULTI_INDEX >
inline typename boost::enable_if<
    boost::mpl::is_sequence< T_MULTI_INDEX >,
    int
>::type
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Linear_Index(T_MULTI_INDEX)
{ return LINEAR_INDEX< T_MULTI_INDEX >::value; }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::MULTI_INDEX_TYPE
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Multi_Index(const int linear_index)
{
    assert(1 <= linear_index && linear_index <= (STATIC_POW_C< WIDTH, DIMENSION >::value));
    return Detail_STATIC_MULTI_INDEX_CUBE::DISPATCH< D, MIN_INDEX, MAX_INDEX >::Multi_Index(linear_index);
}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< class T_LINEAR_INDEX >
inline typename boost::enable_if<
    HAS_ISC_VALUE< T_LINEAR_INDEX >,
    typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::MULTI_INDEX_TYPE
>::type
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
Multi_Index(T_LINEAR_INDEX)
{ return typename MULTI_INDEX< T_LINEAR_INDEX >::type(); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline int
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
operator()(const MULTI_INDEX_TYPE& multi_index) const
{ return Linear_Index(multi_index); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< class T_MULTI_INDEX >
inline typename boost::enable_if<
    boost::mpl::is_sequence< T_MULTI_INDEX >,
    int
>::type
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
operator()(T_MULTI_INDEX) const
{ return LINEAR_INDEX< T_MULTI_INDEX >::value; }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::MULTI_INDEX_TYPE
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
operator()(const int linear_index) const
{ return Multi_Index(linear_index); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
template< class T_LINEAR_INDEX >
inline typename boost::enable_if<
    HAS_ISC_VALUE< T_LINEAR_INDEX >,
    typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::MULTI_INDEX_TYPE
>::type
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
operator()(T_LINEAR_INDEX) const
{ return typename MULTI_INDEX< T_LINEAR_INDEX >::type(); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::iterator
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
begin() const
{ return iterator(BEGIN_TAG()); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::iterator
STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ >::
end() const
{ return iterator(END_TAG()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_STATIC_MULTI_INDEX_CUBE_HPP
