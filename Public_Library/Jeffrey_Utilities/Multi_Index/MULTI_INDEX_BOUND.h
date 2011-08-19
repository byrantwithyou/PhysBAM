//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOUND_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOUND_HPP

#include <cassert>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/bitxor.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int D >
class MULTI_INDEX_BOUND_ITERATOR;

template< int D >
struct MULTI_INDEX_BOUND
{
    static const int DIMENSION = D;
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    MULTI_INDEX_TYPE max_multi_index;

    MULTI_INDEX_BOUND();
    MULTI_INDEX_BOUND(const MULTI_INDEX_TYPE& max_multi_index_);

    static MULTI_INDEX_TYPE Min_Multi_Index();
    const MULTI_INDEX_TYPE& Max_Multi_Index() const;
    int Size() const;

    bool Contains(const MULTI_INDEX_TYPE multi_index) const;
    MULTI_INDEX_TYPE Clamp(MULTI_INDEX_TYPE multi_index) const;

    MULTI_INDEX_TYPE Strides() const;

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

    typedef MULTI_INDEX_BOUND_ITERATOR<D> iterator;
    typedef iterator const_iterator;
    iterator begin() const;
    iterator end() const;
};

template< int D >
inline bool
operator==(const MULTI_INDEX_BOUND<D>& bound1, const MULTI_INDEX_BOUND<D>& bound2)
{ return bound1.max_multi_index == bound2.max_multi_index; }

template< int D >
inline bool
operator!=(const MULTI_INDEX_BOUND<D>& bound1, const MULTI_INDEX_BOUND<D>& bound2)
{ return !operator==(bound1, bound2); }

template< int D >
inline MULTI_INDEX_BOUND<D>
operator+(const MULTI_INDEX_BOUND<D>& bound, const int shift)
{ return MULTI_INDEX_BOUND<D>(bound.max_multi_index + shift); }

template< int D >
inline MULTI_INDEX_BOUND<D>
operator+(const int shift, const MULTI_INDEX_BOUND<D>& bound)
{ return MULTI_INDEX_BOUND<D>(shift + bound.max_multi_index); }

template< int D >
inline MULTI_INDEX_BOUND<D>
operator-(const MULTI_INDEX_BOUND<D>& bound, const int inv_shift)
{ return MULTI_INDEX_BOUND<D>(bound.max_multi_index - inv_shift); }

template< int D >
inline MULTI_INDEX_BOUND<D>
operator*(const MULTI_INDEX_BOUND<D>& bound, const int scaling)
{ return MULTI_INDEX_BOUND<D>(bound.max_multi_index * scaling); }

template< int D >
inline MULTI_INDEX_BOUND<D>
operator*(const int scaling, const MULTI_INDEX_BOUND<D>& bound)
{ return MULTI_INDEX_BOUND<D>(scaling * bound.max_multi_index); }

template< int D >
inline MULTI_INDEX_BOUND<D>
operator/(const MULTI_INDEX_BOUND<D>& bound, const int inv_scaling)
{ return MULTI_INDEX_BOUND<D>(bound.max_multi_index / inv_scaling); }

//#####################################################################
//#####################################################################

template< int D >
inline
MULTI_INDEX_BOUND<D>::
MULTI_INDEX_BOUND()
    : max_multi_index(MULTI_INDEX_TYPE::All_Ones_Vector())
{ }

template< int D >
inline
MULTI_INDEX_BOUND<D>::
MULTI_INDEX_BOUND(const MULTI_INDEX_TYPE& max_multi_index_)
    : max_multi_index(max_multi_index_)
{ }

template< int D >
inline typename MULTI_INDEX_BOUND<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOUND<D>::
Min_Multi_Index()
{ return MULTI_INDEX_TYPE::All_Ones_Vector(); }

template< int D >
inline typename MULTI_INDEX_BOUND<D>::MULTI_INDEX_TYPE const &
MULTI_INDEX_BOUND<D>::
Max_Multi_Index() const
{ return max_multi_index; }

template< int D >
inline int
MULTI_INDEX_BOUND<D>::
Size() const
{ return max_multi_index.Product(); }

template< int D >
inline bool
MULTI_INDEX_BOUND<D>::
Contains(const MULTI_INDEX_TYPE multi_index) const
{
    for(int d = 1; d <= D; ++d)
        if(multi_index[d] < 1 || max_multi_index[d] < multi_index[d])
            return false;
    return true;
}

template< int D >
inline typename MULTI_INDEX_BOUND<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOUND<D>::
Clamp(MULTI_INDEX_TYPE multi_index) const
{
    for(int d = 1; d <= D; ++d) {
        if(multi_index[d] < 1)
            multi_index[d] = 1;
        else if(max_multi_index[d] < multi_index[d])
            multi_index[d] = max_multi_index[d];
    }
    return multi_index;
}

namespace Detail_MULTI_INDEX_BOUND
{

template< int D >
struct DISPATCH;

template<>
struct DISPATCH<1>
{
    static VECTOR<int,1> Strides(const MULTI_INDEX_BOUND<1>& /*this_*/)
    { return VECTOR<int,1>(1); }
    static int Linear_Index(const MULTI_INDEX_BOUND<1>& this_, const VECTOR<int,1>& multi_index)
    { return multi_index[1]; }
    static VECTOR<int,1> Multi_Index(const MULTI_INDEX_BOUND<1>& this_, const int linear_index)
    { return VECTOR<int,1>(linear_index); }
};

template<>
struct DISPATCH<2>
{
    static VECTOR<int,2> Strides(const MULTI_INDEX_BOUND<2>& this_)
    { return VECTOR<int,2>(this_.max_multi_index[2], 1); }
    static int Linear_Index(const MULTI_INDEX_BOUND<2>& this_, const VECTOR<int,2>& multi_index)
    {
        return 1 +                        (multi_index[2] - 1) +
               this_.max_multi_index[2] * (multi_index[1] - 1);
    }
    static VECTOR<int,2> Multi_Index(const MULTI_INDEX_BOUND<2>& this_, const int linear_index)
    {
        return VECTOR<int,2>(
            1 + (linear_index - 1) / this_.max_multi_index[2],
            1 + (linear_index - 1) % this_.max_multi_index[2]
        );
    }
};

template<>
struct DISPATCH<3>
{
    static VECTOR<int,3> Strides(const MULTI_INDEX_BOUND<3>& this_)
    { return VECTOR<int,3>(this_.max_multi_index[2] * this_.max_multi_index[3], this_.max_multi_index[3], 1); }
    static int Linear_Index(const MULTI_INDEX_BOUND<3>& this_, const VECTOR<int,3>& multi_index)
    {
        return 1 +                         (multi_index[3] - 1) +
               this_.max_multi_index[3] * ((multi_index[2] - 1) +
               this_.max_multi_index[2] *  (multi_index[1] - 1));
    }
    static VECTOR<int,3> Multi_Index(const MULTI_INDEX_BOUND<3>& this_, const int linear_index)
    {
        return VECTOR<int,3>(
            1 + (linear_index - 1) / (this_.max_multi_index[3] * this_.max_multi_index[2]),
            1 + (linear_index - 1) / this_.max_multi_index[3] % this_.max_multi_index[2],
            1 + (linear_index - 1) % this_.max_multi_index[3]
        );
    }
};

} // namespace Detail_MULTI_INDEX_BOUND

template< int D >
inline typename MULTI_INDEX_BOUND<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOUND<D>::
Strides() const
{ return Detail_MULTI_INDEX_BOUND::DISPATCH<D>::Strides(*this); }

template< int D >
inline int
MULTI_INDEX_BOUND<D>::
Linear_Index(const MULTI_INDEX_TYPE& multi_index) const
{
#ifndef NDEBUG
    for(int d = 1; d <= D; ++d)
        assert(1 <= multi_index[d] && multi_index[d] <= max_multi_index[d]);
#endif // #ifndef NDEBUG
    return Detail_MULTI_INDEX_BOUND::DISPATCH<D>::Linear_Index(*this, multi_index);
}

template< int D >
inline typename MULTI_INDEX_BOUND<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOUND<D>::
Multi_Index(const int linear_index) const
{
    assert(1 <= linear_index && linear_index <= Size());
    return Detail_MULTI_INDEX_BOUND::DISPATCH<D>::Multi_Index(*this, linear_index);
}

template< int D >
inline int
MULTI_INDEX_BOUND<D>::
operator()(const MULTI_INDEX_TYPE& multi_index) const
{ return Linear_Index(multi_index); }

template< int D >
inline typename MULTI_INDEX_BOUND<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOUND<D>::
operator()(const int linear_index) const
{ return Multi_Index(linear_index); }

template< int D >
inline typename MULTI_INDEX_BOUND<D>::iterator
MULTI_INDEX_BOUND<D>::
begin() const
{ return iterator(*this, BEGIN_TAG()); }

template< int D >
inline typename MULTI_INDEX_BOUND<D>::iterator
MULTI_INDEX_BOUND<D>::
end() const
{ return iterator(*this, END_TAG()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOUND_HPP
