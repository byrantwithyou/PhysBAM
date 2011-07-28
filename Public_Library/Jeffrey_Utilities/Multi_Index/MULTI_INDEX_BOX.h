//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_HPP

#include <cassert>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/bitxor.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOX_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int D >
class MULTI_INDEX_BOX_ITERATOR;

template< int D >
struct MULTI_INDEX_BOX
{
    static const int DIMENSION = D;
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    MULTI_INDEX_TYPE min_multi_index;
    MULTI_INDEX_TYPE max_multi_index;
    MULTI_INDEX_TYPE widths;

    MULTI_INDEX_BOX();
    MULTI_INDEX_BOX(const MULTI_INDEX_TYPE& multi_index);
    MULTI_INDEX_BOX(
        const MULTI_INDEX_TYPE& min_multi_index_,
        const MULTI_INDEX_TYPE& max_multi_index_);

    const MULTI_INDEX_TYPE& Min_Multi_Index() const;
    const MULTI_INDEX_TYPE& Max_Multi_Index() const;
    int Size() const;

    bool Contains(const MULTI_INDEX_TYPE multi_index) const;

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

    typedef MULTI_INDEX_BOX_ITERATOR<D> iterator;
    typedef iterator const_iterator;
    iterator begin() const;
    iterator end() const;
};

template< int D >
inline bool
operator==(const MULTI_INDEX_BOX<D>& box1, const MULTI_INDEX_BOX<D>& box2)
{ return box1.min_multi_index == box2.min_multi_index && box1.max_multi_index == box2.max_multi_index; }

template< int D >
inline bool
operator!=(const MULTI_INDEX_BOX<D>& box1, const MULTI_INDEX_BOX<D>& box2)
{ return !operator==(box1, box2); }

//#####################################################################
//#####################################################################

template< int D >
inline
MULTI_INDEX_BOX<D>::
MULTI_INDEX_BOX()
    : min_multi_index(), // init'ed to 0
      max_multi_index(), // init'ed to 0
      widths(MULTI_INDEX_TYPE::All_Ones_Vector())
{ }

template< int D >
inline
MULTI_INDEX_BOX<D>::
MULTI_INDEX_BOX(const MULTI_INDEX_TYPE& multi_index)
    : min_multi_index(multi_index),
      max_multi_index(multi_index),
      widths(MULTI_INDEX_TYPE::All_Ones_Vector())
{ }

template< int D >
inline
MULTI_INDEX_BOX<D>::
MULTI_INDEX_BOX(
    const MULTI_INDEX_TYPE& min_multi_index_,
    const MULTI_INDEX_TYPE& max_multi_index_)
    : min_multi_index(min_multi_index_),
      max_multi_index(max_multi_index_),
      widths(1 + max_multi_index - min_multi_index)
{ }

template< int D >
inline typename MULTI_INDEX_BOX<D>::MULTI_INDEX_TYPE const &
MULTI_INDEX_BOX<D>::
Min_Multi_Index() const
{ return min_multi_index; }

template< int D >
inline typename MULTI_INDEX_BOX<D>::MULTI_INDEX_TYPE const &
MULTI_INDEX_BOX<D>::
Max_Multi_Index() const
{ return max_multi_index; }

template< int D >
inline int
MULTI_INDEX_BOX<D>::
Size() const
{ return widths.Product(); }

template< int D >
inline bool
MULTI_INDEX_BOX<D>::
Contains(const MULTI_INDEX_TYPE multi_index) const
{
    for(int d = 1; d <= D; ++d)
        if(multi_index[d] < min_multi_index[d] || max_multi_index[d] < multi_index[d])
            return false;
    return true;
}

namespace Detail_MULTI_INDEX_BOX
{

template< int D >
struct DISPATCH;

template<>
struct DISPATCH<1>
{
    static VECTOR<int,1> Strides(const MULTI_INDEX_BOX<1>& /*this_*/)
    { return VECTOR<int,1>(1); }
    static int Linear_Index(const MULTI_INDEX_BOX<1>& this_, const VECTOR<int,1>& multi_index)
    { return 1 + (multi_index[1] - this_.min_multi_index[1]); }
    static VECTOR<int,1> Multi_Index(const MULTI_INDEX_BOX<1>& this_, const int linear_index)
    { return VECTOR<int,1>(this_.min_multi_index[1] + (linear_index - 1)); }
};

template<>
struct DISPATCH<2>
{
    static VECTOR<int,2> Strides(const MULTI_INDEX_BOX<2>& this_)
    { return VECTOR<int,2>(this_.widths[2], 1); }
    static int Linear_Index(const MULTI_INDEX_BOX<2>& this_, const VECTOR<int,2>& multi_index)
    {
        return 1 +               (multi_index[2] - this_.min_multi_index[2]) +
               this_.widths[2] * (multi_index[1] - this_.min_multi_index[1]);
    }
    static VECTOR<int,2> Multi_Index(const MULTI_INDEX_BOX<2>& this_, const int linear_index)
    {
        return VECTOR<int,2>(
            this_.min_multi_index[1] + (linear_index - 1) / this_.widths[2],
            this_.min_multi_index[2] + (linear_index - 1) % this_.widths[2]
        );
    }
};

template<>
struct DISPATCH<3>
{
    static VECTOR<int,3> Strides(const MULTI_INDEX_BOX<3>& this_)
    { return VECTOR<int,3>(this_.widths[2] * this_.widths[3], this_.widths[3], 1); }
    static int Linear_Index(const MULTI_INDEX_BOX<3>& this_, const VECTOR<int,3>& multi_index)
    {
        return 1 +                (multi_index[3] - this_.min_multi_index[3]) +
               this_.widths[3] * ((multi_index[2] - this_.min_multi_index[2]) +
               this_.widths[2] *  (multi_index[1] - this_.min_multi_index[1]));
    }
    static VECTOR<int,3> Multi_Index(const MULTI_INDEX_BOX<3>& this_, const int linear_index)
    {
        return VECTOR<int,3>(
            this_.min_multi_index[1] + (linear_index - 1) / (this_.widths[3] * this_.widths[2]),
            this_.min_multi_index[2] + (linear_index - 1) / this_.widths[3] % this_.widths[2],
            this_.min_multi_index[3] + (linear_index - 1) % this_.widths[3]
        );
    }
};

} // namespace Detail_MULTI_INDEX_BOX

template< int D >
inline typename MULTI_INDEX_BOX<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOX<D>::
Strides() const
{ return Detail_MULTI_INDEX_BOX::DISPATCH<D>::Strides(*this); }

template< int D >
inline int
MULTI_INDEX_BOX<D>::
Linear_Index(const MULTI_INDEX_TYPE& multi_index) const
{
#ifndef NDEBUG
    for(int d = 1; d <= D; ++d)
        assert(min_multi_index[d] <= multi_index[d] && multi_index[d] <= max_multi_index[d]);
#endif // #ifndef NDEBUG
    return Detail_MULTI_INDEX_BOX::DISPATCH<D>::Linear_Index(*this, multi_index);
}

template< int D >
inline typename MULTI_INDEX_BOX<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOX<D>::
Multi_Index(const int linear_index) const
{
    assert(1 <= linear_index && linear_index <= widths.Product());
    return Detail_MULTI_INDEX_BOX::DISPATCH<D>::Multi_Index(*this, linear_index);
}

template< int D >
inline int
MULTI_INDEX_BOX<D>::
operator()(const MULTI_INDEX_TYPE& multi_index) const
{ return Linear_Index(multi_index); }

template< int D >
inline typename MULTI_INDEX_BOX<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOX<D>::
operator()(const int linear_index) const
{ return Multi_Index(linear_index); }

template< int D >
inline typename MULTI_INDEX_BOX<D>::iterator
MULTI_INDEX_BOX<D>::
begin() const
{ return iterator(*this, BEGIN_TAG()); }

template< int D >
inline typename MULTI_INDEX_BOX<D>::iterator
MULTI_INDEX_BOX<D>::
end() const
{ return iterator(*this, END_TAG()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_HPP
