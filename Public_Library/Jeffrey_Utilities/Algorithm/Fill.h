//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FILL_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FILL_HPP

#include <algorithm>

#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T_ARRAY, class T >
inline void
Fill(T_ARRAY& array, const T& x)
{
    typedef typename ARRAY_VALUE< T_ARRAY >::type VALUE_TYPE;
    BOOST_MPL_ASSERT_NOT((boost::is_const< VALUE_TYPE >));
    BOOST_MPL_ASSERT((boost::is_convertible< const T&, VALUE_TYPE >));
    std::fill(array.begin(), array.end(), x);
}

//#####################################################################
//#####################################################################

namespace Detail_Fill
{

template< class T, class U >
struct FILL_HELPER;

} // namespace Detail_Fill

template< class T_ARRAY, class T >
inline void
Fill_MT(
    const unsigned int n_thread,
    T_ARRAY& array, const T& x)
{
    typedef typename ARRAY_VALUE< T_ARRAY >::type VALUE_TYPE;
    BOOST_MPL_ASSERT_NOT((boost::is_const< VALUE_TYPE >));
    BOOST_MPL_ASSERT((boost::is_convertible< const T&, VALUE_TYPE >));
    typedef Detail_Fill::FILL_HELPER< VALUE_TYPE, T > FILL_HELPER_;
    assert(n_thread >= 1);
    For_Each_MT(
        n_thread,
        0, Size(array) - 1,
        FILL_HELPER_(&Front(array), x)
    );
}

namespace Detail_Fill
{

template< class T, class U >
struct FILL_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        FILL_HELPER,
        (( typename T* const, p_array ))
        (( typename U const, x ))
    )
public:
    typedef void result_type;
    void operator()(const int i) const
    { p_array[i] = x; }
};

} // namespace Detail_Fill

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FOR_EACH_HPP
