//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FIND_ALL_IF_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FIND_ALL_IF_HPP

#include <numeric>
#include <vector>

#include <cassert>

#include <boost/mpl/assert.hpp>
#include <boost/thread/thread.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include <Jeffrey_Utilities/Algorithm/Count_If.h>
#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/IDENTITY_FUNCTION.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

template< class T_P, class T_F, class T_ARRAY >
inline void
Find_All_If(
    const int min_index, const int max_index, T_P p, T_F f,
    T_ARRAY& result)
{
    typedef typename ARRAY_VALUE< T_ARRAY >::type VALUE_TYPE;
    BOOST_MPL_ASSERT_NOT((boost::is_const< VALUE_TYPE >));
    BOOST_MPL_ASSERT((boost::is_convertible<
        typename RESULT_OF< T_F ( int ) >::type,
        VALUE_TYPE
    >));
    const int count = Count_If(min_index, max_index, p);
    if(count == 0)
        return;
    const int old_size = static_cast< int >(Size(result));
    Exact_Resize(result, old_size + count, false, true);
    VALUE_TYPE* p_result = &Front(result) + old_size;
    for(int i = min_index; i <= max_index; ++i) {
        if(p(i)) {
            *p_result = f(i);
            ++p_result;
        }
    }
}

template< class P, class T_ARRAY >
inline void
Find_All_If(
    const int min_index, const int max_index, const P& p,
    T_ARRAY& result)
{ Find_All_If(min_index, max_index, p, IDENTITY_FUNCTION(), result); }

//#####################################################################
//#####################################################################

namespace Detail_Find_All_If
{

template< class P, class F, class T >
struct FIND_ALL_IF_HELPER;

} // namespace Detail_Find_All_If

template< class T_P, class T_F, class T_ARRAY >
inline void
Find_All_If_MT(
    const unsigned int n_thread,
    const int min_index, const int max_index, const T_P& p, const T_F& f,
    T_ARRAY& result)
{
    typedef typename ARRAY_VALUE< T_ARRAY >::type VALUE_TYPE;
    BOOST_MPL_ASSERT_NOT((boost::is_const< VALUE_TYPE >));
    BOOST_MPL_ASSERT((boost::is_convertible<
        typename RESULT_OF< T_F ( int ) >::type,
        VALUE_TYPE
    >));
    typedef Detail_Count_If::COUNT_IF_HELPER<T_P> COUNT_IF_HELPER_;
    typedef Detail_Find_All_If::FIND_ALL_IF_HELPER< T_P, T_F, VALUE_TYPE > FIND_ALL_IF_HELPER_;

    const unsigned int n = static_cast< unsigned int >(1 + (max_index - min_index));

    std::vector<int> counts(n_thread, 0); 
    {
        boost::thread_group threads;
        for(unsigned int i = 0; i != n_thread; ++i) {
            const int local_min_index =  min_index      + static_cast<int>((i  ) * n / n_thread);
            const int local_max_index = (min_index - 1) + static_cast<int>((i+1) * n / n_thread);
            threads.create_thread(COUNT_IF_HELPER_(local_min_index, local_max_index, p, counts[i]));
        }
        assert(threads.size() == n_thread);
        threads.join_all();
    }

    const int count = std::accumulate(counts.begin(), counts.end(), 0);
    if(count == 0)
        return;
    const int old_size = static_cast< int >(Size(result));
    Exact_Resize(result, old_size + count, false, true);
    VALUE_TYPE* p_result = &Front(result) + old_size;
    {
        boost::thread_group threads;
        for(unsigned int i = 0; i != n_thread; ++i) {
            const int local_min_index =  min_index      + static_cast<int>((i  ) * n / n_thread);
            const int local_max_index = (min_index - 1) + static_cast<int>((i+1) * n / n_thread);
            threads.create_thread(FIND_ALL_IF_HELPER_(local_min_index, local_max_index, p, f, p_result));
            p_result += counts[i];
        }
        assert(threads.size() == n_thread);
        threads.join_all();
    }
}

template< class P, class T_ARRAY >
inline void
Find_All_If_MT(
    const unsigned int n_thread,
    const int min_index, const int max_index, const P& p,
    T_ARRAY& result)
{ Find_All_If_MT(n_thread, min_index, max_index, p, IDENTITY_FUNCTION(), result); }

namespace Detail_Find_All_If
{

template< class T_P, class T_F, class T >
struct FIND_ALL_IF_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        FIND_ALL_IF_HELPER,
        (( /******/ int const, min_index ))
        (( /******/ int const, max_index ))
        (( typename T_P const, p ))
        (( typename T_F const, f ))
        (( typename T* const, p_result ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        T* p_result_ = p_result;
        for(int i = min_index; i <= max_index; ++i) {
            if(p(i)) {
                *p_result_ = f(i);
                ++p_result_;
            }
        }
    }
};

} // namespace Detail_Find_All_If

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FIND_ALL_IF_HPP
