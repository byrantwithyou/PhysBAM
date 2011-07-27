//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FOR_EACH_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FOR_EACH_HPP

#include <cassert>

#include <boost/thread/thread.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class F >
inline void
For_Each(const int min_index, const int max_index, F f)
{
    for(int i = min_index; i <= max_index; ++i)
        f(i);
}

//#####################################################################
//#####################################################################

namespace Detail_For_Each
{

template< class F >
struct FOR_EACH_HELPER;

} // namespace Detail_For_Each

template< class F >
inline void
For_Each_MT(
    const unsigned int n_thread,
    const int min_index, const int max_index,
    const F& f)
{
    typedef Detail_For_Each::FOR_EACH_HELPER<F> FOR_EACH_HELPER_;
    assert(n_thread >= 1);
    const unsigned int n = static_cast< unsigned int >(1 + (max_index - min_index));
    boost::thread_group threads;
    for(unsigned int i = 0; i != n_thread; ++i) {
        const int local_min_index =  min_index      + static_cast<int>((i  ) * n / n_thread);
        const int local_max_index = (min_index - 1) + static_cast<int>((i+1) * n / n_thread);
        threads.create_thread(FOR_EACH_HELPER_(local_min_index, local_max_index, f));
    }
    assert(threads.size() == n_thread);
    threads.join_all();
}

namespace Detail_For_Each
{

template< class F >
struct FOR_EACH_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        FOR_EACH_HELPER,
        (( /******/ int const, min_index ))
        (( /******/ int const, max_index ))
        (( typename F const, f ))
    )
public:
    typedef void result_type;
    void operator()()
    {
        for(int i = min_index; i <= max_index; ++i)
            f(i);
    }
};

} // namespace Detail_For_Each

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_FOR_EACH_HPP
