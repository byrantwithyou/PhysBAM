//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_COUNT_IF_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_COUNT_IF_HPP

#include <cassert>

#include <numeric>
#include <vector>

#include <boost/thread/thread.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class P >
inline int
Count_If(const int min_index, const int max_index, P p)
{
    int count = 0;
    for(int i = min_index; i <= max_index; ++i)
        if(p(i))
            ++count;
    return count;
}

//#####################################################################
//#####################################################################

namespace Detail_Count_If
{

template< class P >
struct COUNT_IF_HELPER;

} // namespace Detail_Count_If

template< class P >
inline int
Count_If_MT(
    const unsigned int n_thread,
    const int min_index, const int max_index,
    const P& p)
{
    typedef Detail_Count_If::COUNT_IF_HELPER<P> COUNT_IF_HELPER_;
    assert(n_thread >= 1);

    if(n_thread == 1)
        return Count_If(min_index, max_index, p);

    const unsigned int n = static_cast< unsigned int >(1 + (max_index - min_index));
    std::vector<int> counts(n_thread, 0);
    boost::thread_group threads;
    for(unsigned int i = 0; i != n_thread; ++i) {
        const int local_min_index =  min_index      + static_cast<int>((i  ) * n / n_thread);
        const int local_max_index = (min_index - 1) + static_cast<int>((i+1) * n / n_thread);
        threads.create_thread(COUNT_IF_HELPER_(local_min_index, local_max_index, p, counts[i]));
    }
    assert(threads.size() == n_thread);
    threads.join_all();
    return std::accumulate(counts.begin(), counts.end(), 0);
}

namespace Detail_Count_If
{

template< class P >
struct COUNT_IF_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        COUNT_IF_HELPER,
        (( /******/ int const, min_index ))
        (( /******/ int const, max_index ))
        (( typename P const, p ))
        (( /******/ int&, count ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        int count_ = 0;
        for(int i = min_index; i <= max_index; ++i)
            if(p(i))
                ++count_;
        count = count_;
    }
};

} // namespace Detail_Count_If

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_COUNT_IF_HPP
