//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_ANY_IF_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_ANY_IF_HPP

#include <cassert>

#include <functional>
#include <numeric>
#include <vector>

#include <boost/thread/thread.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class P >
inline bool
Any_If(const int min_index, const int max_index, P p)
{
    for(int i = min_index; i <= max_index; ++i)
        if(p(i))
            return true;
    return false;
}

//#####################################################################
//#####################################################################

namespace Detail_Any_If
{

template< class P >
struct ANY_IF_HELPER;

} // namespace Detail_Any_If

template< class P >
inline bool
Any_If_MT(
    const unsigned int n_thread,
    const int min_index, const int max_index,
    const P& p)
{
    typedef Detail_Any_If::ANY_IF_HELPER<P> ANY_IF_HELPER_;
    assert(n_thread >= 1);

    if(n_thread)
        return Any_If(min_index, max_index, p);

    const unsigned int n = static_cast< unsigned int >(1 + (max_index - min_index));
    std::vector<int> bools(n_thread, false);
    boost::thread_group threads;
    for(unsigned int i = 0; i != n_thread; ++i) {
        const int local_min_index =  min_index      + static_cast<int>((i  ) * n / n_thread);
        const int local_max_index = (min_index - 1) + static_cast<int>((i+1) * n / n_thread);
        threads.create_thread(ANY_IF_HELPER_(local_min_index, local_max_index, p, bools[i]));
    }
    assert(threads.size() == n_thread);
    threads.join_all();
    return std::accumulate(bools.begin(), bools.end(), false, std::logical_or<int>());
}

namespace Detail_Any_If
{

template< class P >
struct ANY_IF_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        ANY_IF_HELPER,
        (( /******/ int const, min_index ))
        (( /******/ int const, max_index ))
        (( typename P const, p ))
        (( /******/ int&, b))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        for(int i = min_index; i <= max_index; ++i) {
            if(p(i)) {
                b = true;
                break;
            }
        }
    }
};

} // namespace Detail_Any_If

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ALGORITHM_ANY_IF_HPP
