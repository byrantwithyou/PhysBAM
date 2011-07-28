//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_MIN_DIST_TO_PARAMETRIC_CURVE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_MIN_DIST_TO_PARAMETRIC_CURVE_HPP

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <limits>

#include <boost/math/tools/minima.hpp>
#include <boost/mpl/assert.hpp>

#include <Jeffrey_Utilities/Functional/Post_Compose_Magnitude_Squared.h>
#include <Jeffrey_Utilities/Functional/Post_Compose_Translate.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int N_SUBINTERVAL, class T, int D, class T_PARAMETRIC_CURVE >
inline T
Min_Dist_To_Parametric_Curve_Via_Subintervals(
    const VECTOR<T,D> x,
    const T_PARAMETRIC_CURVE& parametric_curve)
{
    T min_dist2 = std::numeric_limits<T>::infinity();
    for(int i = 0; i != N_SUBINTERVAL; ++i) {
        const T dist2 = boost::math::tools::brent_find_minima(
            Post_Compose_Magnitude_Squared(Post_Compose_Translate(-x, parametric_curve)),
            static_cast<T>(i    ) / N_SUBINTERVAL,
            static_cast<T>(i + 1) / N_SUBINTERVAL,
            std::numeric_limits<T>::digits / 2
        ).second;
        min_dist2 = std::min(min_dist2, dist2);
    }
    return std::sqrt(min_dist2);
}

template< int N_LOCAL_MINIMA, class T, int D, class T_PARAMETRIC_CURVE, int N_CACHED_X >
inline T
Min_Dist_To_Parametric_Curve_Via_Cached_X(
    const VECTOR<T,D> x,
    const T_PARAMETRIC_CURVE& parametric_curve,
    const VECTOR<T,D> (&cached_x)[N_CACHED_X])
{
    BOOST_MPL_ASSERT_RELATION( N_CACHED_X, >=, N_LOCAL_MINIMA );
    std::pair<T,int> min_dist2_i[N_LOCAL_MINIMA];
    std::fill(
        &min_dist2_i[0], &min_dist2_i[N_LOCAL_MINIMA],
        std::pair<T,int>(std::numeric_limits<T>::infinity(), 0)
    );
    for(int i = 0; i != N_CACHED_X; ++i) {
        const T dist2 = (cached_x[i] - x).Magnitude_Squared();
        if(dist2 < min_dist2_i[N_LOCAL_MINIMA-1].first) {
            min_dist2_i[N_LOCAL_MINIMA-1] = std::pair<T,int>(dist2, i);
            std::inplace_merge(
                &min_dist2_i[0],
                &min_dist2_i[N_LOCAL_MINIMA-1],
                &min_dist2_i[N_LOCAL_MINIMA]
            );
        }
    }
    T min_dist2 = std::numeric_limits<T>::infinity();
    for(int j = 0; j != N_LOCAL_MINIMA; ++j) {
        const int min_i = min_dist2_i[j].second;
        for(int k = 0; k != j; ++k)
            if(std::abs(min_dist2_i[k].second - min_i) <= 1)
                goto CONTINUE_FOR_J;
        {
            const T local_min_dist2 = boost::math::tools::brent_find_minima(
                Post_Compose_Magnitude_Squared(Post_Compose_Translate(-x, parametric_curve)),
                static_cast<T>(min_i - 1) / N_CACHED_X,
                static_cast<T>(min_i + 1) / N_CACHED_X,
                std::numeric_limits<T>::digits / 2
            ).second;
            min_dist2 = std::min(min_dist2, local_min_dist2);
        }
        CONTINUE_FOR_J:;
    }
    return std::sqrt(min_dist2);
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_MIN_DIST_TO_PARAMETRIC_CURVE_HPP
