//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_CURVATURE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_CURVATURE_HPP

#include <cmath>

#include <boost/math/special_functions/pow.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D, class T_PHI_OF_INDEX >
inline T Curvature(
    const VECTOR<T,D> dx,
    const VECTOR<int,D> strides,
    const T max_curvature,
    T_PHI_OF_INDEX phi_of_index,
    const int linear_index)
{
    VECTOR<T,D> d1phi; // init'ed to 0
    VECTOR<T,D> d2phi; // init'ed to 0
    const T phi11 = phi_of_index(linear_index);
    for(int d1 = 1; d1 <= D; ++d1) {
        const T one_over_dx1 = 1 / dx[d1];
        const T phi01 = phi_of_index(linear_index - strides[d1]);
        const T phi21 = phi_of_index(linear_index + strides[d1]);
        d1phi[d1] = (phi21 - phi01) * one_over_dx1 / 2;
        d2phi[d1] = ((phi21 - phi11) - (phi11 - phi01)) * (one_over_dx1 * one_over_dx1);
    }
    const T denominator = d1phi.Magnitude();
    if(denominator != 0) {
        T result = 0;
        for(int d1 = 1; d1 < D; ++d1) {
            const T one_over_dx1 = 1 / dx[d1];
            for(int d2 = d1 + 1; d2 <= D; ++d2) {
                const T one_over_dx2 = 1 / dx[d2];
                const T phi00 = phi_of_index(linear_index - strides[d1] - strides[d2]);
                const T phi02 = phi_of_index(linear_index - strides[d1] + strides[d2]);
                const T phi20 = phi_of_index(linear_index + strides[d1] - strides[d2]);
                const T phi22 = phi_of_index(linear_index + strides[d1] + strides[d2]);
                const T d12phi = ((phi22 - phi20) - (phi02 - phi00)) * (one_over_dx1 * one_over_dx2) / 4;
                result += d1phi[d1] * (d1phi[d2] * d12phi - d1phi[d1] * d2phi[d2])
                       +  d1phi[d2] * (d1phi[d1] * d12phi - d1phi[d2] * d2phi[d1]);
            }
        }
        result /= boost::math::pow<3>(denominator);
        if(std::abs(result) <= max_curvature)
            return result;
    }
    return (phi11 <= 0 ? -1 : +1) * max_curvature;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_CURVATURE_HPP
