//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_CONSTANT_VECTORS_IN_NULL_SPACE_H
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_CONSTANT_VECTORS_IN_NULL_SPACE_H

#include <cmath>

#include <limits>

#include <Jeffrey_Utilities/Algorithm/Any_If.h>

namespace PhysBAM
{

namespace Detail_Has_Constant_Vectors_In_Null_Space
{

template< class T, class T_SYSTEM >
struct NONZERO_STENCIL_SUM
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        NONZERO_STENCIL_SUM,
        (( typename T_SYSTEM const &, system ))
        (( /******/ int const, threshold_factor ))
    )
public:
    typedef bool result_type;
    bool operator()(const int index) const
    {
        const T eps = std::numeric_limits<T>::epsilon();
        const T diag = system.Diag(index);
        const int n_nonzero = system.Stencil_N_Nonzero(index);
        const T stencil_sum = system.Stencil_Sum(index);
        return std::abs(stencil_sum) > threshold_factor * n_nonzero * eps * std::abs(diag);
    }
};

} // namespace Detail_Has_Constant_Vectors_In_Null_Space

template< class T, class T_SYSTEM >
inline bool
Has_Constant_Vectors_In_Null_Space(
    const unsigned int n_thread,
    const int n_index,
    const T_SYSTEM& system,
    const int threshold_factor = 1024)
{
    typedef Detail_Has_Constant_Vectors_In_Null_Space::NONZERO_STENCIL_SUM< T, T_SYSTEM > NONZERO_STENCIL_SUM_;
    return !Any_If_MT(
        n_thread,
        1, n_index,
        NONZERO_STENCIL_SUM_(system, threshold_factor)
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_CONSTANT_VECTORS_IN_NULL_SPACE_H
