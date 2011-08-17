//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_MAXIMUM_MAGNITUDE_CONVERGENCE_NORM_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_MAXIMUM_MAGNITUDE_CONVERGENCE_NORM_HPP

#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>

namespace PhysBAM
{

template< class T >
struct MAXIMUM_MAGNITUDE_CONVERGENCE_NORM
{
    typedef T result_type;
    template< class T_VECTOR >
    T operator()(const T_VECTOR& x) const
    { return ARRAYS_COMPUTATIONS::Maximum_Magnitude(x); }
};


} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_MAXIMUM_MAGNITUDE_CONVERGENCE_NORM_HPP
