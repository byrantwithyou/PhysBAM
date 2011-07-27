//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_QUADRATURE_WEIGHT_X_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_QUADRATURE_WEIGHT_X_HPP

namespace PhysBAM
{

template< class T, int D >
struct QUADRATURE_WEIGHT_X
{
    T weight;
    T x[D];
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_QUADRATURE_WEIGHT_X_HPP
