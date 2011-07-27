//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//
// Formulas for 2D were taken from
//     http://en.wikipedia.org/wiki/Gaussian_quadrature
// Formulas for 3D were taken from
//     "Gaussian Quadrature Formulas for Triangles" by G.R. Cowper.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_GAUSSIAN_QUADRATURE_SIMPLEX_RULE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_GAUSSIAN_QUADRATURE_SIMPLEX_RULE_HPP

#include <boost/mpl/assert.hpp>

#include <Jeffrey_Utilities/Math/QUADRATURE_WEIGHT_X.h>

namespace PhysBAM
{

template< class T, int D, int ORDER >
struct GAUSSIAN_QUADRATURE_SIMPLEX_RULE;

#define PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( D, ORDER, N_WEIGHT_X ) \
template< class T > \
struct GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, D, ORDER > \
{ \
    static const int n_weight_x = N_WEIGHT_X; \
    static const QUADRATURE_WEIGHT_X<T,D> weight_x[N_WEIGHT_X]; \
};

//#####################################################################
// struct GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, 1, ORDER >
//#####################################################################

PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 1, 0, 1 )

template< class T, int ORDER >
struct GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, 1, ORDER >
    : GAUSSIAN_QUADRATURE_SIMPLEX_RULE<T,1,0>
{ };

//#####################################################################
// struct GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, 2, ORDER >
//#####################################################################

PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 2, 1, 1 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 2, 3, 2 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 2, 5, 3 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 2, 7, 4 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 2, 9, 5 )

template< class T, int ORDER >
struct GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, 2, ORDER >
    : GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, 2, ORDER + 1 >
{ BOOST_MPL_ASSERT_RELATION( ORDER, <, 9 ); };

//#####################################################################
// struct GAUSSIAN_QUADRATURE_SIMPLEX_RULE< T, 3, ORDER >
//#####################################################################

template< class T >
struct GAUSSIAN_QUADRATURE_SIMPLEX_RULE<T,3,0>
    : GAUSSIAN_QUADRATURE_SIMPLEX_RULE<T,3,1>
{ };

PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 3, 1, 1 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 3, 2, 3 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 3, 3, 4 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 3, 4, 6 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 3, 5, 7 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 3, 6, 12 )
PHYSBAM_GAUSSIAN_QUADRATURE_SIMPLEX_RULE( 3, 7, 13 )

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_GAUSSIAN_QUADRATURE_SIMPLEX_RULE_HPP
