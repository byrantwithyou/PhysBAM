//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_DOT_PRODUCT_INNER_PRODUCT_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_DOT_PRODUCT_INNER_PRODUCT_HPP

#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>

namespace PhysBAM
{

struct DOT_PRODUCT_INNER_PRODUCT
{
    typedef double result_type;
    template< class T_VECTOR >
    double operator()(const T_VECTOR& x, const T_VECTOR& y) const
    { return ARRAYS_COMPUTATIONS::Dot_Product_Double_Precision(x, y); }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_DOT_PRODUCT_INNER_PRODUCT_HPP
