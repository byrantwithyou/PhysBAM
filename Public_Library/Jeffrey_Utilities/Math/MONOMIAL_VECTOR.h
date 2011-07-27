//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_VECTOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_VECTOR_HPP

#include <boost/type_traits/remove_reference.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T_MONOMIAL, int AXIS >
struct MONOMIAL_VECTOR
{
    static const int DIMENSION = T_MONOMIAL::DIMENSION;
    typedef typename T_MONOMIAL::MULTI_POWER_TYPE MULTI_POWER_TYPE;

    MONOMIAL_VECTOR() : m_monomial() { }
    explicit MONOMIAL_VECTOR(const T_MONOMIAL& monomial) : m_monomial(monomial) { }

    template<class> struct result;
    template< class T_THIS, class T_VECTOR1, class T_VECTOR2 >
    struct result< T_THIS ( T_VECTOR1, T_VECTOR2 ) >
    { typedef typename boost::remove_reference< T_VECTOR1 >::type::value_type type; };

    template< class T >
    T operator()(
        const VECTOR< T, DIMENSION >& x,
        const VECTOR< T, DIMENSION >& normal) const
    { return normal[AXIS] * m_monomial(x); }

private:
    const T_MONOMIAL m_monomial;
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_VECTOR_HPP
