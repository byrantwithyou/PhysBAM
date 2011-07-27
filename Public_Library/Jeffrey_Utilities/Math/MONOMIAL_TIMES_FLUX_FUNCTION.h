//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_TIMES_FLUX_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_TIMES_FLUX_FUNCTION_HPP

#include <boost/type_traits/remove_reference.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T_MONOMIAL, class T_F >
struct MONOMIAL_TIMES_FLUX_FUNCTION
{
    static const int DIMENSION = T_MONOMIAL::DIMENSION;
    typedef typename T_MONOMIAL::MULTI_POWER_TYPE MULTI_POWER_TYPE;

    explicit MONOMIAL_TIMES_FLUX_FUNCTION(const T_F& f)
        : m_f(f), m_monomial()
    { }
    MONOMIAL_TIMES_FLUX_FUNCTION(const T_F& f, const T_MONOMIAL& monomial)
        : m_f(f), m_monomial(monomial)
    { }

    template<class> struct result;
    template< class T_THIS, class T_VECTOR1, class T_VECTOR2 >
    struct result< T_THIS ( T_VECTOR1, T_VECTOR2 ) >
    { typedef typename boost::remove_reference< T_VECTOR1 >::type::value_type type; };

    template< class T >
    T operator()(
        const VECTOR< T, DIMENSION >& x,
        const VECTOR< T, DIMENSION >& normal) const
    { return m_monomial(x) * m_f(x, normal); }

private:
    const T_F m_f;
    const T_MONOMIAL m_monomial;
};

template< class T_MONOMIAL, class T_F >
inline MONOMIAL_TIMES_FLUX_FUNCTION< T_MONOMIAL, T_F >
Make_Monomial_Times_Flux_Function(const T_F& f)
{ return MONOMIAL_TIMES_FLUX_FUNCTION< T_MONOMIAL, T_F >(f); }

template< class T_MONOMIAL, class T_F >
inline MONOMIAL_TIMES_FLUX_FUNCTION< T_MONOMIAL, T_F >
Make_Monomial_Times_Flux_Function(const T_MONOMIAL& monomial, const T_F& f)
{ return MONOMIAL_TIMES_FLUX_FUNCTION< T_MONOMIAL, T_F >(monomial, f); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MONOMIAL_TIMES_FLUX_FUNCTION_HPP
