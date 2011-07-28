//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_BETA_GRAD_U_DOT_N_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_BETA_GRAD_U_DOT_N_HPP

#include <boost/type_traits/remove_reference.hpp>
#include <boost/utility/result_of.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T_BETA, class T_GRAD_U >
struct BETA_GRAD_U_DOT_N
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        BETA_GRAD_U_DOT_N,
        (( typename T_BETA const &, beta ))
        (( typename T_GRAD_U const &, grad_u ))
    )
public:
    template<class> struct result;
    template< class T_THIS, class T_X, class T_N >
    struct result< T_THIS ( T_X, T_N ) >
    { typedef typename boost::remove_reference< T_X >::type::value_type type; };

    template< class T, int D >
    T operator()(const VECTOR<T,D>& x, const VECTOR<T,D>& n) const
    { return beta(x) * VECTOR<T,D>::Dot_Product(grad_u(x), n); }
};

template< class T_BETA, class T_GRAD_U >
inline BETA_GRAD_U_DOT_N< T_BETA, T_GRAD_U >
Make_Beta_Grad_U_Dot_N(const T_BETA& beta, const T_GRAD_U& grad_u)
{ return BETA_GRAD_U_DOT_N< T_BETA, T_GRAD_U >(beta, grad_u); }

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_BETA_GRAD_U_DOT_N_HPP
