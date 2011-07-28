//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PRE_COMPOSE_TRANSLATE_X_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PRE_COMPOSE_TRANSLATE_X_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D, class T_F_OF_X_AND_N >
struct PRE_COMPOSE_TRANSLATE_X
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        PRE_COMPOSE_TRANSLATE_X,
        (( typename T_F_OF_X_AND_N const, f_of_x_and_n ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, x0 ))
    )
public:
    typedef T result_type;
    T operator()(const VECTOR<T,D>& x, const VECTOR<T,D>& n) const
    { return f_of_x_and_n(x + x0, n); }
};

template< class T_F_OF_X_AND_N, class T, int D >
inline PRE_COMPOSE_TRANSLATE_X< T, D, T_F_OF_X_AND_N >
Make_Pre_Compose_Translate_X(const T_F_OF_X_AND_N& f_of_x_and_n, const VECTOR<T,D>& x0)
{ return PRE_COMPOSE_TRANSLATE_X< T, D, T_F_OF_X_AND_N >(f_of_x_and_n, x0); }

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PRE_COMPOSE_TRANSLATE_X_HPP
