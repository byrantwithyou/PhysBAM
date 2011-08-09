//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_MAIN_EXAMPLES_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_MAIN_EXAMPLES_HPP

#include <exception>

#include <boost/exception/exception.hpp>
#include <boost/function.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Params/EXAMPLE_PARAMS.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

struct U_BETA_ERROR
    : virtual std::exception, virtual boost::exception
{
    const char* what() const throw ( )
    { return "PhysBAM::Multigrid_Embedded_Poisson::U_BETA_ERROR"; }
};

template< int D >
void Visit_U_Examples_Str(const boost::function< void ( unsigned int, const char* ) >& visitor);
template< int D >
void Visit_Beta_Examples_Str(const boost::function< void ( unsigned int, const char* ) >& visitor);

template< class T, int D >
void
Get_Main_Example(
    const unsigned int u_id, const unsigned int beta_id,
    typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& u,
    typename EXAMPLE_PARAMS<T,D>::VECTOR_FUNCTION_TYPE& grad_u,
    typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& beta,
    typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& f);

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_MAIN_EXAMPLES_HPP
