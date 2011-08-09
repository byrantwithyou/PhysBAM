//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_RAND_MT19937_UNIFORM_REAL_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_RAND_MT19937_UNIFORM_REAL_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T >
struct RAND_MT19937_UNIFORM_REAL
{
    typedef boost::variate_generator<
        boost::mt19937&,
        boost::uniform_real<T>
    > type;
};

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_RAND_MT19937_UNIFORM_REAL_HPP
