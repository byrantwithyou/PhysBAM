//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_FORCE_HELPER
//#####################################################################
#ifndef __MPM_FORCE_HELPER__
#define __MPM_FORCE_HELPER__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Hybrid_Methods/Forces/OLDROYD_CONSTITUTIVE_MODEL.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
namespace PhysBAM{

template<class TV>
class MPM_FORCE_HELPER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_PARTICLES<TV>& particles;
    ARRAY<MATRIX<T,TV::m> > Fn,A;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > Sn;

    MPM_FORCE_HELPER(MPM_PARTICLES<TV>& particles)
        :particles(particles)
    {
    }
};
}
#endif
