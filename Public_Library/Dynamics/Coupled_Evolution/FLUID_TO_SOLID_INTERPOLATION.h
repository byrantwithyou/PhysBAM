//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION
//#####################################################################
#ifndef __FLUID_TO_SOLID_INTERPOLATION__
#define __FLUID_TO_SOLID_INTERPOLATION__
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_BASE.h>

namespace PhysBAM{

template<class TV>
class FLUID_TO_SOLID_INTERPOLATION:public FLUID_TO_SOLID_INTERPOLATION_BASE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    using FLUID_TO_SOLID_INTERPOLATION_BASE<TV>::index_map;

    struct ENTRY
    {
        T w;
        int i;
    };

public:
    T max_dist;
    ARRAY<VECTOR<ARRAY<ENTRY>,TV::m> > entries;
    ARRAY<int> coupled_particles;
    const DEFORMABLE_PARTICLES<TV>& particles;

    FLUID_TO_SOLID_INTERPOLATION(const COLLISION_AWARE_INDEX_MAP<TV>& map,const DEFORMABLE_PARTICLES<TV>& particles_input);
    virtual ~FLUID_TO_SOLID_INTERPOLATION();

//#####################################################################
    void Compute(const int ghost_cells) override;
    void Compute_Weights(const TV& X,int axis,ARRAY<ENTRY>& array);
    void Times_Add(const ARRAY<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const override;
    void Transpose_Times_Add(const GENERALIZED_VELOCITY<TV>& solid_force,ARRAY<T>& fluid_force) const override;
    void Print_Each_Matrix(int n,int fluid_faces,GENERALIZED_VELOCITY<TV>& G) const override;
//#####################################################################
};
}
#endif
