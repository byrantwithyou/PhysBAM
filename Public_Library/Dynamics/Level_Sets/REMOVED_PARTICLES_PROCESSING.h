//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REMOVED_PARTICLES_PROCESSING
//#####################################################################
#ifndef __REMOVED_PARTICLES_PROCESSING__
#define __REMOVED_PARTICLES_PROCESSING__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Log/LOG.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Geometry/Basic_Geometry/ELLIPSOID.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
namespace PhysBAM{

template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class T> class REMOVED_PARTICLES_BLENDER_3D;

template<class T>
class REMOVED_PARTICLES_PROCESSING:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> particles;
    ARRAY<ELLIPSOID<T> > ellipsoids;
    ARRAY<SYMMETRIC_MATRIX<T,3> > metrics;
    REMOVED_PARTICLES_BLENDER_3D<T>* particle_blender;
    T blending_parameter;
    T scale;
    T relative_tolerance,tolerance;
    RANGE<TV> particle_domain;
    int grid_divisions;
    GRID<TV> particle_grid;
    ARRAY<ARRAY<int> ,VECTOR<int,3> > particle_array;
    KD_TREE<TV> particle_tree;
    
    REMOVED_PARTICLES_PROCESSING(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_input);
    REMOVED_PARTICLES_PROCESSING(GRID<TV>& grid_input,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> *,VECTOR<int,3> >& particles_array_input);

private:
    void Initialize(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_input);
    void Initialize(const GRID<TV>& grid,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> >& particles_array);
//#####################################################################
public:
    void Setup_Processing();
    T Phi(const TV& position) const;
    TV Normal(const TV& position) const;
private:
    ELLIPSOID<T> Get_Ellipsoid(const int p) const;
//#####################################################################
};
}
#endif
