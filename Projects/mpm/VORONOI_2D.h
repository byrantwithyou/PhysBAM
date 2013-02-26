//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORONOI_2D
//#####################################################################
#ifndef __VORONOI_2D__
#define __VORONOI_2D__
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
namespace PhysBAM{
template<class T>
class VORONOI_2D
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
public:    
    GEOMETRY_PARTICLES<TV> segment_mesh_particles;
    ARRAY<TV> Xm;
    SEGMENT_MESH segment_mesh;
    ARRAY<ARRAY<TV_INT> > local_voronoi_elements_of_sample_particle;
    ARRAY<ARRAY<int> > local_voronoi_particles_of_sample_particle;

    VORONOI_2D(){}
    ~VORONOI_2D(){}

    void Initialize_With_A_Regular_Grid(const GRID<TV>& grid);
//#####################################################################
};
}
#endif
