//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORONOI_2D
//#####################################################################
#ifndef __VORONOI_2D__
#define __VORONOI_2D__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "MPM_PARTICLES.h"
namespace PhysBAM{
template<class T>
class VORONOI_2D
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
public:
    ARRAY<int> type; // 1-boundary, 10-interior, 100-face
    ARRAY<TV> Xm;
    ARRAY<ARRAY<int> > elements;
    ARRAY<TV> X;
    ARRAY<ARRAY<int> > association;
    ARRAY<TV_INT> segments;
    ARRAY<TV_INT> boundary_segments;
    ARRAY<TRIPLE<int,int,bool> > neighbor_cells;

    VORONOI_2D(){}
    ~VORONOI_2D(){}

    void Initialize_With_A_Regular_Grid_Of_Particles(const GRID<TV>& grid);
    void Initialize_With_A_Triangulated_Area(const TRIANGULATED_AREA<T>& ta);
    void Initialize_With_And_As_A_Triangulated_Area_And_Relocate_Particles_To_Tri_Centers(const TRIANGULATED_AREA<T>& ta,MPM_PARTICLES<TV>& mpm_particles);
    void Initialize_Neighbor_Cells();
    void Build_Association();
    void Build_Segments();
    void Build_Boundary_Segments();
    void Deform_Mesh_Using_Particle_Deformation(const ARRAY_VIEW<TV>& particle_Xm,const ARRAY_VIEW<TV>& particle_X,const ARRAY_VIEW<MATRIX<T,TV::m> >& particle_Fe,const ARRAY_VIEW<MATRIX<T,TV::m> >& particle_Fp,const bool constrain_face_centers=true,const int deformation_averaging_steps=0);
    void Crack(const ARRAY_VIEW<TV>& particle_X,const T threshold);
//#####################################################################
};
}
#endif
