//#####################################################################
// Copyright 2008, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_REGION
//##################################################################### 
#ifndef __FRACTURE_REGION__
#define __FRACTURE_REGION__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Tools/Read_Write/READ_WRITE_FORWARD.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{

template<class T>
class FRACTURE_REGION
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
public:
    typedef int HAS_TYPED_READ_WRITE;

    TRIANGULATED_SURFACE<T>* triangulated_surface;
    LEVELSET_IMPLICIT_OBJECT<TV>* implicit_object;
    MATRIX<T,3> levelset_RS;TV levelset_T;
    MATRIX<T,3> object_RS;TV object_T;
    TV_INT fracture_offset;
    T particle_intersection_thickness;
    bool need_destroy_data;
    FRAME<TV> extra_levelset_frame;
    PARTICLE_PARTITION<TV>* particle_partition;
    
    FRACTURE_REGION(TRIANGULATED_SURFACE<T>* triangulated_surface_input,LEVELSET_IMPLICIT_OBJECT<TV>* implicit_object_input,const bool create_particle_partition);
    virtual ~FRACTURE_REGION();

//#####################################################################
    ARRAY<FRACTURE_REGION<T>*> Intersect_With_Rigid_Body(const FRACTURE_REGION<T>& body,const bool use_particle_optimization,const bool tessellate_region=false);
    T Compute_Volume() const;
    void Compute_Inertial_Properties(const T density,TV& com,T& mass,SYMMETRIC_MATRIX<T,TV::SPIN::m>& inertia) const; 
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
    void Initialize_Particle_Partition();
//#####################################################################
};
}
#endif
