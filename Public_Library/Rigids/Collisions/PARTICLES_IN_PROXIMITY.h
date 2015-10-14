//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace PARTICLES_IN_PROXIMITY
//##################################################################### 
#ifndef __PARTICLES_IN_PROXIMITY__
#define __PARTICLES_IN_PROXIMITY__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAY_VIEW.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/UNION_FIND.h>
#include <Tools/Vectors/VECTOR_1D.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>

namespace PhysBAM{

namespace PARTICLES_IN_PROXIMITY
{
template<class T> ARRAY<int> Convex_Hull(ARRAY<VECTOR<T,2> >& points);
template<class TV> void Particles_In_Proximity(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,
    ARRAY<TV>& locations,ARRAY<TV>& normals,ARRAY<typename TV::SCALAR>& distances,typename TV::SCALAR proximity_distance);
template<class TV> PARTICLE_PARTITION<TV> Build_Particle_Partition(ARRAY<TV>& locations,
    typename TV::SCALAR region_threshold);
template<class TV> void Aggregate_And_Stagger_Convex_Regions(ARRAY<TV>& locations,ARRAY<TV>& normals,
    ARRAY<typename TV::SCALAR>& distances,typename TV::SCALAR region_threshold);
template<class TV>
void All_Particles_In_Proximity(RIGID_BODY<TV>& body_1,RIGID_BODY<TV>& body_2,ARRAY<TV>& locations,
    ARRAY<TV>& normals,ARRAY<typename TV::SCALAR>& distances,typename TV::SCALAR proximity_distance,
    const bool stagger_points=true);
}
}

#endif
