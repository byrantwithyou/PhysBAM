//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EDGE_EDGE_INITIAL_CULL_VISITOR
//##################################################################### 
#include <Tools/Arrays/ARRAY_VIEW.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <Deformables/Collisions_And_Interactions/EDGE_EDGE_INITIAL_CULL_VISITOR.h>
#include <Deformables/Forces/SEGMENT_ADHESION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> void EDGE_EDGE_INITIAL_CULL_VISITOR<TV>::
Store(const int segment1_local_index,const int segment2_local_index)
{
    int segment1_index=mesh1_indices(segment1_local_index),segment2_index=mesh2_indices(segment2_local_index);
    VECTOR<int,2> segment_indices=VECTOR<int,2>(segment1_index,segment2_index).Sorted();
    const VECTOR<int,2> &segment1_nodes=adhesion.curve.mesh.elements(segment_indices[0]),&segment2_nodes=adhesion.curve.mesh.elements(segment_indices[1]);
    SEGMENT_3D<T> segment1(adhesion.curve.particles.X.Subset(segment1_nodes)),segment2(adhesion.curve.particles.X.Subset(segment2_nodes));

    HAIR_ID segment1_hair_index=adhesion.particle_to_spring_id(segment1_nodes[0]),segment2_hair_index=adhesion.particle_to_spring_id(segment2_nodes[0]);
    assert(segment1_hair_index==adhesion.particle_to_spring_id(segment1_nodes[1]) && segment2_hair_index==adhesion.particle_to_spring_id(segment2_nodes[1]));
    if(segment1_hair_index==segment2_hair_index) adhesion.existing_pairs.Set(segment_indices);
}
//####################################################################
namespace PhysBAM{
template void BOX_HIERARCHY<VECTOR<double,3> >::Intersection_List<EDGE_EDGE_INITIAL_CULL_VISITOR<VECTOR<double,3> > >(BOX_HIERARCHY<VECTOR<double,3> > const&,EDGE_EDGE_INITIAL_CULL_VISITOR<VECTOR<double,3> >&,double) const;
template void BOX_HIERARCHY<VECTOR<float,3> >::Intersection_List<EDGE_EDGE_INITIAL_CULL_VISITOR<VECTOR<float,3> > >(BOX_HIERARCHY<VECTOR<float,3> > const&,EDGE_EDGE_INITIAL_CULL_VISITOR<VECTOR<float,3> >&,float) const;
}
