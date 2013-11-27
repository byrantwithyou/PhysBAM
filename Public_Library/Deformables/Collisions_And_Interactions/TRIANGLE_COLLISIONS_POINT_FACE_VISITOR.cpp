//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISIONS_POINT_FACE_VISITOR
//##################################################################### 
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS_POINT_FACE_VISITOR.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV>::
TRIANGLE_COLLISIONS_POINT_FACE_VISITOR(ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,
    const STRUCTURE_INTERACTION_GEOMETRY<TV>& particle_structure,const STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure,
    const TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry,const T collision_thickness,MPI_SOLIDS<TV>* mpi_solids)
    :pairs_internal(pairs_internal),pairs_external(pairs_external),particle_active_indices(particle_structure.active_indices),
    faces(face_structure.Face_Mesh_Object()->mesh.elements),X(geometry.deformable_body_collection.particles.X),
    X_self_collision_free(geometry.X_self_collision_free),
    collision_thickness(collision_thickness),point_box_modified(particle_structure.point_modified),face_box_modified(face_structure.Face_Modified()),
    intersecting_point_face_pairs(geometry.intersecting_point_face_pairs),mpi_solids(mpi_solids)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV>::
~TRIANGLE_COLLISIONS_POINT_FACE_VISITOR()
{}
//#####################################################################
// Function Store
//#####################################################################
template<class TV> void  TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV>::
Store(const int point_index,const int face_index)
{
    const VECTOR<int,d>& face_nodes=faces(face_index);int p=particle_active_indices(point_index);
    if(face_nodes.Contains(p)) return;
    TV dX_average=(T).5*(X(p)-X_self_collision_free(p)+X.Subset(face_nodes).Average()-X_self_collision_free.Subset(face_nodes).Average());
    RANGE<TV> box1=RANGE<TV>::Bounding_Box(X_self_collision_free(p)+dX_average,X(p));
    RANGE<TV> box2=RANGE<TV>::Bounding_Box(X_self_collision_free.Subset(face_nodes))+dX_average;box2.Enlarge_Nonempty_Box_To_Include_Points(X.Subset(face_nodes));
    if(!box1.Intersection(box2,collision_thickness)) return;
    VECTOR<int,d+1> nodes=face_nodes.Insert(p,0);
    if(intersecting_point_face_pairs.Size() && intersecting_point_face_pairs.Contains(nodes)) return;
    if(mpi_solids){
        VECTOR<PARTITION_ID,d+1> processors(mpi_solids->partition_id_from_particle_index.Subset(nodes));
        for(int i=0;i<d;i++)
            if(processors(i)!=processors(d)){
                pairs_external.Append(nodes);
                return;}
        pairs_internal.Append(nodes);}
    else pairs_internal.Append(nodes);
}
//####################################################################
namespace PhysBAM{
template struct TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,1> >;
template struct TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,2> >;
template struct TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,3> >;
template struct TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,1> >;
template struct TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,2> >;
template struct TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,3> >;
template void BOX_HIERARCHY<VECTOR<double,2> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,2> > > >(
    BOX_HIERARCHY<VECTOR<double,2> > const&,BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,2> > >&,double) const;
template void BOX_HIERARCHY<VECTOR<double,2> >::Intersection_List<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,2> > >(
    BOX_HIERARCHY<VECTOR<double,2> > const&,TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,2> >&,double) const;
template void BOX_HIERARCHY<VECTOR<double,3> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,3> > > >(
    BOX_HIERARCHY<VECTOR<double,3> > const&,BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,3> > >&,double) const;
template void BOX_HIERARCHY<VECTOR<double,3> >::Intersection_List<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,3> > >(
    BOX_HIERARCHY<VECTOR<double,3> > const&,TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,3> >&,double) const;
template void BOX_HIERARCHY<VECTOR<float,2> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,2> > > >(
    BOX_HIERARCHY<VECTOR<float,2> > const&,BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,2> > >&,float) const;
template void BOX_HIERARCHY<VECTOR<float,2> >::Intersection_List<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,2> > >(
    BOX_HIERARCHY<VECTOR<float,2> > const&,TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,2> >&,float) const;
template void BOX_HIERARCHY<VECTOR<float,3> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,3> > > >(
    BOX_HIERARCHY<VECTOR<float,3> > const&,BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,3> > >&,float) const;
template void BOX_HIERARCHY<VECTOR<float,3> >::Intersection_List<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,3> > >(
    BOX_HIERARCHY<VECTOR<float,3> > const&,TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,3> >&,float) const;
}
