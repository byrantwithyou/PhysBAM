//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERSECTING_PAIRS_VISITOR
//##################################################################### 
#include <Core/Data_Structures/HASHTABLE.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Intersections/SEGMENT_2D_SEGMENT_2D_INTERSECTION.h>
#include <Geometry/Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <Deformables/Collisions_And_Interactions/INTERSECTING_PAIRS_VISITOR.h>
#include <Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERSECTING_PAIRS_VISITOR<TV>::
INTERSECTING_PAIRS_VISITOR(HASHTABLE<VECTOR<int,d+1> >& intersecting_point_face_pairs,HASHTABLE<VECTOR<int,2*d-2> >& intersecting_edge_edge_pairs,
    const STRUCTURE_INTERACTION_GEOMETRY<TV>& segment_structure,STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure,
    const ARRAY<TV>& X_self_collision_free,const T thickness_over_2,int& count)
    :intersecting_point_face_pairs(intersecting_point_face_pairs),intersecting_edge_edge_pairs(intersecting_edge_edge_pairs),face_structure(face_structure),
    segments(segment_structure.segmented_curve->mesh.elements),faces(face_structure.Face_Mesh_Object()->mesh.elements),
    X_self_collision_free(X_self_collision_free),thickness_over_2(thickness_over_2),count(count)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERSECTING_PAIRS_VISITOR<TV>::
~INTERSECTING_PAIRS_VISITOR()
{}
//#####################################################################
// Function Store_Helper
//#####################################################################
template<class TV> void INTERSECTING_PAIRS_VISITOR<TV>::
Store_Helper(const int segment1,const int segment2,const VECTOR<T,2>* dimension)
{
    const VECTOR<VECTOR<int,2>,2> nodes(segments(segment1),faces(segment2));
    if(INTERSECTION::Intersects(SEGMENT_2D<T>(X_self_collision_free.Subset(nodes[0])),SEGMENT_2D<T>(X_self_collision_free.Subset(nodes[1])),thickness_over_2)) return;
    for(int i=0;i<2;i++){int j=1-i;
        for(int a=0;a<2;a++){
            intersecting_point_face_pairs.Set(VECTOR<int,3>(nodes[i][a],nodes[j][0],nodes[j][1]));
            for(int b=0;b<2;b++) intersecting_edge_edge_pairs.Set(VECTOR<int,2>(nodes[i][a],nodes[j][b]));}}
    count++;
}
//#####################################################################
// Function Store_Helper
//#####################################################################
template<class TV> void INTERSECTING_PAIRS_VISITOR<TV>::
Store_Helper(const int segment,const int triangle,const VECTOR<T,3>* dimension)
{
    const VECTOR<int,2>& segment_nodes=segments(segment);const VECTOR<int,3>& triangle_nodes=faces(triangle);
    if(INTERSECTION::Intersects(SEGMENT_3D<T>(X_self_collision_free.Subset(segment_nodes)),TRIANGLE_3D<T>(X_self_collision_free.Subset(triangle_nodes)),thickness_over_2))
        return;
    if(!face_structure.triangulated_surface->mesh.element_edges) face_structure.triangulated_surface->mesh.Initialize_Element_Edges();
    const VECTOR<int,3>& triangle_edges=(*face_structure.triangulated_surface->mesh.element_edges)(triangle);
    for(int a=0;a<2;a++) intersecting_point_face_pairs.Set(VECTOR<int,4>(segment_nodes[a],triangle_nodes[0],triangle_nodes[1],triangle_nodes[2]));
    for(int a=0;a<3;a++){
        const VECTOR<int,2> other_segment_nodes=segments(triangle_edges[a]); // need to use triangle_edges to get the correct segment orientation
        intersecting_edge_edge_pairs.Set(VECTOR<int,4>(segment_nodes[0],segment_nodes[1],other_segment_nodes[0],other_segment_nodes[1]));
        intersecting_edge_edge_pairs.Set(VECTOR<int,4>(other_segment_nodes[0],other_segment_nodes[1],segment_nodes[0],segment_nodes[1]));}
    count++;
}
//####################################################################
namespace PhysBAM{
template INTERSECTING_PAIRS_VISITOR<VECTOR<float,2> >::INTERSECTING_PAIRS_VISITOR(HASHTABLE<VECTOR<int,3>,void>&,HASHTABLE<VECTOR<int,2>,void>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,2> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,2> >&,ARRAY<VECTOR<float,2>,int> const&,float,int&);
template void INTERSECTING_PAIRS_VISITOR<VECTOR<float,2> >::Store_Helper(int,int,const VECTOR<float,2>* dimension);
template INTERSECTING_PAIRS_VISITOR<VECTOR<float,2> >::~INTERSECTING_PAIRS_VISITOR();
template INTERSECTING_PAIRS_VISITOR<VECTOR<float,3> >::INTERSECTING_PAIRS_VISITOR(HASHTABLE<VECTOR<int,4>,void>&,HASHTABLE<VECTOR<int,4>,void>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,3> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,3> >&,ARRAY<VECTOR<float,3>,int> const&,float,int&);
template void INTERSECTING_PAIRS_VISITOR<VECTOR<float,3> >::Store_Helper(int,int,const VECTOR<float,3>* dimension);
template INTERSECTING_PAIRS_VISITOR<VECTOR<float,3> >::~INTERSECTING_PAIRS_VISITOR();
template INTERSECTING_PAIRS_VISITOR<VECTOR<double,2> >::INTERSECTING_PAIRS_VISITOR(HASHTABLE<VECTOR<int,3>,void>&,HASHTABLE<VECTOR<int,2>,void>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,2> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,2> >&,ARRAY<VECTOR<double,2>,int> const&,double,int&);
template void INTERSECTING_PAIRS_VISITOR<VECTOR<double,2> >::Store_Helper(int,int,const VECTOR<double,2>* dimension);
template INTERSECTING_PAIRS_VISITOR<VECTOR<double,2> >::~INTERSECTING_PAIRS_VISITOR();
template INTERSECTING_PAIRS_VISITOR<VECTOR<double,3> >::INTERSECTING_PAIRS_VISITOR(HASHTABLE<VECTOR<int,4>,void>&,HASHTABLE<VECTOR<int,4>,void>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,3> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,3> >&,ARRAY<VECTOR<double,3>,int> const&,double,int&);
template void INTERSECTING_PAIRS_VISITOR<VECTOR<double,3> >::Store_Helper(int,int,const VECTOR<double,3>* dimension);
template INTERSECTING_PAIRS_VISITOR<VECTOR<double,3> >::~INTERSECTING_PAIRS_VISITOR();
template void BOX_HIERARCHY<VECTOR<double,2> >::Intersection_List<BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<VECTOR<double,2> > > >(
    BOX_HIERARCHY<VECTOR<double,2> > const&,BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<VECTOR<double,2> > >&,double) const;
template void BOX_HIERARCHY<VECTOR<double,2> >::Intersection_List<INTERSECTING_PAIRS_VISITOR<VECTOR<double,2> > >(
    BOX_HIERARCHY<VECTOR<double,2> > const&,INTERSECTING_PAIRS_VISITOR<VECTOR<double,2> >&,double) const;
template void BOX_HIERARCHY<VECTOR<double,3> >::Intersection_List<BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<VECTOR<double,3> > > >(
    BOX_HIERARCHY<VECTOR<double,3> > const&,BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<VECTOR<double,3> > >&,double) const;
template void BOX_HIERARCHY<VECTOR<double,3> >::Intersection_List<INTERSECTING_PAIRS_VISITOR<VECTOR<double,3> > >(
    BOX_HIERARCHY<VECTOR<double,3> > const&,INTERSECTING_PAIRS_VISITOR<VECTOR<double,3> >&,double) const;
template void BOX_HIERARCHY<VECTOR<float,2> >::Intersection_List<BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<VECTOR<float,2> > > >(
    BOX_HIERARCHY<VECTOR<float,2> > const&,BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<VECTOR<float,2> > >&,float) const;
template void BOX_HIERARCHY<VECTOR<float,2> >::Intersection_List<INTERSECTING_PAIRS_VISITOR<VECTOR<float,2> > >(
    BOX_HIERARCHY<VECTOR<float,2> > const&,INTERSECTING_PAIRS_VISITOR<VECTOR<float,2> >&,float) const;
template void BOX_HIERARCHY<VECTOR<float,3> >::Intersection_List<BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<VECTOR<float,3> > > >(
    BOX_HIERARCHY<VECTOR<float,3> > const&,BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<VECTOR<float,3> > >&,float) const;
template void BOX_HIERARCHY<VECTOR<float,3> >::Intersection_List<INTERSECTING_PAIRS_VISITOR<VECTOR<float,3> > >(
    BOX_HIERARCHY<VECTOR<float,3> > const&,INTERSECTING_PAIRS_VISITOR<VECTOR<float,3> >&,float) const;
}
