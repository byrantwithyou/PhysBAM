//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TRIANGULATED_OBJECT
//#####################################################################
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <Deformables/Fracture/EMBEDDED_TRIANGULATED_OBJECT.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EMBEDDED_TRIANGULATED_OBJECT<TV>::
EMBEDDED_TRIANGULATED_OBJECT(T_TRIANGULATED_OBJECT& simplicial_object_input)
    :EMBEDDED_OBJECT<TV,2>(simplicial_object_input)
{
}
//#####################################################################
// Function Node_Separated_By_Embedded_Subelement
//#####################################################################
template<class TV> int EMBEDDED_TRIANGULATED_OBJECT<TV>::
Node_Separated_By_Embedded_Subelement(const int embedded_segment) const
{
    int global_a,global_b;
    embedded_mesh.elements(embedded_segment).Get(global_a,global_b);
    int a=embedded_particles.subset_index_from_point_cloud_index(global_a),b=embedded_particles.subset_index_from_point_cloud_index(global_b);
    int ppa1,ppa2,ppb1,ppb2;
    parent_particles(a).Get(ppa1,ppa2);
    parent_particles(b).Get(ppb1,ppb2);
    if(Is_Parent(ppa1,b)) return ppa1;
    else return ppa2;
}
//#####################################################################
// Function Nodes_Are_Separated_In_Simplex
//#####################################################################
template<class TV> bool EMBEDDED_TRIANGULATED_OBJECT<TV>::
Nodes_Are_Separated_In_Simplex(const int node1,const int node2,const int triangle) const
{
    int embedded_node=Embedded_Particle_On_Segment(node1,node2);
    if(embedded_node<0) return false;
    VECTOR<int,2> segments=Embedded_Subelements_In_Element(triangle);
    int global_embedded_node=embedded_particles.active_indices(embedded_node);
    if(segments[0]<0) return false;
    else if(embedded_mesh.Node_In_Segment(global_embedded_node,segments[0])) return true;
    if(segments[1]<0) return false;
    else if(embedded_mesh.Node_In_Segment(global_embedded_node,segments[1])) return true;
    return false;
}
//#####################################################################
// Function Segment_Is_Broken
//#####################################################################
template<class TV> bool EMBEDDED_TRIANGULATED_OBJECT<TV>::
Segment_Is_Broken(const int node1,const int node2) const
{
    if(!Embedded_Particle_On_Segment(node1,node2)) return false;
    ARRAY<int> triangles_on_edge;simplicial_object.mesh.Triangles_On_Edge(node1,node2,&triangles_on_edge);
    for(int t=0;t<triangles_on_edge.m;t++) if(!Nodes_Are_Separated_In_Simplex(node1,node2,triangles_on_edge(t)))return false;
    return true;
}
//#####################################################################
// Function Number_Of_Edges_With_Embedded_Particles
//#####################################################################
template<class TV> int EMBEDDED_TRIANGULATED_OBJECT<TV>::
Number_Of_Edges_With_Embedded_Particles(const int triangle)
{
    int i,j,k;simplicial_object.mesh.elements(triangle).Get(i,j,k);
    int ij=Embedded_Particle_On_Segment(i,j),ik=Embedded_Particle_On_Segment(i,k),jk=Embedded_Particle_On_Segment(j,k);
    return (ij>0)+(ik>0)+(jk>0);
}
//#####################################################################
// Function Embedded_Node_Common_To_Both_Segments_In_Triangle
//#####################################################################
template<class TV> int EMBEDDED_TRIANGULATED_OBJECT<TV>::
Embedded_Node_Common_To_Both_Segments_In_Triangle(const int triangle)
{
    VECTOR<int,2> emb_segments=Embedded_Subelements_In_Element(triangle);
    assert(emb_segments[0] && emb_segments[1]);
    int global_a,global_b,global_c,global_d;
    embedded_mesh.elements(emb_segments[0]).Get(global_a,global_b);
    embedded_mesh.elements(emb_segments[1]).Get(global_c,global_d);
    return embedded_particles.subset_index_from_point_cloud_index((global_a==global_c || global_a==global_d)?global_a:global_b);
}
//#####################################################################
// Function Isolated_Node
//#####################################################################
template<class TV> int EMBEDDED_TRIANGULATED_OBJECT<TV>::
Isolated_Node(const int triangle)
{
    VECTOR<int,2> emb_segments=Embedded_Subelements_In_Element(triangle);
    assert(emb_segments[0] && !emb_segments[1]);
    return Node_Separated_By_Embedded_Subelement(emb_segments[0]);
}
//#####################################################################
// Function Diamond_Node
//#####################################################################
template<class TV> int EMBEDDED_TRIANGULATED_OBJECT<TV>::
Diamond_Node(const int triangle)
{
    if(Number_Of_Embedded_Subelements_In_Element(triangle)<2) return -1;
    int common_emb_node=Embedded_Node_Common_To_Both_Segments_In_Triangle(triangle);
    int i,j,k;
    simplicial_object.mesh.elements(triangle).Get(i,j,k);
    if(!Is_Parent(i,common_emb_node)) return i;
    if(!Is_Parent(j,common_emb_node)) return j;
    return k;
}
//#####################################################################
// Function Fraction_Of_Triangles_With_N_Cuts
//#####################################################################
template<class TV> auto EMBEDDED_TRIANGULATED_OBJECT<TV>::
Fraction_Of_Triangles_With_N_Cuts(const int n) -> T
{
    int count=0;
    for(int t=0;t<simplicial_object.mesh.elements.m;t++) if(Number_Of_Embedded_Cuts(t) == n) count++;
    return (T)count/(T)simplicial_object.mesh.elements.m;
}
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,2> >;
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,3> >;
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,2> >;
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,3> >;
}
