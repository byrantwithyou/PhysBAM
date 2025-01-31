//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TETRAHEDRALIZED_VOLUME
//#####################################################################
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
EMBEDDED_TETRAHEDRALIZED_VOLUME(TETRAHEDRALIZED_VOLUME<T>& simplicial_object_input)
    :EMBEDDED_OBJECT<TV,3>(simplicial_object_input)
{
}
//#####################################################################
// Function Node_Separated_By_Embedded_Subelement
//#####################################################################
template<class T> int EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Node_Separated_By_Embedded_Subelement(const int embedded_triangle) const
{
    int global_a,global_b,global_c;embedded_mesh.elements(embedded_triangle).Get(global_a,global_b,global_c);
    int a=embedded_particles.subset_index_from_point_cloud_index(global_a),b=embedded_particles.subset_index_from_point_cloud_index(global_b),
        c=embedded_particles.subset_index_from_point_cloud_index(global_c);
    int ppa1,ppa2,ppb1,ppb2,ppc1,ppc2;parent_particles(a).Get(ppa1,ppa2);parent_particles(b).Get(ppb1,ppb2);parent_particles(c).Get(ppc1,ppc2);
    if(Is_Parent(ppa1,b) && Is_Parent(ppa1,c))return ppa1;
    else if(Is_Parent(ppa2,b) && Is_Parent(ppa2,c))return ppa2;
    else return -1;
}
//#####################################################################
// Function Nodes_Are_Separated_In_Simplex
//#####################################################################
template<class T> bool EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Nodes_Are_Separated_In_Simplex(const int node1,const int node2,const int tetrahedron) const
{
    int embedded_node=Embedded_Particle_On_Segment(node1,node2);if(embedded_node<0)return false;
    VECTOR<int,4> triangles=Embedded_Subelements_In_Element(tetrahedron);
    int global_embedded_node=embedded_particles.active_indices(embedded_node);
    for(int i=0;i<4 && triangles[i];i++){
        if(embedded_mesh.Node_In_Triangle(global_embedded_node,triangles[i])) return true;}
    return false;
}
//#####################################################################
// Function Nodes_Are_Materially_Connected_In_Simplex
//#####################################################################
template<class T> bool EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Nodes_Are_Materially_Connected_In_Simplex(const int node1,const int node2,const int simplex) const
{
    return Node_In_Simplex_Is_Material(node1,simplex) && Node_In_Simplex_Is_Material(node2,simplex) && !Nodes_Are_Separated_In_Simplex(node1,node2,simplex);
}
//#####################################################################
// Function Number_Of_Embedded_Cuts
//#####################################################################
template<class T> int EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Number_Of_Embedded_Cuts(const int tetrahedron) const // a quad only counts as one cut
{
    int triangle_count=Number_Of_Embedded_Subelements_In_Element(tetrahedron);
    if(Cut_By_Quad(tetrahedron)) return triangle_count-1;
    else return triangle_count;
}
//#####################################################################
// Function Cut_By_Quad
//#####################################################################
template<class T> bool EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Cut_By_Quad(const int tetrahedron) const
{
    int triangle_count=Number_Of_Embedded_Subelements_In_Element(tetrahedron);
    if(triangle_count <= 1) return false;
    VECTOR<int,4> triangles=Embedded_Subelements_In_Element(tetrahedron);
    if(triangle_count == 2) return !Node_Separated_By_Embedded_Subelement(triangles[0]);
    if(triangle_count == 3) return !Node_Separated_By_Embedded_Subelement(triangles[0]) || !Node_Separated_By_Embedded_Subelement(triangles[1])
                                || !Node_Separated_By_Embedded_Subelement(triangles[2]);
    return true;
} // if four embedded triangles, there must be a quad
//#####################################################################
// Function Segment_Is_Broken
//#####################################################################
template<class T> bool EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Segment_Is_Broken(const int node1,const int node2) const
{
    if(!Embedded_Particle_On_Segment(node1,node2)) return false;
    ARRAY<int> tetrahedrons_on_edge;simplicial_object.mesh.Tetrahedrons_On_Edge(VECTOR<int,2>(node1,node2),tetrahedrons_on_edge);
    for(int t=0;t<tetrahedrons_on_edge.m;t++) if(!Nodes_Are_Separated_In_Simplex(node1,node2,tetrahedrons_on_edge(t))) return false;
    return true;
}
//#####################################################################
// Function Number_Of_Edges_With_Embedded_Particles
//#####################################################################
template<class T> int EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Number_Of_Edges_With_Embedded_Particles(const int tetrahedron) const
{
    int i,j,k,l;simplicial_object.mesh.elements(tetrahedron).Get(i,j,k,l);
    int ij=Embedded_Particle_On_Segment(i,j),ik=Embedded_Particle_On_Segment(i,k),il=Embedded_Particle_On_Segment(i,l),
        jk=Embedded_Particle_On_Segment(j,k),jl=Embedded_Particle_On_Segment(j,l),kl=Embedded_Particle_On_Segment(k,l);
    return (ij>0)+(ik>0)+(il>0)+(jk>0)+(jl>0)+(kl>0);
}
//#####################################################################
// Function Number_Of_Tetrahedra_With_Cuts
//#####################################################################
template<class T> int EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Number_Of_Tetrahedra_With_Cuts() const
{
    int count=0;for(int t=0;t<simplicial_object.mesh.elements.m;t++) if(Number_Of_Embedded_Cuts(t) > 0) count++;
    return count;
}
//#####################################################################
// Function Number_Of_Tetrahedra_With_N_Cuts
//#####################################################################
template<class T> int EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
Number_Of_Tetrahedra_With_N_Cuts(const int n) const
{
    int count=0;for(int t=0;t<simplicial_object.mesh.elements.m;t++) if(Number_Of_Embedded_Cuts(t) == n) count++;
    return count;
}
//#####################################################################
template class EMBEDDED_TETRAHEDRALIZED_VOLUME<float>;
template class EMBEDDED_TETRAHEDRALIZED_VOLUME<double>;
}
