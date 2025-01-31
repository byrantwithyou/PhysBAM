//#####################################################################
// Copyright 2002-2006, Zhaosheng Bao, Robert Bridson, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Neil Molino, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRON_MESH
//##################################################################### 
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/STACK.h>
#include <Core/Data_Structures/UNION_FIND.h>
#include <Core/Math_Tools/cyclic_shift.h>
#include <Core/Math_Tools/exchange_sort.h>
#include <Core/Math_Tools/permutation.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology/TETRAHEDRON_MESH.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
TETRAHEDRON_MESH::
TETRAHEDRON_MESH() // simplest constructor - null mesh
    :segment_mesh(0),triangle_mesh(0),element_edges(0),tetrahedron_faces(0),boundary_mesh(0),node_on_boundary(0),boundary_nodes(0),edge_tetrahedrons(0),triangle_tetrahedrons(0),
    owns_segment_mesh(true)
{}
//#####################################################################
// Constructor
//#####################################################################
TETRAHEDRON_MESH::
TETRAHEDRON_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,4> >& tetrahedron_list)
    :SIMPLEX_MESH<3>(number_nodes_input,tetrahedron_list),segment_mesh(0),triangle_mesh(0),element_edges(0),tetrahedron_faces(0),boundary_mesh(0),node_on_boundary(0),
    boundary_nodes(0),edge_tetrahedrons(0),triangle_tetrahedrons(0),owns_segment_mesh(true)
{}
//#####################################################################
// Constructor
//#####################################################################
TETRAHEDRON_MESH::
TETRAHEDRON_MESH(const TETRAHEDRON_MESH& tetrahedron_mesh)
    :SIMPLEX_MESH<3>(tetrahedron_mesh),segment_mesh(0),triangle_mesh(0),element_edges(0),tetrahedron_faces(0),boundary_mesh(0),node_on_boundary(0),boundary_nodes(0),
    edge_tetrahedrons(0),triangle_tetrahedrons(0),owns_segment_mesh(true)
{}
//#####################################################################
// Destructor
//#####################################################################
TETRAHEDRON_MESH::
~TETRAHEDRON_MESH()
{
    if(owns_segment_mesh) delete segment_mesh;
    delete triangle_mesh;
    delete element_edges;
    delete tetrahedron_faces;
    delete boundary_mesh;
    delete node_on_boundary;
    delete boundary_nodes;
    delete edge_tetrahedrons;
    delete triangle_tetrahedrons;
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
void TETRAHEDRON_MESH::
Delete_Auxiliary_Structures()
{
    SIMPLEX_MESH<3>::Delete_Auxiliary_Structures();
    if(owns_segment_mesh) delete segment_mesh;
    segment_mesh=0;
    delete triangle_mesh;
    triangle_mesh=0;
    delete element_edges;
    element_edges=0;
    delete tetrahedron_faces;
    tetrahedron_faces=0;
    delete boundary_mesh;
    boundary_mesh=0;
    delete node_on_boundary;
    node_on_boundary=0;
    delete boundary_nodes;
    boundary_nodes=0;
    delete edge_tetrahedrons;
    edge_tetrahedrons=0;
    delete triangle_tetrahedrons;
    triangle_tetrahedrons=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures
//#####################################################################
void TETRAHEDRON_MESH::
Refresh_Auxiliary_Structures()
{
    SIMPLEX_MESH<3>::Refresh_Auxiliary_Structures();
    if(segment_mesh) Initialize_Segment_Mesh();
    if(triangle_mesh) Initialize_Triangle_Mesh();
    if(element_edges) Initialize_Element_Edges();
    if(tetrahedron_faces) Initialize_Tetrahedron_Faces();
    if(boundary_mesh) Initialize_Boundary_Mesh();
    if(node_on_boundary) Initialize_Node_On_Boundary();
    if(boundary_nodes) Initialize_Boundary_Nodes();
    if(edge_tetrahedrons) Initialize_Edge_Tetrahedrons();
    if(triangle_tetrahedrons) Initialize_Triangle_Tetrahedrons();
}
//#####################################################################
// Function Initialize_Octahedron_Mesh
//#####################################################################
static inline int Lattice(const int i,const int j,const int k,const int m,const int n,const int p)
{
    return i+m*j+m*n*k;
}
static inline int Half_Lattice(const int i,const int j,const int k,const int m,const int n,const int p)
{
    return m*n*p+i+1+(m+1)*(j+1)+(m+1)*(n+1)*(k+1);
}
void TETRAHEDRON_MESH::
Initialize_Octahedron_Mesh(const int m,const int n,const int p)
{
    Clean_Memory();
    number_nodes=m*n*p+(m+1)*(n+1)*(p+1);elements.Exact_Resize(4*(m-1)*n*p+4*m*(n-1)*p+4*m*n*(p-1));
    int t=0,i,j,k;
    for(i=0;i<m;i++)for(j=0;j<n;j++)for(k=0;k<p-1;k++){ // loop over k-oriented edges in inner cube mesh
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i,j,k+1,m,n,p),Half_Lattice(i-1,j-1,k,m,n,p),Half_Lattice(i,j-1,k,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i,j,k+1,m,n,p),Half_Lattice(i,j-1,k,m,n,p),Half_Lattice(i,j,k,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i,j,k+1,m,n,p),Half_Lattice(i,j,k,m,n,p),Half_Lattice(i-1,j,k,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i,j,k+1,m,n,p),Half_Lattice(i-1,j,k,m,n,p),Half_Lattice(i-1,j-1,k,m,n,p));}
    for(i=0;i<m;i++)for(j=0;j<n-1;j++)for(k=0;k<p;k++){ // loop over j-oriented edge in inner cube mesh
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i,j+1,k,m,n,p),Half_Lattice(i,j,k-1,m,n,p),Half_Lattice(i-1,j,k-1,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i,j+1,k,m,n,p),Half_Lattice(i,j,k,m,n,p),Half_Lattice(i,j,k-1,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i,j+1,k,m,n,p),Half_Lattice(i-1,j,k,m,n,p),Half_Lattice(i,j,k,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i,j+1,k,m,n,p),Half_Lattice(i-1,j,k-1,m,n,p),Half_Lattice(i-1,j,k,m,n,p));}
    for(i=0;i<m-1;i++)for(j=0;j<n;j++)for(k=0;k<p;k++){ // loop over i-oriented edge in inner cube mesh
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j,k,m,n,p),Half_Lattice(i,j-1,k-1,m,n,p),Half_Lattice(i,j,k-1,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j,k,m,n,p),Half_Lattice(i,j,k-1,m,n,p),Half_Lattice(i,j,k,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j,k,m,n,p),Half_Lattice(i,j,k,m,n,p),Half_Lattice(i,j-1,k,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j,k,m,n,p),Half_Lattice(i,j-1,k,m,n,p),Half_Lattice(i,j-1,k-1,m,n,p));}
}
//#####################################################################
// Function Initialize_Cube_Mesh
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Cube_Mesh(const int m,const int n,const int p) // 5 tetrahedrons per cube
{
    Clean_Memory();number_nodes=m*n*p;elements.Exact_Resize(5*(m-1)*(n-1)*(p-1));int t=0;
    for(int i=0;i<m-1;i++)for(int j=0;j<n-1;j++)for(int k=0;k<p-1;k++){
        if((i+j+k)%2 != 0){
            elements(t++).Set(i+m*j+m*n*k,i+1+m*j+m*n*k,i+m*(j+1)+m*n*k,i+m*j+m*n*(k+1));
            elements(t++).Set(i+1+m*j+m*n*k,i+1+m*j+m*n*(k+1),i+1+m*(j+1)+m*n*(k+1),i+m*j+m*n*(k+1));
            elements(t++).Set(i+m*(j+1)+m*n*k,i+1+m*(j+1)+m*n*k,i+1+m*(j+1)+m*n*(k+1),i+1+m*j+m*n*k);
            elements(t++).Set(i+m*(j+1)+m*n*(k+1),i+1+m*(j+1)+m*n*(k+1),i+m*j+m*n*(k+1),i+m*(j+1)+m*n*k);
            elements(t++).Set(i+1+m*j+m*n*k,i+m*j+m*n*(k+1),i+1+m*(j+1)+m*n*(k+1),i+m*(j+1)+m*n*k);}
        else{
            elements(t++).Set(i+m*j+m*n*k,i+1+m*j+m*n*k,i+1+m*(j+1)+m*n*k,i+1+m*j+m*n*(k+1));
            elements(t++).Set(i+m*j+m*n*k,i+m*(j+1)+m*n*k,i+m*(j+1)+m*n*(k+1),i+1+m*(j+1)+m*n*k);
            elements(t++).Set(i+m*(j+1)+m*n*(k+1),i+1+m*j+m*n*(k+1),i+m*j+m*n*(k+1),i+m*j+m*n*k);
            elements(t++).Set(i+m*(j+1)+m*n*(k+1),i+1+m*(j+1)+m*n*(k+1),i+1+m*j+m*n*(k+1),i+1+m*(j+1)+m*n*k);
            elements(t++).Set(i+m*(j+1)+m*n*(k+1),i+m*j+m*n*k,i+1+m*(j+1)+m*n*k,i+1+m*j+m*n*(k+1));}}
}
//#####################################################################
// Function Initialize_Prismatic_Cube_Mesh
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Prismatic_Cube_Mesh(const int m,const int n,const int p) // 6 tetrahedra per cube
{
    Clean_Memory();number_nodes=m*n*p;elements.Exact_Resize(6*(m-1)*(n-1)*(p-1));int t=0;
    for(int i=0;i<m-1;i++)for(int j=0;j<n-1;j++)for(int k=0;k<p-1;k++){
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j,k,m,n,p),Lattice(i+1,j+1,k,m,n,p),Lattice(i+1,j+1,k+1,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j,k,m,n,p),Lattice(i+1,j+1,k+1,m,n,p),Lattice(i+1,j,k+1,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j,k+1,m,n,p),Lattice(i+1,j+1,k+1,m,n,p),Lattice(i,j,k+1,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j+1,k+1,m,n,p),Lattice(i,j+1,k+1,m,n,p),Lattice(i,j,k+1,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j+1,k+1,m,n,p),Lattice(i,j+1,k,m,n,p),Lattice(i,j+1,k+1,m,n,p));
        elements(t++).Set(Lattice(i,j,k,m,n,p),Lattice(i+1,j+1,k,m,n,p),Lattice(i,j+1,k,m,n,p),Lattice(i+1,j+1,k+1,m,n,p));}
}
//#####################################################################
// Function Initialize_Segment_Mesh
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Segment_Mesh()
{
    if(owns_segment_mesh) delete segment_mesh;else owns_segment_mesh=true;
    bool neighbor_nodes_defined=neighbor_nodes!=0;if(!neighbor_nodes) Initialize_Neighbor_Nodes();
    // number of edges = half the sum of the degree (or a little more if there are loop edges)
    int total_degree=0;for(int i=0;i<number_nodes;i++) total_degree+=(*neighbor_nodes)(i).m;
    segment_mesh=new SEGMENT_MESH();
    segment_mesh->number_nodes=number_nodes;
    segment_mesh->elements.Preallocate((total_degree+1)/2); // add one for subtle optimization purposes
    for(int node=0;node<neighbor_nodes->m;node++) for(int k=0;k<(*neighbor_nodes)(node).m;k++) // do nodes in ascending order so that no edges are counted more than once
        if(node<=(*neighbor_nodes)(node)(k)) segment_mesh->elements.Append(VECTOR<int,2>(node,(*neighbor_nodes)(node)(k))); 
    if(!neighbor_nodes_defined){delete neighbor_nodes;neighbor_nodes=0;}
}
//#####################################################################
// Function Initialize_Segment_Mesh_From_Triangle_Mesh
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Segment_Mesh_From_Triangle_Mesh()
{
    if(!triangle_mesh->segment_mesh) PHYSBAM_FATAL_ERROR();
    if(owns_segment_mesh) delete segment_mesh;
    segment_mesh=triangle_mesh->segment_mesh;owns_segment_mesh=false;
}
//#####################################################################
// Function Initialize_Triangle_Mesh
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Triangle_Mesh()
{
    delete triangle_mesh;
    HASHTABLE<VECTOR<int,3> > triangle_list(2*elements.m); // list of triangles currently found
    for(int t=0;t<elements.m;t++){
        VECTOR<int,4> sorted_nodes=elements(t).Sorted();
        for(int i=0;i<sorted_nodes.m;i++){VECTOR<int,3> triangle=sorted_nodes.Remove_Index(i);
            triangle_list.Set(triangle);}}
    triangle_mesh=new TRIANGLE_MESH();
    triangle_list.Get_Keys(triangle_mesh->elements);
    triangle_mesh->number_nodes=number_nodes;
}
//#####################################################################
// Function Initialize_Element_Edges
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Element_Edges()
{
    delete element_edges;element_edges=new ARRAY<VECTOR<int,6> >(elements.m);
    if(!segment_mesh) Initialize_Segment_Mesh(); // edges only makes sense when referring to a segment mesh
    // incident_segments make the SEGMENT_MESH::Segment() function faster
    bool incident_segments_defined=segment_mesh->incident_elements!=0;
    if(!segment_mesh->incident_elements) segment_mesh->Initialize_Incident_Elements();
    // for each tetrahedron, make the lists of edges
    for(int t=0;t<elements.m;t++){
        int i,j,k,l;elements(t).Get(i,j,k,l);
        (*element_edges)(t).Set(segment_mesh->Segment(i,j),segment_mesh->Segment(j,k),segment_mesh->Segment(k,i),
            segment_mesh->Segment(i,l),segment_mesh->Segment(j,l),segment_mesh->Segment(k,l));}
    if(!incident_segments_defined){delete segment_mesh->incident_elements;segment_mesh->incident_elements=0;}
}
//#####################################################################
// Function Initialize_Tetrahedron_Faces
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Tetrahedron_Faces()
{
    delete tetrahedron_faces;tetrahedron_faces=new ARRAY<VECTOR<int,4> >(elements.m);
    if(!triangle_mesh) Initialize_Triangle_Mesh(); //triangles only make sense when referring to a triangle mesh
    bool triangle_tetrahedrons_defined=triangle_tetrahedrons!=0;if(!triangle_tetrahedrons_defined) Initialize_Triangle_Tetrahedrons();
    for(int tri=0;tri<triangle_tetrahedrons->m;tri++){
        int tri_i,tri_j,tri_k;triangle_mesh->elements(tri).Get(tri_i,tri_j,tri_k);
        for(int face=0;face<2;face++){
            int tet=(*triangle_tetrahedrons)(tri)(face),tet_i,tet_j,tet_k,tet_l,local_tet_index_for_face=1;
            if(tet>=0){
                elements(tet).Get(tet_i,tet_j,tet_k,tet_l);
                if(tet_i==tri_i||tet_i==tri_j||tet_i==tri_k)
                    if(tet_j==tri_i||tet_j==tri_j||tet_j==tri_k)
                        if(tet_k==tri_i||tet_k==tri_j||tet_k==tri_k)
                            if(tet_l==tri_i||tet_l==tri_j||tet_l==tri_k) local_tet_index_for_face=0;
                            else local_tet_index_for_face=4;
                        else local_tet_index_for_face=3;
                    else local_tet_index_for_face=2;
                else local_tet_index_for_face=1;
                (*tetrahedron_faces)(tet)(local_tet_index_for_face)=tri;}}}
    if(!triangle_tetrahedrons_defined){delete triangle_tetrahedrons;triangle_tetrahedrons=0;}
}
//#####################################################################
// Function Number_Of_Tetrahedrons_Across_Face
//#####################################################################
int TETRAHEDRON_MESH::
Number_Of_Tetrahedrons_Across_Face(const int tetrahedron,const int node1,const int node2,const int node3) const
{
    int count=0;
    for(int t=0;t<(*incident_elements)(node1).m;t++){
        int tetrahedron2=(*incident_elements)(node1)(t); // tetrahedron in question
        int i,j,k,l;elements(tetrahedron2).Get(i,j,k,l);
        if(tetrahedron2 != tetrahedron){
            if(i == node2 && (j == node3 || k == node3 || l == node3)) count++;
            else if(j == node2 && (i == node3 || k == node3 || l == node3)) count++;
            else if(k == node2 && (i == node3 || j == node3 || l == node3)) count++;
            else if(l == node2 && (i == node3 || j == node3 || k == node3)) count++;}}
    return count;
}
//#####################################################################
// Function Initialize_Boundary_Mesh
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Boundary_Mesh()
{
    delete boundary_mesh;boundary_mesh=new TRIANGLE_MESH();
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    for(int t=0;t<elements.m;t++){
        int i,j,k,l;elements(t).Get(i,j,k,l);
        if(Number_Of_Tetrahedrons_Across_Face(t,i,k,j) == 0)boundary_mesh->elements.Append(VECTOR<int,3>(i,k,j));
        if(Number_Of_Tetrahedrons_Across_Face(t,i,j,l) == 0)boundary_mesh->elements.Append(VECTOR<int,3>(i,j,l));
        if(Number_Of_Tetrahedrons_Across_Face(t,i,l,k) == 0)boundary_mesh->elements.Append(VECTOR<int,3>(i,l,k));
        if(Number_Of_Tetrahedrons_Across_Face(t,j,k,l) == 0)boundary_mesh->elements.Append(VECTOR<int,3>(j,k,l));}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
    boundary_mesh->number_nodes=number_nodes; // use the same number of nodes as in the tetrahedron mesh
}
//#####################################################################
// Function Initialize_Node_On_Boundary
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Node_On_Boundary()
{
    delete node_on_boundary;node_on_boundary=new ARRAY<bool>(number_nodes);
    bool boundary_mesh_defined=boundary_mesh!=0;if(!boundary_mesh_defined) Initialize_Boundary_Mesh();
    for(int t=0;t<boundary_mesh->elements.m;t++){
        int i,j,k;boundary_mesh->elements(t).Get(i,j,k);
        (*node_on_boundary)(i)=(*node_on_boundary)(j)=(*node_on_boundary)(k)=true;}
    if(!boundary_mesh_defined){delete boundary_mesh;boundary_mesh=0;}
}
//#####################################################################
// Function Initialize_Boundary_Nodes
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Boundary_Nodes()
{
    delete boundary_nodes;boundary_nodes=new ARRAY<int>;
    bool boundary_mesh_defined=boundary_mesh!=0;if(!boundary_mesh_defined) Initialize_Boundary_Mesh();
    for(int t=0;t<boundary_mesh->elements.m;t++){
        int i,j,k;boundary_mesh->elements(t).Get(i,j,k);
        boundary_nodes->Append(i);boundary_nodes->Append(j);boundary_nodes->Append(k);}
    if(!boundary_mesh_defined){delete boundary_mesh;boundary_mesh=0;}

    boundary_nodes->Sort();
    Prune_Duplicates(*boundary_nodes);
}
//#####################################################################
// Function Initialize_Edge_Tetrahedrons
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Edge_Tetrahedrons()
{  
    delete edge_tetrahedrons;
    if(!segment_mesh) Initialize_Segment_Mesh(); // edges only makes sense when referring to a segment mesh
    edge_tetrahedrons=new ARRAY<ARRAY<int> >(segment_mesh->elements.m);
    bool element_edges_defined=element_edges!=0;if(!element_edges_defined) Initialize_Element_Edges();
    for(int t=0;t<elements.m;t++) for(int i=0;i<6;i++) (*edge_tetrahedrons)((*element_edges)(t)(i)).Append(t);
    for(int i=0;i<segment_mesh->elements.m;i++) (*edge_tetrahedrons)(i).Compact();
    if(!element_edges_defined){delete element_edges;element_edges=0;}
}
//#####################################################################
// Function Initialize_Triangle_Tetrahedrons
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Triangle_Tetrahedrons()
{  
    delete triangle_tetrahedrons;
    if(!triangle_mesh) Initialize_Triangle_Mesh();
    triangle_tetrahedrons=new ARRAY<VECTOR<int,2> >(triangle_mesh->elements.m);
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements_defined) Initialize_Incident_Elements();
    for(int t=0;t<triangle_mesh->elements.m;t++){
         const VECTOR<int,3>& triangle_nodes=triangle_mesh->elements(t);int i=triangle_nodes[0];
         for(int tet=0;tet<(*incident_elements)(i).m;tet++){
             const VECTOR<int,4>& original_nodes=elements((*incident_elements)(i)(tet));
             VECTOR<int,4> nodes;
             int p=0;for(;p<24;p++){
                 nodes=permute_four(original_nodes,p);
                 if(nodes.Remove_Index(1)==triangle_nodes) break;}
             if(p<24) (*triangle_tetrahedrons)(t)(permutation_of_four_is_even(p)?1:0)=(*incident_elements)(i)(tet);}}
    // if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Number_Of_Boundary_Tetrahedrons
//#####################################################################
int TETRAHEDRON_MESH::
Number_Of_Boundary_Tetrahedrons()
{
    bool adjacent_elements_defined=adjacent_elements!=0;if(!adjacent_elements) Initialize_Adjacent_Elements();
    int number=0;for(int t=0;t<elements.m;t++) if((*adjacent_elements)(t).m != 4) number++;
    if(!adjacent_elements_defined){delete adjacent_elements;adjacent_elements=0;}
    return number;
}
//#####################################################################
// Function Number_Of_Interior_Tetrahedrons
//#####################################################################
int TETRAHEDRON_MESH::
Number_Of_Interior_Tetrahedrons()
{  
    bool adjacent_elements_defined=adjacent_elements!=0;if(!adjacent_elements) Initialize_Adjacent_Elements();   
    int number=0;for(int t=0;t<elements.m;t++) if((*adjacent_elements)(t).m == 4) number++;  
    if(!adjacent_elements_defined){delete adjacent_elements;adjacent_elements=0;}
    return number;
}
//#####################################################################
// Function Edge_Neighbours
//#####################################################################
bool TETRAHEDRON_MESH::
Edge_Neighbors(const int tet1,const int tet2) const
{
    int i1,j1,k1,l1;elements(tet1).Get(i1,j1,k1,l1);
    int i2,j2,k2,l2;elements(tet2).Get(i2,j2,k2,l2);
    return (i1==i2 || i1==j2 || i1==k2 || i1==l2)+(j1==i2 || j1==j2 || j1== k2 || j1==l2)+(k1==i2 || k1==j2 || k1==k2 || k1==l2)+(l1==i2 || l1==j2 || l1==k2 || l1==l2) >= 2;
}
//#####################################################################
// Function Number_Of_Edge_Neighbors
//#####################################################################
int TETRAHEDRON_MESH::
Number_Of_Edge_Neighbors(const int segment) const
{
    assert(segment_mesh); // edges only makes sense when referring to a segment mesh
    if(edge_tetrahedrons) return (*edge_tetrahedrons)(segment).m;
    else{
        assert(incident_elements);
        int count=0,node1,node2;segment_mesh->elements(segment).Get(node1,node2);
        for(int t=0;t<(*incident_elements)(node1).m;t++) if(Node_In_Tetrahedron(node2,(*incident_elements)(node1)(t))) count++; 
        return count;}
}
//#####################################################################
// Function Initialize_Segment_Mesh_Of_Subset
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Segment_Mesh_Of_Subset(SEGMENT_MESH& segment_mesh_of_subset,const ARRAY<bool>& subset) const
{
    segment_mesh_of_subset.Clean_Memory();int t,i,j,k,l;
    ARRAY<ARRAY<int> > higher_neighbors(number_nodes); // figure out all higher-numbered neighbors first
    for(t=0;t<elements.m;t++) if(subset(t)){
        elements(t).Get(i,j,k,l);
        if(i < j) higher_neighbors(i).Append_Unique(j);else higher_neighbors(j).Append_Unique(i);
        if(i < k) higher_neighbors(i).Append_Unique(k);else higher_neighbors(k).Append_Unique(i);
        if(i < l) higher_neighbors(i).Append_Unique(l);else higher_neighbors(l).Append_Unique(i);
        if(j < k) higher_neighbors(j).Append_Unique(k);else higher_neighbors(k).Append_Unique(j);
        if(j < l) higher_neighbors(j).Append_Unique(l);else higher_neighbors(l).Append_Unique(j);
        if(k < l) higher_neighbors(k).Append_Unique(l);else higher_neighbors(l).Append_Unique(k);}
    int number_segments=0;for(i=0;i<number_nodes;i++) number_segments+=higher_neighbors(i).m;
    segment_mesh_of_subset.elements.Exact_Resize(number_segments);
    t=0;for(i=0;i<number_nodes;i++) for(j=0;j<higher_neighbors(i).m;j++) segment_mesh_of_subset.elements(t++).Set(i,higher_neighbors(i)(j));
    segment_mesh_of_subset.number_nodes=number_nodes;
}
//#####################################################################
// Function Initialize_Boundary_Mesh_Of_Subset
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Boundary_Mesh_Of_Subset(TRIANGLE_MESH& boundary_mesh_of_subset,const ARRAY<bool>& subset)
{
    boundary_mesh_of_subset.Clean_Memory();
    boundary_mesh_of_subset.number_nodes=number_nodes;
    bool adjacent_elements_defined=(adjacent_elements!=0);if(!adjacent_elements) Initialize_Adjacent_Elements();
    // check tetrahedra for out-of-subset neighbors 
    for(int t=0;t<elements.m;t++) if(subset(t)){
        int i,j,k,l;elements(t).Get(i,j,k,l);
        VECTOR<ARRAY<int>,4> adjacent_tets_per_face;
        for(int p=0;p<(*adjacent_elements)(t).m;p++)
            if(!Node_In_Tetrahedron(i,(*adjacent_elements)(t)(p))) adjacent_tets_per_face(0).Append((*adjacent_elements)(t)(p));
            else if(!Node_In_Tetrahedron(j,(*adjacent_elements)(t)(p))) adjacent_tets_per_face(1).Append((*adjacent_elements)(t)(p));
            else if(!Node_In_Tetrahedron(k,(*adjacent_elements)(t)(p))) adjacent_tets_per_face(2).Append((*adjacent_elements)(t)(p));
            else adjacent_tets_per_face(3).Append((*adjacent_elements)(t)(p));
        for(int p=0;p<4;p++){
            bool lowest_index_tet=true,boundary_face=false;
            // For non-manifold meshes, add a boundary face only for the lowest indexed incident tet
            if(!adjacent_tets_per_face(p).m) boundary_face=true;
            for(int q=0;q<adjacent_tets_per_face(p).m;q++)
                if(!subset(adjacent_tets_per_face(p)(q))) boundary_face=true;
                else if(adjacent_tets_per_face(p)(q)<t) lowest_index_tet=false;
            if(boundary_face&&lowest_index_tet) switch(p){
                case 0:boundary_mesh_of_subset.elements.Append(VECTOR<int,3>(j,k,l));break;
                case 1:boundary_mesh_of_subset.elements.Append(VECTOR<int,3>(i,l,k));break;
                case 2:boundary_mesh_of_subset.elements.Append(VECTOR<int,3>(i,j,l));break;
                case 3:boundary_mesh_of_subset.elements.Append(VECTOR<int,3>(i,k,j));break;}}}
    if(!adjacent_elements_defined){delete adjacent_elements;adjacent_elements=0;}
}
//#####################################################################
// Function Add_Triangle_Or_Subtriangles
//#####################################################################
static bool Add_Triangle_Or_Subtriangles(ARRAY<VECTOR<int,3> >& triangle_list,const TRIANGLE_MESH& triangle_mesh,const SEGMENT_MESH& segment_mesh,ARRAY<int>& segment_midpoints,const VECTOR<int,3>& nodes)
{
    int i,j,k;nodes.Get(i,j,k);
    int ij=segment_mesh.Segment(i,j),jk=segment_mesh.Segment(j,k),ik=segment_mesh.Segment(i,k);
    if(ij>=0) ij=segment_midpoints(ij);
    if(jk>=0) jk=segment_midpoints(jk);
    if(ik>=0) ik=segment_midpoints(ik);
    if(ij<0 && jk<0 && ik<0){triangle_list.Append(VECTOR<int,3>(i,j,k));return triangle_mesh.Triangle(i,j,k)!=0;}
    ARRAY<VECTOR<int,3> > subtriangles;
    if(ij>=0 && jk<0 && ik<0){subtriangles.Append(VECTOR<int,3>(i,ij,k));subtriangles.Append(VECTOR<int,3>(ij,j,k));}
    else if(ij<0 && jk>=0 && ik<0){subtriangles.Append(VECTOR<int,3>(i,j,jk));subtriangles.Append(VECTOR<int,3>(i,jk,k));}
    else if(ij<0 && jk<0 && ik>=0){subtriangles.Append(VECTOR<int,3>(i,j,ik));subtriangles.Append(VECTOR<int,3>(ik,j,k));}
    else if(ij>=0 && jk>=0 && ik>=0){subtriangles.Append(VECTOR<int,3>(i,ij,ik));subtriangles.Append(VECTOR<int,3>(ij,j,jk));subtriangles.Append(VECTOR<int,3>(ik,jk,k));subtriangles.Append(VECTOR<int,3>(jk,ik,ij));}
    else PHYSBAM_FATAL_ERROR(); // cannot have only 2 midpoints
    ARRAY<VECTOR<int,3> > candidate_triangle_list;
    bool subtriangle_in_triangle_mesh=false;
    for(int t=0;t<subtriangles.m;t++) if(Add_Triangle_Or_Subtriangles(candidate_triangle_list,triangle_mesh,segment_mesh,segment_midpoints,subtriangles(t))) subtriangle_in_triangle_mesh=true;
    if(subtriangle_in_triangle_mesh) triangle_list.Append_Elements(candidate_triangle_list);else triangle_list.Append(nodes);
    return subtriangle_in_triangle_mesh || triangle_mesh.Triangle(i,j,k);
}
//#####################################################################
// Function Initialize_Boundary_Mesh_With_T_Junctions
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Boundary_Mesh_With_T_Junctions(TRIANGLE_MESH& boundary_mesh_with_t_junctions,const ARRAY<int>& t_junctions,const ARRAY<VECTOR<int,2> >& t_junction_parents)
{
    assert(t_junctions.m==t_junction_parents.m);
    boundary_mesh_with_t_junctions.Clean_Memory();boundary_mesh_with_t_junctions.number_nodes=number_nodes;
    bool incident_tetrahedrons_defined=(incident_elements!=0);if(!incident_tetrahedrons_defined) Initialize_Incident_Elements();
    bool triangle_mesh_defined=(triangle_mesh!=0);if(!triangle_mesh_defined) Initialize_Triangle_Mesh();
    bool incident_triangles_defined=(triangle_mesh->incident_elements!=0);if(!incident_triangles_defined) triangle_mesh->Initialize_Incident_Elements();
    bool segment_mesh_defined=(triangle_mesh->segment_mesh!=0);if(!segment_mesh_defined) triangle_mesh->Initialize_Segment_Mesh();
    bool incident_segments_defined=(triangle_mesh->segment_mesh->incident_elements!=0);if(!incident_segments_defined) triangle_mesh->segment_mesh->Initialize_Incident_Elements();

    SEGMENT_MESH extended_segment_mesh(*triangle_mesh->segment_mesh);
    ARRAY<int> segment_midpoints(extended_segment_mesh.elements.m);
    for(int i=0;i<t_junctions.m;i++){
        int s=triangle_mesh->segment_mesh->Simplex(t_junction_parents(i));
        if(s>=0) segment_midpoints(s)=t_junctions(i);
        else{
            extended_segment_mesh.elements.Append(t_junction_parents(i));
            segment_midpoints.Append(t_junctions(i));
            extended_segment_mesh.number_nodes=max(extended_segment_mesh.number_nodes,t_junctions(i),t_junction_parents(i)(0),t_junction_parents(i)(1));}}
    extended_segment_mesh.Initialize_Incident_Elements();
    for(int tet=0;tet<elements.m;tet++) for(int e=0;e<4;e++){
        int i,j,k,l;
        switch(e){ // get edge (i,j,k) with proper orientation
          case 0: elements(tet).Get(i,k,j,l);break;
          case 1: elements(tet).Get(i,j,l,k);break;
          case 2: elements(tet).Get(i,l,k,j);break;
          case 3:default: elements(tet).Get(l,i,j,k);}
        if(Number_Of_Tetrahedrons_Across_Face(tet,i,j,k)==1) continue;
        Add_Triangle_Or_Subtriangles(boundary_mesh_with_t_junctions.elements,*triangle_mesh,extended_segment_mesh,segment_midpoints,VECTOR<int,3>(i,j,k));}

    HASHTABLE<VECTOR<int,3>,int> triangles_hash_table(boundary_mesh_with_t_junctions.elements.m);
    for(int tri=0;tri<boundary_mesh_with_t_junctions.elements.m;tri++){
        VECTOR<int,3>& triangle=boundary_mesh_with_t_junctions.elements(tri);
        VECTOR<int,3> sorted_triangle=triangle.Sorted();
        int tri2;
        if(!triangles_hash_table.Get(sorted_triangle,tri2))
            triangles_hash_table.Insert(sorted_triangle,tri);
        else{
            VECTOR<int,3>& triangle2=boundary_mesh_with_t_junctions.elements(tri2);
            if(triangle2.x){
                if(!TRIANGLE_MESH::Equivalent_Oriented_Triangles(triangle.x,triangle.y,triangle.z,triangle2.x,triangle2.y,triangle2.z)) triangle2.x=0; // the two triangles cancel out
                else PHYSBAM_FATAL_ERROR();}
            triangle.x=0;}} // mark triangle for deletion
    boundary_mesh_with_t_junctions.Delete_Elements_With_Missing_Nodes();

    if(!incident_segments_defined){delete triangle_mesh->segment_mesh->incident_elements;triangle_mesh->segment_mesh->incident_elements=0;}
    if(!segment_mesh_defined){delete triangle_mesh->segment_mesh;triangle_mesh->segment_mesh=0;}
    if(!incident_triangles_defined){delete triangle_mesh->incident_elements;triangle_mesh->incident_elements=0;}
    if(!triangle_mesh_defined){delete triangle_mesh;triangle_mesh=0;}
    if(!incident_tetrahedrons_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Initialize_Bending_Tetrahedrons
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Bending_Tetrahedrons(TRIANGLE_MESH& triangle_mesh)
{
    if(!triangle_mesh.adjacent_elements) triangle_mesh.Initialize_Adjacent_Elements();
    
    for(int t=0;t<triangle_mesh.elements.m;t++){
        VECTOR<int,3> nodes=triangle_mesh.elements(t);
        for(int a=0;a<(*triangle_mesh.adjacent_elements)(t).m;a++){
            int s=(*triangle_mesh.adjacent_elements)(t)(a);
            if(s>t){
                const VECTOR<int,3>& other_nodes=triangle_mesh.elements(s);
                if(other_nodes.Contains(nodes[0])){
                    cyclic_shift(nodes[0],nodes[1],nodes[2]);
                    if(other_nodes.Contains(nodes[0])) cyclic_shift(nodes[0],nodes[1],nodes[2]);}
                elements.Append(VECTOR<int,4>(nodes[0],nodes[1],nodes[2],triangle_mesh.Other_Node(nodes[1],nodes[2],s)));}}}
}
//#####################################################################
// Function Mark_Face_Connected_Component_Incident_On_A_Node
//#####################################################################
// implements a depth first search to find marked connected component
void TETRAHEDRON_MESH::
Mark_Face_Connected_Component_Incident_On_A_Node(const int node,const int tetrahedron_index_in_incident_elements,ARRAY<bool>& marked) const
{
    assert(incident_elements); // without this tetrahedron_index_in_incident_tetraherons makes no sense
    assert(adjacent_elements); // too expensive to do this every time, since this is a recursive function
    int tetrahedron=(*incident_elements)(node)(tetrahedron_index_in_incident_elements);
    if(marked.m != (*incident_elements)(node).m) marked.Resize((*incident_elements)(node).m);
    marked(tetrahedron_index_in_incident_elements)=true;
    for(int i=0;i<(*adjacent_elements)(tetrahedron).m;i++){
        int adjacent_tetrahedron=(*adjacent_elements)(tetrahedron)(i);
        if(Node_In_Tetrahedron(node,adjacent_tetrahedron)){int index;
            if((*incident_elements)(node).Find(adjacent_tetrahedron,index) && !marked(index)) Mark_Face_Connected_Component_Incident_On_A_Node(node,index,marked);}}
}
//#####################################################################
// Function Tetrahedrons_Across_Face
//#####################################################################
int TETRAHEDRON_MESH::
Tetrahedrons_Across_Face(const int tetrahedron,const int node1,const int node2,const int node3,ARRAY<int>& tetrahedrons_across_face) const
{
    assert(incident_elements);
    for(int k=0;k<(*incident_elements)(node1).m;k++){
        int tetrahedron2=(*incident_elements)(node1)(k); // tetrahedron in question
        if(tetrahedron2 != tetrahedron){
            int i,j,k,l;elements(tetrahedron2).Get(i,j,k,l);
            if(i == node2 && (j == node3 || k == node3 || l == node3)) tetrahedrons_across_face.Append(tetrahedron2);
            else if(j == node2 && (i == node3 || k == node3 || l == node3)) tetrahedrons_across_face.Append(tetrahedron2);
            else if(k == node2 && (i == node3 || j == node3 || l == node3)) tetrahedrons_across_face.Append(tetrahedron2);
            else if(l == node2 && (i == node3 || j == node3 || k == node3)) tetrahedrons_across_face.Append(tetrahedron2);}}
    return tetrahedrons_across_face.m;
}
//#####################################################################
// Function Identify_Connected_Components
//#####################################################################
// labels each face connected components with a unique id
void TETRAHEDRON_MESH::
Identify_Face_Connected_Components(ARRAY<int>& label)
{
    STACK<int> flood_fill_stack;
    flood_fill_stack.Preallocate(elements.m);
    bool adjacent_elements_defined=adjacent_elements!=0;
    if(!adjacent_elements_defined)Initialize_Adjacent_Elements();
    label.Resize(elements.m,init_all,0);
    int id=0;for(int t=0;t<elements.m;t++) if(!label(t)){
        id++;label(t)=id;
        flood_fill_stack.Push(t);
        while(!flood_fill_stack.Empty()){
            int top=flood_fill_stack.Pop();
            for(int i=0;i<(*adjacent_elements)(top).m;i++){
                int current_tetrahedron=(*adjacent_elements)(top)(i);
                if(!label(current_tetrahedron)){label(current_tetrahedron)=id;flood_fill_stack.Push(current_tetrahedron);}}}}
    if(!adjacent_elements_defined){delete adjacent_elements;adjacent_elements=0;}
}
//#####################################################################
// Function Identify_Edge_Node_Connected_Components
//#####################################################################
// labels each node connected components with a unique id
void TETRAHEDRON_MESH::
Identify_Edge_Connected_Components(ARRAY<int>& label)
{
    STACK<int> flood_fill_stack;flood_fill_stack.Preallocate(elements.m);
    bool neighbor_nodes_defined=neighbor_nodes!=0;if(!neighbor_nodes_defined)Initialize_Neighbor_Nodes();
    label.Resize(number_nodes,init_all,0);
    int id=0;for(int p=0;p<number_nodes;p++) if(!label(p)){
        id++;label(p)=id;flood_fill_stack.Push(p);
        while(!flood_fill_stack.Empty()){
            int top=flood_fill_stack.Pop();
            for(int n=0;n<(*neighbor_nodes)(top).m;n++){
                int current_node=(*neighbor_nodes)(top)(n);
                if(!label(current_node)){label(current_node)=id;flood_fill_stack.Push(current_node);}}}}  
    if(!neighbor_nodes_defined){delete neighbor_nodes;neighbor_nodes=0;}
}
//#####################################################################
// Function Assert_Consistent
//#####################################################################
bool TETRAHEDRON_MESH::
Assert_Consistent() const
{
    if(segment_mesh){assert(segment_mesh->number_nodes==number_nodes);assert(segment_mesh->Assert_Consistent());}
    if(triangle_mesh){assert(triangle_mesh->number_nodes==number_nodes);assert(triangle_mesh->Assert_Consistent());}
    if(element_edges){assert(segment_mesh);assert(element_edges->m==elements.m);}
    if(tetrahedron_faces){assert(triangle_mesh);assert(tetrahedron_faces->m==elements.m);}
    if(boundary_mesh){assert(boundary_mesh->number_nodes==number_nodes);assert(boundary_mesh->Assert_Consistent());}
    if(node_on_boundary) assert(node_on_boundary->m==number_nodes);
    if(edge_tetrahedrons){assert(segment_mesh);assert(edge_tetrahedrons->m==segment_mesh->elements.m);}
    if(triangle_tetrahedrons){assert(triangle_mesh);assert(triangle_tetrahedrons->m==triangle_mesh->elements.m);}
    return SIMPLEX_MESH<3>::Assert_Consistent();
}
//#####################################################################
// Function Set_Number_Nodes
//#####################################################################
void TETRAHEDRON_MESH::
Set_Number_Nodes(const int number_nodes_input)
{
    SIMPLEX_MESH<3>::Set_Number_Nodes(number_nodes_input);
    if(segment_mesh&&owns_segment_mesh) segment_mesh->Set_Number_Nodes(number_nodes); // TODO: clean this up once auxiliary structures are handled in a nice way
    if(triangle_mesh) triangle_mesh->Set_Number_Nodes(number_nodes);
    if(boundary_mesh) boundary_mesh->Set_Number_Nodes(number_nodes);
    if(node_on_boundary) node_on_boundary->Resize(number_nodes);
}
//#####################################################################
// Function Initialize_Swept_Mesh
//#####################################################################
void TETRAHEDRON_MESH::
Initialize_Swept_Mesh(const TRIANGLE_MESH& tri_mesh,int layers)
{
    for(int l=0,particle_offset=0;l<layers;l++){
        for(int e=0;e<tri_mesh.elements.m;e++){
            VECTOR<int,3> e0=tri_mesh.elements(e)+particle_offset,e1=e0+tri_mesh.number_nodes;
            int a=e0.Arg_Min(),c=e0.Arg_Max(),b=3-a-c;
            VECTOR<int,4> t0=e0.Append(e1(a)),t1(t0);
            t1(a)=e1(b);
            VECTOR<int,4> t2=t1;
            t2(b)=e1(c);
            exchange(t1(0),t1(1));
            elements.Append(t0);
            elements.Append(t1);
            elements.Append(t2);}
        particle_offset+=tri_mesh.number_nodes;}
    Set_Number_Nodes(tri_mesh.number_nodes*(layers+1));
}
//#####################################################################
