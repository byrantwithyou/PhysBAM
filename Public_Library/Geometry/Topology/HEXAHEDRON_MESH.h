//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ron Fedkiw, Neil Molino, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRON_MESH
//#####################################################################
#ifndef __HEXAHEDRON_MESH__
#define __HEXAHEDRON_MESH__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{ 

class SEGMENT_MESH;
class TRIANGLE_MESH;

class HEXAHEDRON_MESH
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef VECTOR<int,8> ELEMENT_TYPE;
    enum {dimension=3};

    int number_nodes; // number of nodes in the mesh
    ARRAY<VECTOR<int,8> > elements; // array of 8 indices for each hexahedron - hexahedron(j,i) is j'th index in hexahedron i
    ARRAY<ARRAY<int> >* incident_elements; // for each node, list of neighboring hexahedrons that contain it
    ARRAY<ARRAY<int> >* adjacent_elements; // for each hexahedron, list of (up to 6) adjacent neighboring hexahedrons
    ARRAY<VECTOR<int,4> >* faces; // 4 indices for each face in the mesh
    ARRAY<bool>* node_on_boundary; // node_on_boundary(i) is true if node i is on the boundary
    ARRAY<int>* boundary_nodes;
    ARRAY<VECTOR<int,2> >* face_hexahedrons; // for each face, list of incident hexahedrons
    TRIANGLE_MESH* boundary_mesh;
    static const int face_indices[6][4];
    static const int edge_indices[12][2];

    HEXAHEDRON_MESH();
    HEXAHEDRON_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,8> >& hexahedron_list);
    ~HEXAHEDRON_MESH();

private:
    void operator=(const HEXAHEDRON_MESH&); // use Initialize_Mesh instead
public:

    void Clean_Memory()
    {elements.Clean_Memory();Delete_Auxiliary_Structures();}

    void Refresh_Auxiliary_Structures()
    {if(incident_elements) Initialize_Incident_Elements();if(adjacent_elements) Initialize_Adjacent_Elements();
    if(faces) Initialize_Faces();
    if(node_on_boundary) Initialize_Node_On_Boundary();
    if(boundary_mesh) Initialize_Boundary_Mesh();
    if(face_hexahedrons) Initialize_Face_Hexahedrons();}

    void Initialize_Mesh(const int number_nodes_input,const ARRAY<VECTOR<int,8> >& hexahedron_list)
    {Clean_Memory();number_nodes=number_nodes_input;elements=hexahedron_list;}

    void Initialize_Mesh_With_Particle_Offset(const HEXAHEDRON_MESH& mesh,const int particle_offset)
    {elements.Resize(mesh.elements.m,no_init);
    for(int t=0;t<elements.m;t++) elements(t)=mesh.elements(t)+particle_offset;
    number_nodes=0;} // TODO: currently leaves number_nodes uninitialized

    bool Node_In_Hexahedron(const int node,const int hexahedron) const
    {int p1,p2,p3,p4,p5,p6,p7,p8;elements(hexahedron).Get(p1,p2,p3,p4,p5,p6,p7,p8);
    return p1 == node || p2 == node || p3 == node || p4 == node || p5 == node || p6 == node || p7 == node || p8 == node;}

    void Initialize_Cube_Mesh(const int m,const int n,const int p)
    {Clean_Memory();number_nodes=m*n*p;elements.Resize((m-1)*(n-1)*(p-1));
    int t=0;for(int k=0;k<p-1;k++) for(int j=0;j<n-1;j++) for(int i=0;i<m-1;i++){t++;
        elements(t)(0)=i+m*n*(k-1)+m*(j-1);elements(t)(1)=i+m*n*k+m*(j-1);elements(t)(2)=i+m*n*(k-1)+m*j;elements(t)(3)=i+m*n*k+m*j;elements(t)(4)=(i+1)+m*n*(k-1)+m*(j-1);
        elements(t)(5)=(i+1)+m*n*k+m*(j-1);elements(t)(6)=(i+1)+m*n*(k-1)+m*j;elements(t)(7)=(i+1)+m*n*k+m*j;}}

    template<class T_CONNECTIVITY> void Add_Connectivity(T_CONNECTIVITY& connectivity) const
    {for(int t=0;t<elements.m;t++) connectivity.Union(elements(t));}

    template<class RW>
    void Read(std::istream& input)
    {Clean_Memory();int backward_compatible;Read_Binary<RW>(input,number_nodes,backward_compatible);Read_Binary<RW>(input,elements);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,number_nodes,8);Write_Binary<RW>(output,elements);}

//#####################################################################
    void Delete_Auxiliary_Structures();
    void Initialize_Neighbor_Nodes();
    void Initialize_Incident_Elements();
    void Initialize_Adjacent_Elements();
private:
    void Find_And_Append_Adjacent_Elements(const int hexahedron,const int node1,const int node2,const int node3,const int node4);
public:
    void Initialize_Faces();
    void Initialize_Node_On_Boundary();
    void Initialize_Boundary_Nodes();
    int Number_Of_Hexahedrons_Across_Face(const int hexahedron,const int node1,const int node2,const int node3,const int node4) const;
    int Delete_Hexahedrons_With_Missing_Nodes(); // returns the number deleted
    void Delete_Hexahedrons(const ARRAY<int>& deletion_list);
    void Initialize_Boundary_Mesh();
    void Initialize_Face_Hexahedrons();
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const;
    void Set_Number_Nodes(const int number_nodes_input);
    void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const;
//#####################################################################
};   
}
#endif
