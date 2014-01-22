//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRALIZED_VOLUME
//##################################################################### 
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Geometry/Basic_Geometry/HEXAHEDRON.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> HEXAHEDRALIZED_VOLUME<T>::
HEXAHEDRALIZED_VOLUME()
    :MESH_OBJECT<TV,HEXAHEDRON_MESH>(*new HEXAHEDRON_MESH,*new GEOMETRY_PARTICLES<TV>),hexahedron_list(0),tetrahedralized_volume(0),triangulated_surface(0),hierarchy(0)
{
    this->need_destroy_mesh=this->need_destroy_particles=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> HEXAHEDRALIZED_VOLUME<T>::
HEXAHEDRALIZED_VOLUME(HEXAHEDRON_MESH& mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,HEXAHEDRON_MESH>(mesh_input,particles_input),hexahedron_list(0),tetrahedralized_volume(0),triangulated_surface(0),hierarchy(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> HEXAHEDRALIZED_VOLUME<T>::
~HEXAHEDRALIZED_VOLUME()
{
    Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Clean_Memory()
{
    MESH_OBJECT<TV,HEXAHEDRON_MESH>::Clean_Memory();
    delete hexahedron_list;hexahedron_list=0;
    delete tetrahedralized_volume;tetrahedralized_volume=0;
    delete triangulated_surface;triangulated_surface=0;
    delete hierarchy;hierarchy=0;
}
//#####################################################################
// Function Update_Hexahedron_List
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Update_Hexahedron_List()
{
    if(!hexahedron_list) hexahedron_list=new ARRAY<HEXAHEDRON<T> >(mesh.elements.m);
    for(int t=0;t<mesh.elements.m;t++){
        int p1,p2,p3,p4,p5,p6,p7,p8;mesh.elements(t).Get(p1,p2,p3,p4,p5,p6,p7,p8);
        (*hexahedron_list)(t).x0=particles.X(p1);(*hexahedron_list)(t).x1=particles.X(p2);
        (*hexahedron_list)(t).x2=particles.X(p3);(*hexahedron_list)(t).x3=particles.X(p4);
        (*hexahedron_list)(t).x4=particles.X(p5);(*hexahedron_list)(t).x5=particles.X(p6);
        (*hexahedron_list)(t).x6=particles.X(p7);(*hexahedron_list)(t).x7=particles.X(p8);}
}
//#####################################################################
// Funcion Initialize_Tetrahedralized_Volume
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Initialize_Tetrahedralized_Volume()
{
    mesh.Initialize_Faces();mesh.Initialize_Face_Hexahedrons();
    ARRAY<int> face_particle_indices(mesh.faces->m);ARRAY<int> hex_particle_indices(mesh.elements.m);ARRAY<VECTOR<int,4> > tetrahedron_list;
    //add node in the center of each hex
    for(int h=0;h<mesh.elements.m;h++){
        ARRAY<int> p(8);mesh.elements(h).Get(p(0),p(1),p(2),p(3),p(4),p(5),p(6),p(7));
        TV hex_center;for(int i=0;i<8;i++) hex_center+=particles.X(p(i));hex_center*=(T).125;
        particles.X(particles.Add_Element())=hex_center;hex_particle_indices(h)=particles.Size();}
    //add node in the center of each boundary face
    for(int f=0;f<mesh.faces->m;f++){
        int h=(*mesh.face_hexahedrons)(f)(1),node1,node2,node3,node4;(*mesh.faces)(f).Get(node1,node2,node3,node4);
        if(h<0){particles.X(particles.Add_Element())=(T).25*(particles.X(node1)+particles.X(node2)+particles.X(node3)+particles.X(node4));face_particle_indices(f)=particles.Size();}}
    //for each face, add in four tets from the associated octahedron
    for(int f=0;f<mesh.faces->m;f++){
        int h_outward,h_inward,h1,h2,p1,p2,p3,p4,ph_outward,ph_inward;(*mesh.faces)(f).Get(p1,p2,p3,p4);
        //find which hexahedron the face is outwardly oriented with
        h1=(*mesh.face_hexahedrons)(f)(0);h2=(*mesh.face_hexahedrons)(f)(1);h_outward=h1;h_inward=h2;
        if(h2){
            h_outward=h2;h_inward=h1;
            for(int k=0;k<6;k++){// loop over faces of h1
                int i1=mesh.face_indices[k][0],i2=mesh.face_indices[k][1],i3=mesh.face_indices[k][2],i4=mesh.face_indices[k][3];
                if(p1 == mesh.elements(h1)(i1) && p2 == mesh.elements(h1)(i2) &&
                    p3 == mesh.elements(h1)(i3) && p4 == mesh.elements(h1)(i4)){h_outward=h1;h_inward=h2;break;}}}
        ph_outward=hex_particle_indices(h_outward);if(h_inward == 0) ph_inward=face_particle_indices(f); else ph_inward=hex_particle_indices(h_inward);
        tetrahedron_list.Append(VECTOR<int,4>(p1,ph_inward,p2,ph_outward));tetrahedron_list.Append(VECTOR<int,4>(p2,ph_inward,p3,ph_outward));
        tetrahedron_list.Append(VECTOR<int,4>(p3,ph_inward,p4,ph_outward));tetrahedron_list.Append(VECTOR<int,4>(p4,ph_inward,p1,ph_outward));}
    if(tetrahedralized_volume) delete &(tetrahedralized_volume->mesh);delete tetrahedralized_volume;
    tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    tetrahedralized_volume->mesh.Initialize_Mesh(particles.Size(),tetrahedron_list);
}
//#####################################################################
// Funcion Initialize_Cube_Mesh_And_Particles
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Initialize_Cube_Mesh_And_Particles(const GRID<TV>& grid)
{
    particles.Delete_All_Elements();
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    mesh.Initialize_Cube_Mesh(m,n,mn);
    particles.Preallocate(m*n*mn);
    for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++)for(int i=0;i<m;i++) particles.X(particles.Add_Element())=grid.X(TV_INT(i,j,ij));
}
//#####################################################################
// Funcion Total_Volume
//#####################################################################
template<class T> T HEXAHEDRALIZED_VOLUME<T>::
Total_Volume() const
{
    T volume=0;
    for(int h=0;h<mesh.elements.m;h++){int p1,p2,p3,p4,p5,p6,p7,p8;mesh.elements(h).Get(p1,p2,p3,p4,p5,p6,p7,p8);
        volume+=HEXAHEDRON<T>::Signed_Volume(particles.X(p1),particles.X(p2),particles.X(p3),particles.X(p4),particles.X(p5),particles.X(p6),particles.X(p7),particles.X(p8));}
    return volume;
}
//#####################################################################
// Funcion Initialize_Triangulated_Surface
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Initialize_Triangulated_Surface()
{
    mesh.Initialize_Boundary_Mesh();
    triangulated_surface=new TRIANGULATED_SURFACE<T>(*mesh.boundary_mesh,particles);
}
//#####################################################################
namespace PhysBAM{
template class HEXAHEDRALIZED_VOLUME<float>;
template class HEXAHEDRALIZED_VOLUME<double>;
}
