//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <cstring>
#include <fstream>

#define big_number 100000

using namespace PhysBAM;
using namespace std;

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{

    if (argc!=4){cout<<"Usage: remove_cut_triangles small_mesh big_mesh big_color"<<endl;exit(-1);}

    TRIANGLE_MESH small_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > small_particles;
    TRIANGULATED_SURFACE<float> small_surf(small_mesh,small_particles);
    cout<<"Reading small mesh..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],small_surf);
    small_surf.Update_Triangle_List();
    small_surf.Initialize_Triangle_Hierarchy();

    TRIANGLE_MESH big_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > big_particles;
    TRIANGULATED_SURFACE<float> big_surf(big_mesh,big_particles);
    cout<<"Reading big mesh..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[2],big_surf);

    ARRAY<VECTOR_3D<float> > vertex_colors, alt_vertex_colors;
    cout<<"Reading color files..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[3],vertex_colors);

    big_surf.Update_Triangle_List();
    big_mesh.Initialize_Boundary_Mesh();
    big_mesh.boundary_mesh->Initialize_Connected_Segments();
    cout<<"Connected boundaries: "<<(*big_mesh.boundary_mesh->connected_segments).m<<" ";
    for(int i=0;i<(*big_mesh.boundary_mesh->connected_segments).m;i++)
        cout<<(*big_mesh.boundary_mesh->connected_segments)(i).m<<" ";
    cout<<endl; 

    VECTOR_3D<float> center(-0.001,-0.0285,0.08);VECTOR_3D<float> edges(0.055,0.0055,0.09);
    BOX_3D<float> box(center-(float).5*edges,center+(float).5*edges);

    ARRAY<int> big_particles_inside;
    ARRAY<int> small_particles_inside;

    cout <<"Finding interior particles"<<endl;
    for(int i=0;i<big_particles.number;i++)
        if(box.Inside(big_particles.X(i)))big_particles_inside.Append(i);
    for(int i=0;i<small_particles.number;i++)
        if(box.Inside(small_particles.X(i)))small_particles_inside.Append(i);
    
    ARRAY<int> big_to_small(big_particles.number);
   
    cout <<"Finding nearest corresponding particle for each small particle..."<<endl; //This process should use sorting... 
    for(int i=0;i<small_particles_inside.m;i++){
        int c_index=0;float c_dist=big_number;
        for(int j=0;j<big_particles_inside.m;j++){
            float new_dist=(small_particles.X(small_particles_inside(i))-big_particles.X(big_particles_inside(j))).Magnitude();
            if(new_dist<c_dist&&big_to_small(big_particles_inside(j))==0){c_index=j; c_dist=new_dist;}}
        assert(c_index); big_to_small(big_particles_inside(c_index))=i;}
    
    int cut=0;
    cout<<"Removing unassociated triangles..."<<endl;
    for(int i=big_mesh.triangles.m;i>0;i--){
        int j,k,l;big_mesh.triangles.Get(i,j,k,l);
        if(!box.Inside(big_particles.X(j))||!box.Inside(big_particles.X(k))||!box.Inside(big_particles.X(l)))continue;
        if(big_to_small(j)==0||big_to_small(k)==0||big_to_small(l)==0){big_mesh.triangles.Remove_Index_Lazy(i);cut++;}}

    ARRAY<int> condensation;
    big_surf.Discard_Valence_Zero_Particles_And_Renumber(condensation);
    cout<<"Condensing to "<<big_particles.number<<"..."<<endl;
    alt_vertex_colors.Resize(big_particles.number,true,true);
    for(int i=0;i<vertex_colors.m;i++)
        if(condensation(i)) alt_vertex_colors(condensation(i))=vertex_colors(i);
 
    big_mesh.Initialize_Adjacent_Triangles();
    big_mesh.Initialize_Boundary_Mesh();
    big_mesh.boundary_mesh->Initialize_Connected_Segments();
    cout<<"Connected boundaries: "<<(*big_mesh.boundary_mesh->connected_segments).m<<" ";
    for(int i=0;i<(*big_mesh.boundary_mesh->connected_segments).m;i++)
        cout<<(*big_mesh.boundary_mesh->connected_segments)(i).m<<" ";
    cout<<endl;

    cout<<"Post-processing nasty triangles..."<<endl;
    for(int i=big_mesh.triangles.m;i>0;i--){
        int j,k,l;big_mesh.triangles.Get(i,j,k,l);
        if(!box.Inside(big_particles.X(j))||!box.Inside(big_particles.X(k))||!box.Inside(big_particles.X(l)))continue;
        if((*big_mesh.adjacent_triangles)(i).m<2){big_mesh.triangles.Remove_Index_Lazy(i);cut++;}}

    condensation.Clean_Memory();
    big_surf.Discard_Valence_Zero_Particles_And_Renumber(condensation);
    cout<<"Condensing..."<<big_particles.number<<"..."<<endl;
    vertex_colors.Resize(big_particles.number,true,true);
    for(int i=0;i<alt_vertex_colors.m;i++)
        if(condensation(i)) vertex_colors(condensation(i))=alt_vertex_colors(i);

    big_mesh.Initialize_Boundary_Mesh();
    big_mesh.boundary_mesh->Initialize_Connected_Segments();
    cout<<"Connected boundaries: "<<(*big_mesh.boundary_mesh->connected_segments).m<<" ";
    for(int i=0;i<(*big_mesh.boundary_mesh->connected_segments).m;i++)
        cout<<(*big_mesh.boundary_mesh->connected_segments)(i).m<<" ";
    cout<<endl;

    cout<<"Cut "<<cut<<" triangles"<<endl;

    cout<<"Writing output..."<<endl;
    FILE_UTILITIES::Write_To_File<float>("new_removed.tri",big_surf);
    FILE_UTILITIES::Write_To_File<float>("new_color.col",vertex_colors);
    return 0;
} 

