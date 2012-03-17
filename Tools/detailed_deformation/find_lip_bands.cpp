//#####################################################################
// Copyright 2005-2006, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <cstring>
#include <fstream>
#include "COMPARISON_FUNCTIONS.h"

using namespace PhysBAM;
using namespace std;


//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    if (argc!=3&&argc!=2){cout<<"Usage: patch_mouth unpatched_mesh [bot]"<<endl;exit(-1);}

    TRIANGLE_MESH tri_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tri_particles;
    TRIANGULATED_SURFACE<float> tri_surf(tri_mesh,tri_particles);
    cout<<"Reading mesh..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],tri_surf); 

    cout<<"Initializing incident triangles..."<<endl;
    tri_mesh.Initialize_Incident_Triangles();

    cout<<"Initialize boundary mesh..."<<endl;
    tri_mesh.Initialize_Boundary_Mesh();
    cout<<"Initialize ordered loop segments..."<<endl;
    tri_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();
    cout<<tri_mesh.boundary_mesh->ordered_loop_nodes->m<<" ordered loop segments found."<<endl;

    int smallest_i=ARRAY<ARRAY<int> >::argcomp(*tri_mesh.boundary_mesh->ordered_loop_nodes,compare_array_min_m<ARRAY<ARRAY<int> > >);
    int smallest_m=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i).m;
    cout<<"The smallest boundary loop ("<<smallest_i<<") has "<<smallest_m<<" segments..."<<endl;

    cout<<"Finding edge vertices..."<<endl;
    int x_min_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i),compare_vector_min_x<VECTOR_3D<float> >);
    int x_max_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i),compare_vector_max_x<VECTOR_3D<float> >);

    cout<<x_min_i<<" "<<tri_particles.X((*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i)(x_min_i))<<endl;
    cout<<x_max_i<<" "<<tri_particles.X((*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i)(x_max_i))<<endl;
    
    ARRAY<int> deletion_list;
    for(int i=0;i<tri_mesh.triangles.m;i++){deletion_list.Append(i);}

    string out_name="lip_band_";
    
    if(argc==3){
    for(int i=(x_max_i==smallest_m?1:x_max_i+1);i!=x_min_i;i=i%smallest_m+1){
        int node=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i)(i);
        for(int j=0;j<(*tri_mesh.incident_triangles)(node).m;j++){int index=deletion_list.Find((*tri_mesh.incident_triangles)(node)(j));
        if(index)deletion_list.Remove_Index(index);}}out_name+="bot";}
    else{cout<<"Processing top edge from left to right..."<<endl;
    for(int i=(x_min_i==smallest_m?1:x_min_i+1);i!=x_max_i;i=i%smallest_m+1){
        int node=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i)(i);
        for(int j=0;j<(*tri_mesh.incident_triangles)(node).m;j++){int index=deletion_list.Find((*tri_mesh.incident_triangles)(node)(j));
        if(index)deletion_list.Remove_Index(index);}}out_name+="top";}

    tri_mesh.Delete_Triangles(deletion_list);
    tri_surf.Refresh_Auxiliary_Structures();
    tri_surf.Discard_Valence_Zero_Particles_And_Renumber();

    cout<<"Writing output mesh..."<<endl;
    FILE_UTILITIES::Write_To_File<float>(out_name+".tri",tri_surf);

    cout<<"Updating bounding box..."<<endl;
    tri_surf.Update_Bounding_Box();

    cout<<"Flattening..."<<endl;
    for(int i=0;i<tri_particles.number;i++)tri_particles.X(i).y=tri_surf.bounding_box->Center().y;
    tri_surf.Refresh_Auxiliary_Structures();

    cout<<"Writing flat output mesh..."<<endl;
    FILE_UTILITIES::Write_To_File<float>(out_name+"_flat.tri",tri_surf);

    return 0;
} 

