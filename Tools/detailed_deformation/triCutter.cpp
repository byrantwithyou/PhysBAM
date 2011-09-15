//#####################################################################
// Copyright 2005-2006, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################r
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>

using namespace PhysBAM;
using namespace std;

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{

    if (argc!=2){cout<<"Useage triCutter mesh"<<endl;exit(-1);}
    
    TRIANGLE_MESH tri_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tri_particles;
    TRIANGULATED_SURFACE<float> tri_surf(tri_mesh,tri_particles);
    cout<<"Reading mesh..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],tri_surf);

    // Nose
    VECTOR_3D<float> center(0,0.01,0.13);VECTOR_3D<float> edges(0.04,0.03,0.04);      
    BOX_3D<float> box(center-(float).5*edges,center+(float).5*edges);   

    
    // Both eyes 
    //VECTOR_3D<float> center(0,0.04,0.10);VECTOR_3D<float> edges(0.22,0.024,0.05);      
    //BOX_3D<float> box(center-(float).5*edges,center+(float).5*edges);   
    
    // Right eye
    //VECTOR_3D<float> center(0.055,0.04,0.10);VECTOR_3D<float> edges(0.105,0.024,0.05);      
    //BOX_3D<float> box(center-(float).5*edges,center+(float).5*edges);   

    // Left eye
    //VECTOR_3D<float> center(-0.055,0.04,0.10);VECTOR_3D<float> edges(0.105,0.024,0.05);      
    //BOX_3D<float> box(center-(float).5*edges,center+(float).5*edges);   
    

    ARRAY<int> deletion_list;
    
    cout<<"Checking triangles..."<<endl;
    for(int i=1;i<=tri_mesh.triangles.m;i++){
        int j,k,l;tri_mesh.triangles.Get(i,j,k,l);
        if(!box.Inside(tri_particles.X(j))||!box.Inside(tri_particles.X(k))||!box.Inside(tri_particles.X(l))) deletion_list.Append(i);}

    cout<<"Triangles before deletion: "<<tri_mesh.triangles.m<<endl;
    tri_mesh.Delete_Triangles(deletion_list); 
    cout<<"Triangles after deletion: "<<tri_mesh.triangles.m<<endl;

    ARRAY<int> condensation;
    cout<<"Renumbering and discarding particles..."<<endl;
    tri_surf.Discard_Valence_Zero_Particles_And_Renumber(condensation);


    cout<<"Writing output..."<<endl;
    FILE_UTILITIES::Write_To_File<float>("cut_mesh.tri",tri_surf);
    FILE_UTILITIES::Write_To_File<float>("cut_condensation",condensation);

    return 0;

}
