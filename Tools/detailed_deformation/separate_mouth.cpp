//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
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

    if (argc!=3){cout<<"Usage: separate_mouth input_mesh projection_volume"<<endl;exit(-1);}

    TRIANGLE_MESH tri_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tri_particles;
    TRIANGULATED_SURFACE<float> tri_surf(tri_mesh,tri_particles);
    cout<<"Reading mesh..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],tri_surf); 

    TETRAHEDRON_MESH tet_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tet_particles;
    TETRAHEDRALIZED_VOLUME<float> tet_vol(tet_mesh,tet_particles);
    cout<<"Reading volume..."<<endl;        
    FILE_UTILITIES::Read_From_File<float>(argv[2],tet_vol);

    cout<<"Initializing triangulated surface..."<<endl;
    tet_vol.Initialize_Triangulated_Surface();
    cout<<"Update vertex normals...."<<endl;
    tri_surf.Update_Vertex_Normals();
    cout<<"Initializing hierarchies..."<<endl;
    tet_vol.triangulated_surface->Update_Triangle_List();
    tet_vol.triangulated_surface->Initialize_Triangle_Hierarchy();

    tri_surf.avoid_normal_interpolation_across_sharp_edges=false;
    
    cout<<"Projecting particles..."<<endl;
    for(int i=1;i<tri_particles.number;i++)
    tri_particles.X(i)=tet_vol.triangulated_surface->Oriented_Surface(tri_particles.X(i),(*tri_surf.vertex_normals)(1,i),0.001,0.001,0,0);

    cout<<endl<<"Writing output mesh..."<<endl;
    FILE_UTILITIES::Write_To_File<float>("projected_mouth.tri",tri_surf);

    return 0;
} 

