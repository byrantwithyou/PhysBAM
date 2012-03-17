//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>

#define lower_cutoff .12

using namespace PhysBAM;
using namespace std;

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    if (argc!=4){cout<<"Usage: create_mesh_from_offsets tetrahedral_volume triangle_mesh offset_file"<<::endl;exit(-1);}

    TETRAHEDRON_MESH tet_mesh;
    TRIANGLE_MESH tri_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tet_particles;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tri_particles;
    TETRAHEDRALIZED_VOLUME<float> tet_vol(tet_mesh,tet_particles);
    TRIANGULATED_SURFACE<float> tri_surf(tri_mesh,tri_particles);

    cout<<"Reading files..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],tet_vol);
    FILE_UTILITIES::Read_From_File<float>(argv[2],tri_surf); 

    cout<<"Updating bounding boxes..."<<endl;
    tet_vol.Update_Bounding_Box();
    cout<<"Initializing triangulated surface..."<<endl;
    tet_vol.Initialize_Triangulated_Surface();
    cout<<"Updating triangle list..."<<endl;
    tet_vol.triangulated_surface->Update_Triangle_List();
    cout<<"Updating tetrahedron list..."<<endl;
    tet_vol.Update_Tetrahedron_List();

    ARRAY<PAIR<int,VECTOR_3D<float> > > offset_particles(tri_particles.number);

    cout<<"Reading coord file..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[3],offset_particles);cout<<"Processing..."<<endl;
    
    for(int p=0;p<tri_particles.number;p++)
        if(tri_particles.X(p).y>-lower_cutoff)tri_particles.X(p)=(*tet_vol.tetrahedron_list)(offset_particles(p).x).Point_From_Barycentric_Coordinates(offset_particles(p).y);

    FILE_UTILITIES::Write_To_File<float>("deformed_mesh_coord.tri",tri_surf);
    
    return 0;
}
//#################################################################
