//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
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
    if (argc!=3){cout<<"Usage: create correspondence tetrahedral_volume triangle_mesh"<<::endl;exit(-1);}
    
    TETRAHEDRON_MESH tet_mesh;
    TRIANGLE_MESH tri_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tet_particles;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tri_particles;
    TETRAHEDRALIZED_VOLUME<float> tet_vol(tet_mesh,tet_particles);
    TRIANGULATED_SURFACE<float> tri_surf(tri_mesh,tri_particles);

    cout<<"Reading files..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],tet_vol);
    FILE_UTILITIES::Read_From_File<float>(argv[2],tri_surf);

    cout<<"Updating tetrahedron list..."<<endl;
    tet_vol.Update_Tetrahedron_List();
    cout<<"Updating tetrahedron hierarcy..."<<endl;
    tet_vol.Initialize_Tetrahedron_Hierarchy();
    cout<<"Initializing triangulated surface..."<<endl;
    tet_vol.Initialize_Triangulated_Surface();
    cout<<"Updating triangle list..."<<endl;
    tet_vol.triangulated_surface->Update_Triangle_List();
    tri_surf.Update_Triangle_List();
    cout<<"Updating triangle hierarcy..."<<endl;
    tet_vol.triangulated_surface->Initialize_Triangle_Hierarchy();
    cout<<"Initialize incident tetrahedrons..."<<endl;
    tet_vol.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();
    cout<<"Initialize triangle mesh..."<<endl;
    tet_vol.tetrahedron_mesh.Initialize_Triangle_Mesh();
    cout<<"Initialize incident triangles..."<<endl;
    tet_vol.tetrahedron_mesh.triangle_mesh->Initialize_Incident_Triangles();
    tri_mesh.Initialize_Incident_Triangles();
    cout<<"Updating bounding boxes..."<<endl;
    tet_vol.Update_Bounding_Box();tri_surf.Update_Bounding_Box();
    cout<<"Updating vertex normals..."<<endl;
    tri_surf.Update_Vertex_Normals();

    cout<<"tet_vol.particles.number="<<tet_vol.particles.number<<endl;
    cout<<"tet_vol.tetrahedron_mesh.tetrahedrons.m="<<tet_vol.tetrahedron_mesh.tetrahedrons.m<<endl;
    cout<<"tet bounding_box size:"<<tet_vol.bounding_box->Edge_Lengths()<<endl;
    cout<<"tri bounding box size:"<<tri_surf.bounding_box->Edge_Lengths()<<endl;
    cout<<"tet surf triangles.m="<<tet_vol.triangulated_surface->triangle_mesh.triangles.m<<endl;
    cout<<"tri surf triangle.m="<<tri_surf.triangle_mesh.triangles.m<<endl;
    cout<<"tet surf particles.number="<<tet_particles.number<<endl;
    cout<<"tri surf particles.number="<<tri_particles.number<<endl;
    cout<<"Processing"<<flush;

    ARRAY<PAIR<int,VECTOR_3D<float> > > offset_particles(tri_particles.number);int nearest_triangle;

    tet_vol.triangulated_surface->avoid_normal_interpolation_across_sharp_edges=false;

    int one_hundredth=tri_particles.number/100;
    for(int p=0;p<tri_particles.number;p++){
        VECTOR_3D<float> proj=tet_vol.triangulated_surface->Oriented_Surface(tri_particles.X(p),(*tri_surf.vertex_normals)(1,p),(float)0.001,(float)0.001,&nearest_triangle,0);
        int i,j,k;tet_vol.triangulated_surface->triangle_mesh.triangles.Get(nearest_triangle,i,j,k);ARRAY<int> tets_on_face;tet_vol.tetrahedron_mesh.Tetrahedrons_On_Face(i,j,k,&tets_on_face);
        offset_particles(p).x=tets_on_face(1);offset_particles(p).y=(*tet_vol.tetrahedron_list)(tets_on_face(1)).Barycentric_Coordinates(tri_particles.X(p));
        tri_particles.X(p)=proj;if(p%one_hundredth==0)cout<<"."<<flush;}

    cout<<endl<<"Writing to file..."<<endl;FILE_UTILITIES::Write_To_File<float>("barycentrics.coord",offset_particles);
    FILE_UTILITIES::Write_To_File<float>("projected_surf.tri",tri_surf);

    return 0;
}
//#################################################################
