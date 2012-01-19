//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <cstring>
#include <fstream>
#include "COMPARISON_FUNCTIONS.h"

using namespace PhysBAM;
using namespace std;

#define levels 4

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{

    if (argc!=3){cout<<"Usage: project_triangles_towards_mouth input_mesh projection_volume"<<endl;exit(-1);}

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
    cout<<"Initializing triangle list..."<<endl;
    tet_vol.triangulated_surface->Update_Triangle_List();
    cout<<"Update vertex normals...."<<endl;
    tri_surf.Update_Vertex_Normals();
        
    cout<<"Initializing hierarchies..."<<endl;
    tet_vol.triangulated_surface->Initialize_Triangle_Hierarchy();

    cout<<"Initializing neighbor nodes..."<<endl;
    tri_mesh.Initialize_Neighbor_Nodes();

    cout<<"Initialize boundary mesh..."<<endl;
    tri_mesh.Initialize_Boundary_Mesh();
    cout<<"Initialize ordered loop segments..."<<endl;
    tri_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();

    cout<<tri_mesh.boundary_mesh->ordered_loop_nodes->m<<" ordered loop segments found."<<endl;

    tri_surf.avoid_normal_interpolation_across_sharp_edges=false;


    for(int i=1;i<=(*tri_mesh.boundary_mesh->ordered_loop_nodes).m;i++)
        cout<<i<<" "<<(*tri_mesh.boundary_mesh->ordered_loop_nodes)(i).m<<endl;

    //int smallest_i=ARRAY<ARRAY<int> >::argcomp(*tri_mesh.boundary_mesh->ordered_loop_nodes,compare_array_min_m<ARRAY<ARRAY<int> > >);
    //int smallest_m=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i).m;

    for(int smallest_i=0;smallest_i<2;smallest_i++){
    
        int smallest_m=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i).m;
        
        ARRAY<int> curr_level_nodes,next_level_nodes,visited_nodes;
        
        for(int i=0;i<smallest_m;i++) curr_level_nodes.Append((*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i)(i));

        ARRAY<float> weights;
        weights.Append(0);weights.Append(0.5);weights.Append(0.75);weights.Append(0.95);
        if(levels!=weights.m){cout<<"Bad weights..."<<endl;exit(-1);}

        for(int i=0;i<levels;i++){
            cout<<"Level "<<i<<" with "<<curr_level_nodes.m<<" nodes"<<endl;
            visited_nodes.Append_Elements(curr_level_nodes);
            for(int j=0;j<curr_level_nodes.m;j++){
                for(int k=1;k<=(*tri_mesh.neighbor_nodes)(curr_level_nodes(j)).m;k++)
                    if(!visited_nodes.Find((*tri_mesh.neighbor_nodes)(curr_level_nodes(j))(k)))next_level_nodes.Append((*tri_mesh.neighbor_nodes)(curr_level_nodes(j))(k));
                tri_particles.X(curr_level_nodes(j))=(float)weights(i)*tri_particles.X(curr_level_nodes(j))+(float)(1.0-weights(i))*
                    tet_vol.triangulated_surface->Oriented_Surface(tri_particles.X(curr_level_nodes(j)),(*tri_surf.vertex_normals)(1,curr_level_nodes(j)),0.001,0.001,0,0);}
            curr_level_nodes.Clean_Memory();curr_level_nodes=next_level_nodes;next_level_nodes.Clean_Memory();}
    }

    cout<<endl<<"Writing output mesh..."<<endl;
    FILE_UTILITIES::Write_To_File<float>("full_mesh_projected_mouth.tri",tri_surf);

    return 0;
} 

