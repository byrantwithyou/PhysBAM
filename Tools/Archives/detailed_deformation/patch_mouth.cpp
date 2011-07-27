//#####################################################################
// Copyright 2005-2006, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################r
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>
#include "COMPARISON_FUNCTIONS.h"

using namespace PhysBAM;
using namespace std;

#define levels 7
#define iterations 10

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{

    if (argc!=4){cout<<"Usage: patch_mouth hires_tri tet_vol color_file"<<endl;exit(-1);}

    ARRAY<VECTOR_3D<float> > vertex_colors;
    cout<<"Reading color file..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[3],vertex_colors);

    TRIANGLE_MESH tri_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tri_particles;
    TRIANGULATED_SURFACE<float> tri_surf(tri_mesh,tri_particles);
    cout<<"Reading mesh..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],tri_surf); 

    cout<<"Initialize boundary mesh..."<<endl;
    tri_mesh.Initialize_Boundary_Mesh();
    cout<<"Initialize connected segments..."<<endl;
    tri_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();
    cout<<tri_mesh.boundary_mesh->ordered_loop_nodes->m<<" connected segments found."<<endl;
    cout<<"Update vertex normals...."<<endl;
    tri_surf.Update_Vertex_Normals();

    int smallest_i=ARRAY<ARRAY<int> >::argcomp(*tri_mesh.boundary_mesh->ordered_loop_nodes,compare_array_min_m<ARRAY<ARRAY<int> > >);
    int smallest_m=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i).m;
    cout<<"The smallest boundary loop ("<<smallest_i<<") has "<<smallest_m<<" segments..."<<endl;

    TETRAHEDRON_MESH tet_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tet_particles;
    TETRAHEDRALIZED_VOLUME<float> tet_vol(tet_mesh,tet_particles);
    cout<<"Reading volume..."<<endl;        
    FILE_UTILITIES::Read_From_File<float>(argv[2],tet_vol);

    cout<<"Initializing triangulated surface..."<<endl;
    tet_vol.Initialize_Triangulated_Surface();

    ARRAY<int> deletion_list;

    VECTOR_3D<float> center(0,-0.0299357,0.12);VECTOR_3D<float> edges(0.2,0.06,0.15);      
    BOX_3D<float> box(center-(float).5*edges,center+(float).5*edges);

    cout<<"Process bounding box..."<<endl;
    for(int i=1;i<=tet_vol.triangulated_surface->triangle_mesh.triangles.m;i++){
        int j,k,l;tet_vol.triangulated_surface->triangle_mesh.triangles.Get(i,j,k,l);
        if(!box.Inside(tet_particles.X(j))||!box.Inside(tet_particles.X(k))||!box.Inside(tet_particles.X(l))) deletion_list.Append(i);}


    deletion_list.Remove_All();
    tet_vol.triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();
    tet_vol.triangulated_surface->Refresh_Auxiliary_Structures();
   
    cout<<"Initializing hierarchies..."<<endl;
    tet_vol.triangulated_surface->Update_Triangle_List();
    tet_vol.triangulated_surface->Initialize_Triangle_Hierarchy();
    cout<<"Initializing adjacent triangles..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Initialize_Adjacent_Triangles();

    ARRAY<int> triangle_marks(tet_vol.triangulated_surface->triangle_mesh.triangles.m);
    
    cout<<"Processing boundary nodes..."<<endl;
    for(int i=1;i<=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i).m;i++){
        int node=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i)(i),t;
        tet_vol.triangulated_surface->Oriented_Surface(tri_particles.X(node),(*tri_surf.vertex_normals)(1,node),0.001,0.001,&t,0);
        triangle_marks(t)=t;}

    cout<<"Walking mesh..."<<endl;
    ARRAY<int> queue;queue.Append_Unique(1);
    while(queue.m){
        int curr_triangle;queue.Pop(curr_triangle);
        deletion_list.Append(curr_triangle);
        for(int i=1;i<=(*tet_vol.triangulated_surface->triangle_mesh.adjacent_triangles)(curr_triangle).m;i++){
            int cand_triangle=(*tet_vol.triangulated_surface->triangle_mesh.adjacent_triangles)(curr_triangle)(i);
            if(!triangle_marks(cand_triangle)&&!deletion_list.Find(cand_triangle))queue.Append_Unique(cand_triangle);}}
    
    cout<<"Deleting and refreshing..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Delete_Triangles(deletion_list);
    tet_vol.triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();

    cout<<"Initialize clipped tet boundary mesh..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Initialize_Boundary_Mesh();
    cout<<"Initialize clipped tet boundary loops segments..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->Initialize_Connected_Segments();
    cout<<tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments->m<<" connected segments found."<<endl;

    smallest_i=ARRAY<ARRAYS<int> >::argcomp(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments,compare_array_max_m<ARRAY<ARRAYS<int> > >);
    smallest_m=(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments)(smallest_i).m;
    cout<<"The largest boundary segment loop("<<smallest_i<<") has "<<smallest_m<<" segments..."<<endl;


    deletion_list.Remove_All();
    cout<<"Initialize incident triangles..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Initialize_Incident_Triangles();

    ARRAY<int> nodes;
    for(int i=1;i<=(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments)(smallest_i).m;i++)for(int j=1;j<=2;j++)
    nodes.Append_Unique((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments)(smallest_i)(j,i));

    cout<<"Processing nodes..."<<endl;
    for(int i=1;i<=nodes.m;i++)for(int j=1;j<=(*tet_vol.triangulated_surface->triangle_mesh.incident_triangles)(nodes(i)).m;j++)
            deletion_list.Append((*tet_vol.triangulated_surface->triangle_mesh.incident_triangles)(nodes(i))(j));
    
    cout<<"Deleting and refreshing..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Delete_Triangles(deletion_list);
    tet_vol.triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();

    cout<<"Re-Initialize clipped tet boundary mesh..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Initialize_Boundary_Mesh();
    cout<<"Re-Initialize clipped tet ordered loops..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();
    cout<<tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes->m<<" connected segments found."<<endl;

    int tet_i=ARRAY<ARRAY<int> >::argcomp(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes,compare_array_min_m<ARRAY<ARRAY<int> > >);
    int tet_m=(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i).m;
    cout<<"The tet boundary loop ("<<tet_i<<") has "<<tet_m<<" segments..."<<endl;

    cout<<"Finding tet edge vertices..."<<endl;
    int tet_x_min_i=ARRAY<VECTOR_3D<float> >::argcomp(tet_particles.X.array,(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i),
        compare_vector_min_x<VECTOR_3D<float> >);
    int tet_x_max_i=ARRAY<VECTOR_3D<float> >::argcomp(tet_particles.X.array,(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i),
        compare_vector_max_x<VECTOR_3D<float> >);

    cout<<tet_x_min_i<<" "<<tet_particles.X((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(tet_x_min_i))<<endl;
    cout<<tet_x_max_i<<" "<<tet_particles.X((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(tet_x_max_i))<<endl;
    
    cout<<"Re-Initialize tri boundary mesh..."<<endl;
    tri_mesh.Initialize_Boundary_Mesh();
    cout<<"Re-Initialize tri ordered loops..."<<endl;
    tri_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();
    cout<<tri_mesh.boundary_mesh->ordered_loop_nodes->m<<" connected segments found."<<endl;

    int tri_i=ARRAY<ARRAY<int> >::argcomp(*tri_mesh.boundary_mesh->ordered_loop_nodes,compare_array_min_m<ARRAY<ARRAY<int> > >);
    int tri_m=((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)).m;
    cout<<"The smallest boundary loop ("<<tri_i<<") has "<<tri_m<<" segments..."<<endl;

    cout<<"Finding tri edge vertices..."<<endl;    
    int tri_x_min_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i),compare_vector_min_x<VECTOR_3D<float> >);
    int tri_x_max_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i),compare_vector_max_x<VECTOR_3D<float> >);

    cout<<tri_x_min_i<<" "<<tri_particles.X((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(tri_x_min_i))<<endl;
    cout<<tri_x_max_i<<" "<<tri_particles.X((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(tri_x_max_i))<<endl;

    ARRAY<int>triangle_map(tet_vol.triangulated_surface->triangle_mesh.triangles.m),particle_map(tet_particles.number);
    
    VECTOR_3D<float> interior_color=vertex_colors((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)((tri_x_min_i+tri_m/2)%tri_m));
    cout<<"Interior color is: "<<interior_color<<endl;
    cout<<"Adding tet triangles and particles to tri mesh..."<<endl;
    for(int i=1;i<=tet_vol.triangulated_surface->triangle_mesh.triangles.m;i++){
        int j,k,l;tet_vol.triangulated_surface->triangle_mesh.triangles.Get(i,j,k,l);
        if(!particle_map(l)){
            tri_particles.Add_Particle();particle_map(l)=tri_particles.number;tri_particles.X(tri_particles.number)=tet_particles.X(l);tri_mesh.number_nodes++;vertex_colors.Append(interior_color);}
        if(!particle_map(k)){
            tri_particles.Add_Particle();particle_map(k)=tri_particles.number;tri_particles.X(tri_particles.number)=tet_particles.X(k);tri_mesh.number_nodes++;vertex_colors.Append(interior_color);}
        if(!particle_map(j)){
            tri_particles.Add_Particle();particle_map(j)=tri_particles.number;tri_particles.X(tri_particles.number)=tet_particles.X(j);tri_mesh.number_nodes++;vertex_colors.Append(interior_color);}
        tri_mesh.triangles.Append(particle_map(j),particle_map(k),particle_map(l));triangle_map(i)=tri_mesh.triangles.m;}

    cout<<"Connecting triangles..."<<endl; 
    ARRAY<int> old_nodes;
    int curr_j=tet_x_min_i,i=tri_x_min_i;float curr_pos,next_pos;
    do{
        int node1=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(i);
        int node2=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)((i==1?tri_m:i-1));
        old_nodes.Append_Unique(node1);old_nodes.Append_Unique(node2);
        float avg=0.5*(tri_particles.X(node1).x+tri_particles.X(node2).x);
        curr_pos=tet_particles.X((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j)).x;
        next_pos=tet_particles.X((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j%tet_m+1)).x;
        tri_mesh.triangles.Append(node1,node2,particle_map((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j)));
        if(fabs(avg-curr_pos)>fabs(avg-next_pos)){           
            int old_j=curr_j;curr_j=curr_j%tet_m+1;
            tri_mesh.triangles.Append(node2,particle_map((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j)),
                                              particle_map((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(old_j)));}
        i=(i==1?tri_m:i-1);
    } while(i!=tri_x_min_i);
    
    cout<<"Removing junk and cleaning..."<<endl;
    ARRAY<int> condensation;
    tri_surf.Discard_Valence_Zero_Particles_And_Renumber(condensation);
    tri_mesh.Initialize_Neighbor_Nodes();

    nodes.Clean_Memory();
    cout<<"Getting new nodes for triangle old boundary"<<endl;
    for(int i=1;i<=old_nodes.m;i++)nodes.Append(condensation(old_nodes(i)));
    int start=1,end=nodes.m;
    for(int j=1;j<=levels;j++){
        for(int i=start;i<=end;i++) nodes.Append_Unique_Elements((*tri_mesh.neighbor_nodes)(nodes(i)));
        start=end;end=nodes.m;}

    cout<<"Smoothing geometry and colors..."<<endl; 
    ARRAY<VECTOR_3D<float> > new_positions(nodes.m);
    ARRAY<VECTOR_3D<float> > new_colors(nodes.m);
    for(int k=1;k<=iterations;k++){
        cout<<"Iteration "<<k<<endl;
        for(int i=1;i<nodes.m;i++){
            new_colors(i)=VECTOR_3D<float>();new_positions(i)=VECTOR_3D<float>();float total_weight=0;
            for(int j=1;j<=(*tri_mesh.neighbor_nodes)(nodes(i)).m;j++){
                new_positions(i)+=tri_particles.X((*tri_mesh.neighbor_nodes)(nodes(i))(j));
                float weight=(tri_particles.X((*tri_mesh.neighbor_nodes)(nodes(i))(j))-tri_particles.X(nodes(i))).Magnitude();
                new_colors(i)+=weight*vertex_colors((*tri_mesh.neighbor_nodes)(nodes(i))(j));total_weight+=weight;}
            new_colors(i)*=(1.0/total_weight);
            new_positions(i)*=(1.0/(*tri_mesh.neighbor_nodes)(nodes(i)).m);}
        for(int i=1;i<nodes.m;i++){
            tri_particles.X(nodes(i))=new_positions(i);vertex_colors(nodes(i))=new_colors(i);}
    }

    cout<<endl<<"Writing output mesh..."<<endl;
    FILE_UTILITIES::Write_To_File<float>("patched_mouth.tri",tri_surf);
    FILE_UTILITIES::Write_To_File<float>("patched_mouth.col",vertex_colors);

    return 0;
}
