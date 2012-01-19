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

#define levels 4
#define iterations 15 

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{

    if (argc!=4){cout<<"Usage: patch_eyes hires_tri tet_vol color_file"<<endl;exit(-1);}

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

    for(int i=0;i<4;i++)
        cout<<"Boundary loop ("<<i<<") has "<<(*tri_mesh.boundary_mesh->ordered_loop_nodes)(i).m<<" segments..."<<endl;

    TETRAHEDRON_MESH tet_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tet_particles;
    TETRAHEDRALIZED_VOLUME<float> tet_vol(tet_mesh,tet_particles);
    cout<<"Reading volume..."<<endl;        
    FILE_UTILITIES::Read_From_File<float>(argv[2],tet_vol);

    cout<<"Initializing triangulated surface..."<<endl;
    tet_vol.Initialize_Triangulated_Surface();

    ARRAY<int> deletion_list;
    
    VECTOR_3D<float> center(0,0.04,0.10);VECTOR_3D<float> edges(0.22,0.024,0.05);
    BOX_3D<float> box(center-(float).5*edges,center+(float).5*edges);  

    VECTOR_3D<float> left_center(-.055,0.04,0.10);VECTOR_3D<float> left_edges(0.10,0.024,0.05);
    BOX_3D<float> left_box(left_center-(float).5*left_edges,left_center+(float).5*left_edges);
    
    cout<<"Process bounding box..."<<endl;
    for(int i=0;i<tet_vol.triangulated_surface->triangle_mesh.triangles.m;i++){
        int j,k,l;tet_vol.triangulated_surface->triangle_mesh.triangles.Get(i,j,k,l);
        if((!left_box.Inside(tet_particles.X(j))||!left_box.Inside(tet_particles.X(k))||!left_box.Inside(tet_particles.X(l))))deletion_list.Append(i);}

    tet_vol.triangulated_surface->triangle_mesh.Delete_Triangles(deletion_list);

    deletion_list.Remove_All();
    tet_vol.triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();
    tet_vol.triangulated_surface->Refresh_Auxiliary_Structures();

    cout<<"Initializing hierarchies..."<<endl;
    tet_vol.triangulated_surface->Update_Triangle_List();
    tet_vol.triangulated_surface->Initialize_Triangle_Hierarchy();
    cout<<"Initializing adjacent triangles..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Initialize_Adjacent_Triangles();
    tri_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();
    tet_vol.triangulated_surface->triangle_mesh.Initialize_Neighbor_Triangles();

    ARRAY<int> triangle_marks(tet_vol.triangulated_surface->triangle_mesh.triangles.m);

    int smallest_i=1;
    
    cout<<"Processing boundary nodes..."<<endl;
    for(int i=1;i<=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i).m;i++){
        int node=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(smallest_i)(i),t;
        tet_vol.triangulated_surface->Oriented_Surface(tri_particles.X(node),(*tri_surf.vertex_normals)(1,node),0.1,0.1,&t,0);
        triangle_marks(t)=t;
        for(int j=1;j<=(*tet_vol.triangulated_surface->triangle_mesh.neighbor_triangles)(t).m;j++){
            triangle_marks((*tet_vol.triangulated_surface->triangle_mesh.neighbor_triangles)(t)(j))=(*tet_vol.triangulated_surface->triangle_mesh.neighbor_triangles)(t)(j);}}

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

    ARRAY<int> nodes;
    ARRAY<int> old_nodes;
    
    cout<<"Initialize clipped tet boundary loops segments..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Initialize_Boundary_Mesh();
    tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->Initialize_Connected_Segments();
    cout<<tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments->m<<" connected segments found."<<endl;
    
    for(int w=1;w<=(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments).m;w++)
        cout<<"The boundary segment loop("<<w<<") has "<<(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments)(w).m<<" segments..."<<endl;
    
    int largest_i=1;
    
    cout<<"Initialize incident triangles..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Initialize_Incident_Triangles();
    tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->Initialize_Connected_Segments();
    
    nodes.Clean_Memory();
    for(int i=1;i<=(*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments)(largest_i).m;i++)for(int j=0;j<2;j++)
        nodes.Append_Unique((*tet_vol.triangulated_surface->triangle_mesh.boundary_mesh->connected_segments)(largest_i)(j,i));
    
    deletion_list.Clean_Memory();
    cout<<"Processing nodes..."<<endl;
    for(int i=0;i<nodes.m;i++)for(int j=1;j<=(*tet_vol.triangulated_surface->triangle_mesh.incident_triangles)(nodes(i)).m;j++)
        deletion_list.Append((*tet_vol.triangulated_surface->triangle_mesh.incident_triangles)(nodes(i))(j));
    
    cout<<"Deleting..."<<endl;
    tet_vol.triangulated_surface->triangle_mesh.Delete_Triangles(deletion_list);
    cout<<"Refreshing..."<<endl;
    tet_vol.triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();
    
    tri_particles.Set_Array_Buffer_Size(1000);

    ARRAY<int> particle_map(tet_particles.number);
    VECTOR_3D<float> interior_color(.88,.67,.63);
    cout<<"Interior color is: "<<interior_color<<endl;
    cout<<"Adding tet triangles and particles to tri mesh..."<<endl;
    for(int i=0;i<tet_vol.triangulated_surface->triangle_mesh.triangles.m;i++){
        int j,k,l;tet_vol.triangulated_surface->triangle_mesh.triangles.Get(i,j,k,l);
        if(!particle_map(l)){
            tri_particles.Add_Particle();particle_map(l)=tri_particles.number;tri_particles.X(tri_particles.number)=tet_particles.X(l);tri_mesh.number_nodes++;old_nodes.Append_Unique(tri_particles.number);}
        if(!particle_map(k)){
            tri_particles.Add_Particle();particle_map(k)=tri_particles.number;tri_particles.X(tri_particles.number)=tet_particles.X(k);tri_mesh.number_nodes++;old_nodes.Append_Unique(tri_particles.number);}
        if(!particle_map(j)){
            tri_particles.Add_Particle();particle_map(j)=tri_particles.number;tri_particles.X(tri_particles.number)=tet_particles.X(j);tri_mesh.number_nodes++;old_nodes.Append_Unique(tri_particles.number);}
        tri_mesh.triangles.Append(particle_map(j),particle_map(k),particle_map(l));}

    ARRAY<VECTOR_3D<float> > vertex_append(tri_particles.number-vertex_colors.m);
    for(int i=0;i<vertex_append.m;i++)vertex_append(i)=interior_color;
    
    vertex_colors.Append_Elements(vertex_append);

    cout<<"Re-Initialize tri ordered loops..."<<endl;
    tri_mesh.Initialize_Boundary_Mesh();
    tri_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();

    //Bad stuff below
    int tri_i=1,tet_i=5;

    int tri_m=((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)).m;
    int tri_x_min_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i),compare_vector_min_x<VECTOR_3D<float> >);
    int tri_x_max_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i),compare_vector_max_x<VECTOR_3D<float> >);

    int tet_m=((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)).m;
    int tet_x_min_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i),compare_vector_min_x<VECTOR_3D<float> >);
    int tet_x_max_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i),compare_vector_max_x<VECTOR_3D<float> >);

    for(int w=1;w<=(*tri_mesh.boundary_mesh->ordered_loop_nodes).m;w++)
        cout<<"The boundary segment loop("<<w<<") has "<<(*tri_mesh.boundary_mesh->ordered_loop_nodes)(w).m<<" segments..."<<endl;
    
    float curr_pos,next_pos;int i=tri_x_min_i,curr_j=tet_x_min_i;
    
    do{
        int node1=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(i);
        int node2=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(i==1?tri_m:i-1);
        float avg=0.5*(tri_particles.X(node1).x+tri_particles.X(node2).x);
        curr_pos=tri_particles.X((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j)).x;
        next_pos=tri_particles.X((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j==tet_m?1:curr_j+1)).x;
        tri_mesh.triangles.Append(node1,node2,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j));
        if(fabs(avg-curr_pos)>fabs(avg-next_pos)||(next_pos<curr_pos)){
            int old_j=curr_j;
            curr_j=(curr_j==tet_m?1:curr_j+1);
            tri_mesh.triangles.Append(node2,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j),(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(old_j));}
        i=(i==1?tri_m:i-1);
    } while(i!=tri_x_max_i);

    i=tri_x_max_i,curr_j=tet_x_max_i;
    do{
        int node1=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(i);
        int node2=(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(i==1?tri_m:i-1);
        float avg=0.5*(tri_particles.X(node1).x+tri_particles.X(node2).x);
        curr_pos=tri_particles.X((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j)).x;
        next_pos=tri_particles.X((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j==tet_m?1:curr_j+1)).x;
        tri_mesh.triangles.Append(node1,node2,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j));
        if(fabs(avg-curr_pos)>fabs(avg-next_pos)||(curr_pos<next_pos)){
            int old_j=curr_j;
            curr_j=(curr_j==tet_m?1:curr_j+1);
            tri_mesh.triangles.Append(node2,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(curr_j),(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tet_i)(old_j));}
        i=(i==1?tri_m:i-1);
    } while(i!=tri_x_min_i);

    cout<<"Removing junk and cleaning..."<<endl;
    ARRAY<int> condensation;
    tri_surf.Discard_Valence_Zero_Particles_And_Renumber(condensation);
    
    cout<<"Re-Initialize tri ordered loops..."<<endl;
    tri_mesh.Initialize_Boundary_Mesh();
    tri_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();
    for(int w=1;w<=(*tri_mesh.boundary_mesh->ordered_loop_nodes).m;w++)
        cout<<"The boundary segment loop("<<w<<") has "<<(*tri_mesh.boundary_mesh->ordered_loop_nodes)(w).m<<" segments..."<<endl;

    //Bad stuff below
    tri_i=6;
    tri_mesh.triangles.Append((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(2),(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(1),
                                      (*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(3));

    
    tri_i=4;
    tri_x_min_i=ARRAY<VECTOR_3D<float> >::argcomp(tri_particles.X.array,(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i),compare_vector_max_x<VECTOR_3D<float> >);
    i=tri_x_min_i;int triangles=0;
    do{
        int next1=(i==7?1:i+1),next2=(i==6?1:(i==7)?2:i+2);
        tri_mesh.triangles.Append((*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(tri_x_min_i),(*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(next1),
                                          (*tri_mesh.boundary_mesh->ordered_loop_nodes)(tri_i)(next2));
        i=i%7+1;triangles++;
        if(triangles==5)break;
    }while(1);
    
    tri_mesh.Initialize_Neighbor_Nodes();
    nodes.Clean_Memory();
    cout<<"Getting new nodes for triangle old boundary"<<endl;
    for(int i=0;i<old_nodes.m;i++)nodes.Append(condensation(old_nodes(i)));
    int start=1,end=nodes.m;
    for(int j=0;j<levels;j++){
        for(int i=start;i<=end;i++)nodes.Append_Unique_Elements((*tri_mesh.neighbor_nodes)(nodes(i)));
        start=end;end=nodes.m;}

    cout<<"Smoothing geometry and colors..."<<endl; 
    ARRAY<VECTOR_3D<float> > new_positions(nodes.m);
    ARRAY<VECTOR_3D<float> > new_colors(nodes.m);
    for(int k=0;k<iterations;k++){
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
            tri_particles.X(nodes(i))=new_positions(i);vertex_colors(nodes(i))=new_colors(i);}}
    
    cout<<endl<<"Orienting..."<<endl;
    tri_mesh.Make_Orientations_Consistent();

    cout<<"Re-Initialize tri ordered loops..."<<endl;
    tri_mesh.Initialize_Boundary_Mesh();
    tri_mesh.boundary_mesh->Initialize_Ordered_Loop_Nodes();
    for(int w=1;w<=(*tri_mesh.boundary_mesh->ordered_loop_nodes).m;w++)
        cout<<"The boundary segment loop("<<w<<") has "<<(*tri_mesh.boundary_mesh->ordered_loop_nodes)(w).m<<" segments..."<<endl;

    cout<<endl<<"Writing output mesh..."<<endl;
    FILE_UTILITIES::Write_To_File<float>("patched_right_eye.tri",tri_surf);
    FILE_UTILITIES::Write_To_File<float>("patched_right_eye.col",vertex_colors);

    return 0;
}
