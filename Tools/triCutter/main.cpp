//#####################################################################
// Copyright 2005, Melody Wu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>
using namespace PhysBAM;
using namespace std;

//#################################################################
// Function MakeCuts
//#################################################################
template<class T> void MakeCuts(std::string tri_filename,std::string new_tri_filename, std::string color_filename)
{
    BOX_3D<float> box = BOX_3D<float>(-0.02,0.022,-0.01,0.018,0.0,0.15);
    //BOX_3D<float> box = BOX_3D<float>(-0.10,0.11,-0.17,-0.1,-0.11,0.15);
    ARRAY<int> Outside_Triangles;
    ARRAY<int> condensation_mapping;
    int polygon_count = 0;
    ARRAY<VECTOR_3D<T> > old_vertex_color_array;
    ARRAY<VECTOR_3D<T> > new_vertex_color_array;
    
    TRIANGULATED_SURFACE<T>* triangulated_surface=0;
    //TRIANGULATED_SURFACE<T>* plane_triangulated_surface=0;
    FILE_UTILITIES::Create_From_File<T>(tri_filename,triangulated_surface);
    //FILE_UTILITIES::Create_From_File<T>("ground.tri.gz",plane_triangulated_surface);
    std::cout<<"done reading from file"<<std::endl;
/*
    triangulated_surface->triangle_mesh.Initialize_Segment_Mesh();
    triangulated_surface->triangle_mesh.Initialize_Boundary_Mesh();
    triangulated_surface->triangle_mesh.boundary_mesh->Initialize_Connected_Segments();
    SEGMENT_MESH* boundary_mesh = triangulated_surface->triangle_mesh.boundary_mesh;
    ARRAY<ARRAYS<int> >* connected_segments = boundary_mesh->connected_segments;
    int largest_array = 1, largest_size = 0;
    std::cout<<"largest array = "<<largest_array<<" size = "<<largest_size<<std::endl;
    std::cout<<"size of connected_segments = "<<connected_segments->Max_Size()<<std::endl;
    for(int i=1; i<=connected_segments->m; i++){
        ARRAYS<int> segment_array;
        connected_segments->Get(i,segment_array);
        std::cout<<"i = "<<i<<" size = "<<segment_array.m<<std::endl;
        if(segment_array.m > largest_size){
            largest_size = segment_array.m;
            largest_array = i;}}
    std::cout<<"largest array = "<<largest_array<<" size = "<<largest_size<<std::endl;
    ARRAYS<int> largest_segment_array;
    connected_segments->Get(largest_array,largest_segment_array);
    FILE_UTILITIES::Write_To_File<T>("neck_boundary_segments",largest_segment_array);
    
    RANDOM_NR3 random;
    random.Set_Seed(-1);

    float minimum_sum_of_squares = 1000000.0;
    VECTOR_3D<T> best_normal_vector = VECTOR_3D<T>(-0.0314169,0.980131,0.195845);
    //PLANE<T> best_plane;
    float d;

    //for(int i=0; i<1000000; i++){
    //VECTOR_3D<T> normal = random.Get_Uniform_Vector(VECTOR_3D<T>(-0.2,0,0),VECTOR_3D<T>(0.2,1,0.5));//VECTOR_3D<T>(0,1,0);
    //normal = normal.Normalized();
    //std::cout<<"normal vector= "<<normal<<std::endl;
    VECTOR_3D<T> point_on_plane;
    VECTOR_3D<T> point_on_curve;
    d = 100.0;
    //find d
    for(int i=1; i<largest_segment_array.m; i++){
        int vertex1, vertex2;
        largest_segment_array.Get(i,vertex1,vertex2);
        VECTOR_3D<T> point1 = triangulated_surface->particles.X(vertex1);
        VECTOR_3D<T> point2 = triangulated_surface->particles.X(vertex2);
        point_on_curve = point1;
        T dot_product1 = point1.Dot_Product(point1,best_normal_vector);
        T dot_product2 = point2.Dot_Product(point2,best_normal_vector);
        
        if(dot_product1<d){
            d=dot_product1;
            point_on_plane = point1 - best_normal_vector*0.001;}
        if(dot_product2<d){
            d=dot_product2;
            point_on_plane = point2 - best_normal_vector*0.001;}}
    
    //make plane
    PLANE<T> best_plane = PLANE<T>(best_normal_vector,point_on_plane);
*/
    //generate a measure of how good the plane is - sum of squares?
    /*float this_plane_sum_of_squares = 0.0;
    for(int i=1; i<largest_segment_array.m; i++){
        int vertex1, vertex2;
        largest_segment_array.Get(i,vertex1,vertex2);
        VECTOR_3D<T> point1 = triangulated_surface->particles.X(vertex1);
        VECTOR_3D<T> point2 = triangulated_surface->particles.X(vertex2);
        float distance1 = (point1 - plane.Surface(point1)).Magnitude_Squared();
        float distance2 = (point2 - plane.Surface(point2)).Magnitude_Squared();
        this_plane_sum_of_squares = this_plane_sum_of_squares + distance1 + distance2;
    }
    if(this_plane_sum_of_squares<minimum_sum_of_squares){
        minimum_sum_of_squares = this_plane_sum_of_squares;
        std::cout<<"min sum of squares = "<<this_plane_sum_of_squares<<std::endl;
        best_normal_vector = normal;
        best_plane = plane;
    }
}

    std::cout<<"best normal vector = "<<best_normal_vector<<std::endl;
    //std::cout<<"minimum sum of squares = "<<minimum_sum_of_squares<<std::endl;
    //project pts onto plane
    for(int i=1; i<largest_segment_array.m; i++){
        int vertex1, vertex2;
        largest_segment_array.Get(i,vertex1,vertex2);
        VECTOR_3D<T> point1 = triangulated_surface->particles.X(vertex1);
        VECTOR_3D<T> point2 = triangulated_surface->particles.X(vertex2);
        triangulated_surface->particles.X(vertex1) = best_plane.Surface(point1);
        triangulated_surface->particles.X(vertex2) = best_plane.Surface(point2);}

    triangulated_surface->particles.Increase_Array_Size(3);
    int index1 = triangulated_surface->particles.Add_Particle();
    triangulated_surface->particles.X(index1)=best_plane.Surface(VECTOR_3D<T>(-0.5,-0.2,0.5));
    std::cout<<"index 1 = "<<triangulated_surface->particles.X(index1)<<std::endl;
    int index2 = triangulated_surface->particles.Add_Particle();
    triangulated_surface->particles.X(index2)=best_plane.Surface(VECTOR_3D<T>(0.5,-0.2,0.5));
    std::cout<<"index 2 = "<<triangulated_surface->particles.X(index2)<<std::endl;
    int index3 = triangulated_surface->particles.Add_Particle();
    triangulated_surface->particles.X(index3)=best_plane.Surface(VECTOR_3D<T>(0,-0.2,-0.5));
    std::cout<<"index 3 = "<<triangulated_surface->particles.X(index3)<<std::endl;
    int triangle_count = triangulated_surface->triangle_mesh.triangles.m;
    triangulated_surface->triangle_mesh.triangles.Exact_Resize(3,triangle_count+largest_segment_array.m);
    triangle_count++;
    triangulated_surface->triangle_mesh.triangles.Set(triangle_count,index1,index2,index3);
    
    triangulated_surface->particles.Increase_Array_Size(1);
    old_vertex_color_array.Resize(triangulated_surface->triangle_mesh.number_nodes + 1);
    int index = triangulated_surface->particles.Add_Particle();
    triangulated_surface->particles.X(index)=VECTOR_3D<T>(0,-0.161,0);
    int triangle_count = triangulated_surface->triangle_mesh.triangles.m;
    triangulated_surface->triangle_mesh.triangles.Exact_Resize(3,triangle_count+largest_segment_array.m);
    triangle_count++;
    float y = -0.161;
    for(int j=1; j<largest_segment_array.m; j++){
        int vertex1, vertex2;
        largest_segment_array.Get(j,vertex1,vertex2);

        VECTOR_3D<T> point1 = triangulated_surface->particles.X(vertex1);
        triangulated_surface->particles.X(vertex1) = VECTOR_3D<T>(point1.x,y,point1.z);
        VECTOR_3D<T> point2 = triangulated_surface->particles.X(vertex2);
        triangulated_surface->particles.X(vertex2) = VECTOR_3D<T>(point2.x,y,point2.z);

        triangulated_surface->triangle_mesh.triangles.Set(triangle_count,vertex2,vertex1,index);
        triangle_count++;}*/
    //FILE_UTILITIES::Write_To_File<T>("projected_neck.tri",*triangulated_surface);

    int number_of_particles = triangulated_surface->particles.number;
    ARRAY<bool> Particle_Inside_Box = ARRAY<bool>(number_of_particles+1);
    for(int i=1; i<=number_of_particles; i++){
        bool flag = box.Inside(triangulated_surface->particles.X(i),0);
        Particle_Inside_Box.Set(i,flag);}
    std::cout<<"total polygons = "<<triangulated_surface->triangle_mesh.triangles.m<<std::endl;
    for(int j=1; j<=triangulated_surface->triangle_mesh.triangles.m; j++){
        int v1, v2, v3;
        triangulated_surface->triangle_mesh.triangles.Get(j,v1,v2,v3);
        if(!(Particle_Inside_Box(v1)&&Particle_Inside_Box(v2)&&Particle_Inside_Box(v3))){
            Outside_Triangles.Append(j);}
        else{polygon_count++;}}
    Outside_Triangles.Compact();
    triangulated_surface->triangle_mesh.triangles.Remove_Indices(Outside_Triangles.array);
    triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber(condensation_mapping);
    std::cout<<"condensation mapping size = "<<condensation_mapping.m<<std::endl;
    std::cout<<"number of condensed particles = "<<triangulated_surface->particles.number<<std::endl;
    FILE_UTILITIES::Read_From_File<T>(color_filename,old_vertex_color_array);
    std::cout<<"old color array size = "<<old_vertex_color_array.m<<std::endl;
    std::cout<<"old new color array size ="<<new_vertex_color_array.m<<std::endl;
    new_vertex_color_array.Resize(polygon_count);    
    std::cout<<"new color array size ="<<new_vertex_color_array.m<<std::endl;
    //std::cout<<"GOT HERE."<<std::endl;
    int highest_mapping_count = 0;
    for(int i=1; i<=condensation_mapping.m; i++){
        int new_vertex = condensation_mapping(i);
        //std::cout<<i<<" "<<new_vertex;
        if(new_vertex!=0){
            if(new_vertex>highest_mapping_count) highest_mapping_count = new_vertex;
            VECTOR_3D<T> temp = old_vertex_color_array(i);
            new_vertex_color_array(new_vertex) = temp;}}
    std::cout<<"mapping count = "<<highest_mapping_count<<std::endl;
    std::cout<<"polygon count= "<<polygon_count<<std::endl;
    FILE_UTILITIES::Write_To_File<T>(new_tri_filename+".tri",*triangulated_surface);
    std::cout<<"new surface written to file."<<std::endl;
    FILE_UTILITIES::Write_To_File<T>("condensation_mapping",condensation_mapping);
    std::cout<<"condensation mapping written to file."<<std::endl;
    FILE_UTILITIES::Write_To_File<T>(new_tri_filename+".col",new_vertex_color_array);
    std::cout<<"colors written to file"<<std::endl;
    std::string temp;
    std::cin>>temp;
}

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    PARSE_ARGS args;
    args.Set_Extra_Arguments(3,"<tri filename> <new filename> <old color filename>");
    args.Parse(argc,argv);
    std::cout<<"num extra "<<args.Num_Extra_Args()<<std::endl;
    std::string tri_filename = args.Extra_Arg(1);
    std::string new_filename = args.Extra_Arg(2);
    std::string vertex_color_filename = args.Extra_Arg(3);
    MakeCuts<float>(tri_filename,new_filename,vertex_color_filename);
    return 0;
}
//#################################################################
