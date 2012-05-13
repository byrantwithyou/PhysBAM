#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace PhysBAM;

PARSE_ARGS parse_args;

template<class T>
void Do_It()
{
    int num_surfaces;
    std::string surface_name;
    FRAME_3D<T> input_frame;
    T scaling_factor;
    FRAME_3D<T> subtract_frame,align_frame;

    std::string output_filename=parse_args.Get_String_Value("-o");

    TRIANGULATED_SURFACE<T>* merged_surface=TRIANGULATED_SURFACE<T>::Create();

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(parse_args.Extra_Arg(0),false);
    *input>>num_surfaces;
    for(int i=0;i<num_surfaces;i++){
        *input>>surface_name>>scaling_factor>>input_frame;
        std::cout << "Reading " << surface_name << std::endl;
        TRIANGULATED_SURFACE<T> *surface=0;FILE_UTILITIES::Create_From_File<T>(surface_name+".tri.gz",surface);
        surface->Rescale(scaling_factor);
        surface->triangle_mesh.number_nodes=surface->particles.number;
        if(i==1){
            subtract_frame=input_frame;
            RIGID_BODY_3D<T> rigid_body;FILE_UTILITIES::Read_From_File<T>(surface_name+".rgd.gz",rigid_body);
            align_frame=rigid_body.frame;}
        for(int p=0;p<surface->particles.number;p++) surface->particles.X(p)=align_frame*subtract_frame.Inverse()*input_frame*surface->particles.X(p);
        surface->Discard_Valence_Zero_Particles_And_Renumber();
        int old_number_of_particles=merged_surface->particles.number,old_number_of_triangles=merged_surface->triangle_mesh.triangles.m;
        merged_surface->particles.Add_Particles(surface->particles.number);
        for(int p=0;p<surface->particles.number;p++) merged_surface->particles.Copy_Particle(surface->particles,p,old_number_of_particles+p);
        merged_surface->triangle_mesh.triangles.Resize(3,merged_surface->triangle_mesh.triangles.m+surface->triangle_mesh.triangles.m);
        for(int t=0;t<surface->triangle_mesh.triangles.m;t++){int node1,node2,node3;surface->triangle_mesh.triangles.Get(t,node1,node2,node3);
            merged_surface->triangle_mesh.triangles.Set(old_number_of_triangles+t,node1+old_number_of_particles,node2+old_number_of_particles,node3+old_number_of_particles);}}

    merged_surface->triangle_mesh.number_nodes=merged_surface->particles.number;
    std::cout << "Writing " << output_filename << std::endl;
    FILE_UTILITIES::Write_To_File<T>(output_filename,*merged_surface);
}

int main(int argc,char *argv[])
{
    parse_args.Add_String_Argument("-o","merged.tri");
    parse_args.Set_Extra_Arguments(-1,"<tri transform file>");
    parse_args.Parse(argc,argv);
    if(parse_args.Num_Extra_Args()<1) return 1;
    Do_It<float>();
}
