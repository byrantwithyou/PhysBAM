#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
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
    std::string output_filename=parse_args.Get_String_Value("-o");

    TRIANGULATED_SURFACE<T>* merged_surface=TRIANGULATED_SURFACE<T>::Create();

    FRAME_3D<T> main_frame;

    for(int i=0;i<parse_args.Num_Extra_Args();i++){
        std::cout << "Reading " << parse_args.Extra_Arg(i) << std::endl;
        std::string filename=parse_args.Extra_Arg(i);
        TRIANGULATED_SURFACE<T> *surface=0;FILE_UTILITIES::Create_From_File<T>(filename,surface);
        // TODO: MAKE OPTIONAL!!
        std::string rgd_filename=FILE_UTILITIES::Get_Basename(filename)+".rgd";
        FRAME_3D<T> xform;
        if(FILE_UTILITIES::File_Exists(rgd_filename)){
            std::cout << "Reading rgd " << rgd_filename << std::endl;
            RIGID_BODY_3D<T> rigid_body;FILE_UTILITIES::Read_From_File<T>(rgd_filename,rigid_body);
            if(i==1) main_frame=rigid_body.frame;
            else xform=main_frame.Inverse()*rigid_body.frame;}
        surface->Discard_Valence_Zero_Particles_And_Renumber();
        int old_number_of_particles=merged_surface->particles.number,old_number_of_triangles=merged_surface->triangle_mesh.triangles.m;
        merged_surface->particles.Add_Particles(surface->particles.number);
        for(int p=0;p<surface->particles.number;p++){
            surface->particles.X(p)=xform*surface->particles.X(p);
            merged_surface->particles.Copy_Particle(surface->particles,p,old_number_of_particles+p);}
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
    parse_args.Set_Extra_Arguments(-1,"<tri file> <tri file> ...");
    parse_args.Parse(argc,argv);
    if(parse_args.Num_Extra_Args()<1) return 1;
    Do_It<float>();
}
