#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <fstream>
#include <sstream>
#include "../../Personal_Libraries/Eran_Library/PARSE_UTILS.h"
#include <string.h>

const std::string DEFAULT_SIMMDIR = "../../Projects/articulated_rigid_bodies/UpperExtremityModel/bones";
const std::string DEFAULT_BONEDIR = "../../Public_Data/Rigid_Bodies/Converted_SIMM_Bones/";

using namespace PhysBAM;
using namespace std;

int main(int argc,char** argv)
{
    std::string bone_directory,simm_directory,filename,physbam_bone_name;
    bool binary=false;
    
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-p","","PhysBAM bone directory","PhysBAM bone directory");
    parse_args.Add_String_Argument("-s","","SIMM bone directory","SIMM bone directory");
    parse_args.Add_String_Argument("-b","","PhysBAM bone name","PhysBAM bone name");
    parse_args.Add_Option_Argument("-g","Binary file");
    parse_args.Set_Extra_Arguments(1, "<muscle_file>");
    parse_args.Parse(argc,argv);

    bone_directory=DEFAULT_BONEDIR;
    simm_directory=DEFAULT_SIMMDIR;

    if(parse_args.Is_Value_Set("-p")) bone_directory=parse_args.Get_String_Value("-p");
    if(parse_args.Is_Value_Set("-s")) simm_directory=parse_args.Get_String_Value("-s");
    if(parse_args.Is_Value_Set("-b")) physbam_bone_name=parse_args.Get_String_Value("-b");
    binary=parse_args.Get_Option_Value("-g");

    if (parse_args.Num_Extra_Args() != 1) parse_args.Print_Usage(true);
    else filename = parse_args.Extra_Arg(1);
    physbam_bone_name=filename;
    if(parse_args.Is_Value_Set("-b")) physbam_bone_name=parse_args.Get_String_Value("-b");

    std::cout<<"Making triangulated surface from SIMM .asc bone file\n";
    // Make triangulated surface from simm bone
    std::string fullfile=simm_directory+"/"+filename+".asc";std::cout<<"opening: "<<fullfile<<std::endl;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(fullfile,binary);
    std::string norm;int num_vertices,num_polygons;float minX,minY,minZ,maxX,maxY,maxZ;
    *input>>norm>>num_vertices>>num_polygons>>minX>>minY>>minZ>>maxX>>maxY>>maxZ; // for bounding box
    std::cout<<"verts: "<<num_vertices<<", polygons"<<num_polygons<<", bounding box: "<<minX<<":"<<maxX<<","<<minY<<":"<<maxY<<","<<minZ<<":"<<maxZ<<std::endl;
    ARRAY<VECTOR_3D<float> > vertices(num_vertices);ARRAYS<VECTOR_3D<float> > vertex_normals(1,num_vertices);float x,y,z;
    for(int i=0;i<num_vertices;i++){*input>>vertices(i).x>>vertices(i).y>>vertices(i).z>>x>>y>>z;vertex_normals(1,i)=VECTOR_3D<float>(x,y,z);}
    int num_sides,s1,s2,s3;ARRAYS<int> triangles(3,num_polygons);
    for(int p=0;p<num_polygons;p++){
        *input>>num_sides;
        if(num_sides==3){*input>>s1>>s2>>s3;triangles.Set(p,s1+1,s2+1,s3+1);}
        else{
            *input>>s1>>s2>>s3;triangles.Set(p,s1+1,s2+1,s3+1);
            for(int s=0;s<num_sides-3;s++){s2=s3;*input>>s3;triangles.Append(s2+1,s3+1,s1+1);}}}
    std::cout<<"num_polygons: "<<triangles.m<<std::endl;
    TRIANGULATED_SURFACE<float>* triangulated_surface=new TRIANGULATED_SURFACE<float>(*new TRIANGLE_MESH(triangles),*new SOLIDS_PARTICLES<float,VECTOR_3D<float> >);
    triangulated_surface->particles.Add_Particles(num_vertices);for(int i=0;i<num_vertices;i++)triangulated_surface->particles.X(i)=vertices(i);
    triangulated_surface->vertex_normals=&vertex_normals;

    //write out triangulated surface
    FILE_UTILITIES::Write_To_File<float>(bone_directory+"/"+physbam_bone_name+".tri",*triangulated_surface);

    return 0;
}
