#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void Convert(const std::string& input_filename,const std::string& output_filename)
{
    typedef VECTOR<T,3> TV;

    std::istream* input=Safe_Open_Input(input_filename,false);
    char buffer[2048];
    ARRAY<TV> vertices;
    ARRAY<VECTOR<int,3> > triangles;
    do{
        input->getline(buffer,2048);
        if(!strlen(buffer)) continue;
        if(buffer[0]=='#') continue;
        else if(buffer[0]=='v'){
            VECTOR<double,3> v;
            sscanf(buffer+2,"%lf %lf %lf",&v.x,&v.y,&v.z);
            vertices.Append((TV)v);}
        else if(buffer[0]=='n' && buffer[1]=='v') continue;
        else if(buffer[0]=='f'){
            VECTOR<int,3> f;
            sscanf(buffer+2,"%d %d %d",&f.x,&f.y,&f.z);
            f-=1;
            triangles.Append(f);}
    }while(!input->eof());
    delete input;

    DEFORMABLE_PARTICLES<TV>& particles=*new DEFORMABLE_PARTICLES<TV>;
    TRIANGULATED_SURFACE<T>* triangulated_surface=TRIANGULATED_SURFACE<T>::Create(particles);

    for(int t=0;t<triangles.m;t++)
        triangulated_surface->mesh.elements.Append(triangles(t));
    triangulated_surface->particles.Preallocate(vertices.m);
    for(int v=0;v<vertices.m;v++)
        triangulated_surface->particles.X(triangulated_surface->particles.Add_Element())=vertices(v);
    triangulated_surface->Update_Number_Nodes();
    Write_To_File<RW>(output_filename,*triangulated_surface);
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false,compute_using_doubles=false;
    std::string input_filename,output_filename;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-compute_using_doubles",&compute_using_doubles,"compute_using_doubles");
    parse_args.Extra(&input_filename,"obj file","obj file to convert");
    parse_args.Extra(&output_filename,"tri file","output tri file name");
    parse_args.Parse();

    if(!File_Extension_Matches_Ignoring_Compression_Suffix(output_filename,"tri",false)){
        std::cerr<<"Not a tri file: "<<output_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(compute_using_doubles){
            std::cout<<"COMPUTING USING DOUBLES!"<<std::endl;
            Convert<double,float>(input_filename,output_filename);
        }else{Convert<float,float>(input_filename,output_filename);}}
    else Convert<double,double>(input_filename,output_filename);
}
