#include <Tools/Log/LOG.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void Convert(const std::string& input_filename,const std::string& output_filename)
{
    typedef VECTOR<T,3> TV;

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_filename,false);
    char buffer[2048];
    ARRAY<TV> vertices;
    ARRAY<VECTOR<int,4> > tetrahedra;
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
            VECTOR<int,4> f;
            sscanf(buffer+2,"%d %d %d %d",&f[0],&f[1],&f[2],&f[3]);
            tetrahedra.Append(f);}
    }while(!input->eof());
    delete input;

    DEFORMABLE_PARTICLES<TV>& particles=*new DEFORMABLE_PARTICLES<TV>;
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);

    for(int t=0;t<tetrahedra.m;t++)
        tetrahedralized_volume->mesh.elements.Append(tetrahedra(t));
    tetrahedralized_volume->particles.Preallocate(vertices.m);
    for(int v=0;v<vertices.m;v++)
        tetrahedralized_volume->particles.X(tetrahedralized_volume->particles.Add_Element())=vertices(v);
    tetrahedralized_volume->Update_Number_Nodes();
    FILE_UTILITIES::Write_To_File<RW>(output_filename,*tetrahedralized_volume);
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
    parse_args.Extra(&output_filename,"tet file","output tet file name");
    parse_args.Parse();

    if(!FILE_UTILITIES::Is_Tet_File(output_filename)){
        std::cerr<<"Not a tet file: "<<output_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(compute_using_doubles){
            std::cout<<"COMPUTING USING DOUBLES!"<<std::endl;
            Convert<double,float>(input_filename,output_filename);
        }else{Convert<float,float>(input_filename,output_filename);}}
    else Convert<double,double>(input_filename,output_filename);
}
