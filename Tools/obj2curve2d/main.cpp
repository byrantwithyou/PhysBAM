#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void Convert(const std::string& input_filename,const std::string& output_filename)
{
    typedef VECTOR<T,2> TV;

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_filename,false);
    char buffer[2048];
    ARRAY<TV> vertices;
    ARRAY<VECTOR<int,2> > segments;
    do{
        input->getline(buffer,2048);
        if(!strlen(buffer)) continue;
        if(buffer[0]=='#') continue;
        else if(buffer[0]=='v'){
            VECTOR<double,2> v;
            sscanf(buffer+2,"%lf %lf",&v.x,&v.y);
            vertices.Append((TV)v);}
        else if(buffer[0]=='n' && buffer[1]=='v') continue;
        else if(buffer[0]=='f'){
            VECTOR<int,2> f;
            sscanf(buffer+2,"%d %d",&f.x,&f.y);
            segments.Append(f);}
    }while(!input->eof());
    delete input;

    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    SEGMENTED_CURVE_2D<T>* segmented_curve_2d=SEGMENTED_CURVE_2D<T>::Create(particles);

    for(int t=0;t<segments.m;t++)
        segmented_curve_2d->mesh.elements.Append(segments(t));
    segmented_curve_2d->particles.array_collection->Preallocate(vertices.m);
    for(int v=0;v<vertices.m;v++)
        segmented_curve_2d->particles.X(segmented_curve_2d->particles.array_collection->Add_Element())=vertices(v);
    segmented_curve_2d->Update_Number_Nodes();
    FILE_UTILITIES::Write_To_File<RW>(output_filename,*segmented_curve_2d);
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    Initialize_Read_Write_General_Structures();

    bool type_double=false;

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float","data should be in float format");
    parse_args.Add_Option_Argument("-double","data should be in double format");
    parse_args.Add_Option_Argument("-compute_using_doubles");
    parse_args.Set_Extra_Arguments(1,"<obj file>","<obj file> obj file to convert");
    parse_args.Set_Extra_Arguments(2,"<curve2d file>","<curve2d file> output curve2d file name");

    parse_args.Parse(argc, argv);

    std::string input_filename=parse_args.Extra_Arg(1);
    std::string output_filename=parse_args.Extra_Arg(2);

    if(parse_args.Get_Option_Value("-float")) type_double=false;
    if(parse_args.Get_Option_Value("-double")) type_double=true;

    if(!FILE_UTILITIES::Is_Curve2D_File(output_filename)){
        std::cerr<<"Not a curve2d file: "<<output_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(parse_args.Get_Option_Value("-compute_using_doubles")){
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
            std::cout<<"COMPUTING USING DOUBLES!"<<std::endl;
            Convert<double,float>(input_filename,output_filename);
#else
            std::cerr<<"Double support not enabled."<<std::endl;exit(1);
#endif
        }else{Convert<float,float>(input_filename,output_filename);}}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double,double>(input_filename,output_filename);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
