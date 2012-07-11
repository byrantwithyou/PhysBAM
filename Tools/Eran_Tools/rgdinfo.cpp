#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T> void Do_It(const char *filename,PARSE_ARGS &parse_args)
{
    std::auto_ptr<std::istream> input=FILE_UTILITIES::Safe_Open_Input(filename);

    for(int i=1;;i++){
        RIGID_BODY_3D<T> rigid_body(*new RIGID_BODY_PARTICLES<VECTOR<T,3> >);
        rigid_body.template Read<T>(*input);
        if(!*input) break;

        cout<<"Rigid body "<<i<<std::endl;
        cout<<"mass = "<<rigid_body.Mass().mass<<endl;
        cout<<"inertia_tensor = "<<rigid_body.Mass().inertia_tensor<<endl;
        cout<<"surface_roughness = "<<rigid_body.surface_roughness<<endl;
        cout<<"position = "<<rigid_body.Frame().t<<endl;
        cout<<"orientation = "<<rigid_body.Frame().r<<endl;
        cout<<"velocity = "<<rigid_body.Twist().linear<<endl;
        cout<<"angular_momentum = "<<rigid_body.Twist().angular<<endl;

        if (parse_args.Get_Option_Value("-w")){
            cout<<endl;
            cout<<"New values:"<<endl;
            if (parse_args.Is_Value_Set("-m")){
                rigid_body.Mass().mass=(T)parse_args.Get_Double_Value("-m");
                cout<<"mass = "<<rigid_body.Mass().mass<<endl;}
            if (parse_args.Is_Value_Set("-i")){
                VECTOR<T,3> moments_of_inertia(parse_args.Get_Vector_3D_Value("-i"));
                rigid_body.Mass().inertia_tensor=DIAGONAL_MATRIX<T,3>(moments_of_inertia.x,moments_of_inertia.y,moments_of_inertia.z);
                cout<<"inertia_tensor = "<<rigid_body.Mass().inertia_tensor<<endl;}
            if (parse_args.Is_Value_Set("-p")){
                rigid_body.Frame().t=(VECTOR<T,3>)parse_args.Get_Vector_3D_Value("-p");
                cout<<"position = "<<rigid_body.Frame().t<<endl;}
            cout<<"WRITING TO FILE..."<<endl;
            FILE_UTILITIES::template Write_To_File<T>(filename,rigid_body);}}
}

int main(int argc,char *argv[])
{
    bool type_double=false;

    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-w","write new values to file");
    parse_args.Add_Double_Argument("-m",1,"mass");
    parse_args.Add_Vector_3D_Argument("-i",VECTOR<double,3>(1,1,1),"moments of inertia");
    parse_args.Add_Vector_3D_Argument("-p",VECTOR<double,3>(),"position");
    parse_args.Set_Extra_Arguments(1,"<filename>");

    parse_args.Parse();

    std::string filename=parse_args.Extra_Arg(0);

#if 0
    if (!Is_Rgd_File(filename)){
        cerr<<"Not a rgd file: "<<filename<<endl;
        return EXIT_FAILURE;}
#endif
    if (!File_Exists(filename)){
        cerr<<"Input file "<<filename<<" does not exist!"<<endl;
        return EXIT_FAILURE;}

    if (parse_args.Get_Option_Value("-double")) type_double=true;
    if (parse_args.Get_Option_Value("-float")) type_double=false;

    cout<<"Filename: "<<filename<<" ["<<((type_double)?"double":"float")<<"]"<<endl;

    if(!type_double) Do_It<float>(filename.c_str(),parse_args);
    else{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Do_It<double>(filename.c_str(),parse_args);
#else
        std::cerr<<"Compiled without double support"<<std::endl;
#endif
    }
}
