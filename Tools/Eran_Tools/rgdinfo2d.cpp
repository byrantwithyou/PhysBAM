#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T> void Do_It(const char *filename, PARSE_ARGS &parse_args)
{
    RIGID_BODY_2D<T> rigid_body;

    if (filename && File_Exists(filename))
    {
        FILE_UTILITIES::template Read_From_File<T>(filename,rigid_body);

        cout << "mass = " << rigid_body.mass << endl;
        cout << "moment_of_inertia = " << rigid_body.moment_of_inertia << endl;
        cout << "surface_roughness = " << rigid_body.surface_roughness << endl;
        cout << "position = " << rigid_body.position << endl;
        cout << "orientation = " << rigid_body.orientation << endl;
        cout << "velocity = " << rigid_body.velocity << endl;
        cout << "angular_momentum = " << rigid_body.angular_momentum << endl;
    }
    else
    {
        if (parse_args.Get_Option_Value("-w"))
            cout << "Creating new file " << filename << endl;
        else
        {
            cerr << "File doesn't exist: " << filename << endl;
            cerr << "Use -w to write new file" << endl;
            exit(-1);
        }
    }

    if (parse_args.Get_Option_Value("-w"))
    {
        cout << endl;
        cout << "New values:" << endl;
        if (parse_args.Is_Value_Set("-m"))
        {
            rigid_body.mass = (T)parse_args.Get_Double_Value("-m");
            cout << "mass = " << rigid_body.mass << endl;
        }
        if (parse_args.Is_Value_Set("-i"))
        {
            rigid_body.moment_of_inertia = parse_args.Get_Double_Value("-i");
            cout << "moment_of_inertia = " << rigid_body.moment_of_inertia << endl;
        }
        if (parse_args.Is_Value_Set("-p"))
        {
            rigid_body.position = parse_args.Get_Vector_2D_Value("-p");
            cout << "position = " << rigid_body.position << endl;
        }
        cout << "WRITING TO FILE..." << endl;

        FILE_UTILITIES::template Write_To_File<T>(filename,rigid_body);
    }
}

int main(int argc, char *argv[])
{
    bool type_double = false;
    char filename[256];

    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-w", "write new values to file");
    parse_args.Add_Double_Argument("-m", 1, "mass");
    parse_args.Add_Double_Argument("-i", 1, "moment of inertia");
    parse_args.Add_Vector_2D_Argument("-p", VECTOR<double,2>(0,0), "position");
    parse_args.Set_Extra_Arguments(1, "<filename>");

    int extraarg = parse_args.Parse(argc, argv);

    if (extraarg < argc)
        strcpy(filename, argv[extraarg]);
    else
        return -1;

    if (!Is_Rgd2D_File(filename))
    {
        cerr << "Not a rgd2d file: " << filename << endl;
        return -1;
    }

    if (parse_args.Get_Option_Value("-double")) type_double = true;
    if (parse_args.Get_Option_Value("-float")) type_double = false;

    cout << "Filename: " << filename << " [" << ((type_double)?"double":"float") << "]" << endl;

    if(!type_double) Do_It<float>(filename, parse_args);
    else{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Do_It<double>(filename, parse_args);
#else
        std::cerr << "Compiled without double support" << std::endl;
#endif
    }
}
