#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <fstream>
#include "../../Personal_Libraries/Eran_Library/POLYGONAL_SURFACE.h"
#include "RIGID_BODY_BUILDER.h"

using namespace PhysBAM;
using namespace std;
using namespace FILE_UTILITIES;

template<class T> void Build_It(const char *input_filename, PARSE_ARGS &parse_args)
{
    RIGID_BODY_BUILDER<T> builder;

    if (parse_args.Is_Value_Set("-m"))
        builder.Set_Mass((T)parse_args.Get_Double_Value("-m"));
    else if (parse_args.Is_Value_Set("-d"))
        builder.Set_Density((T)parse_args.Get_Double_Value("-d"));
    builder.Set_Surface_Roughness((T)parse_args.Get_Double_Value("-r"));

    if (parse_args.Is_Value_Set("-o"))
        cout  << "Weird" << endl;

    bool thin_shell = parse_args.Get_Option_Value("-thin_shell");

    if (!Is_Ply_File(input_filename) && !Is_Tri_File(input_filename))
    {
        cerr << "Error: " << input_filename << " is not a .tri or .ply file" << endl;
        exit(-1);
    }

    ifstream input(input_filename);
    if (!input)
    {
        cerr << "Error reading file " << input_filename << endl;
        exit(-1);
    }
    
    TRIANGULATED_SURFACE<T> *triangulated_surface = 0;

    if (Is_Ply_File(input_filename))
    {
        cout << "Reading in polygonal surface" << endl;
        POLYGONAL_SURFACE *polygonal_surface = new POLYGONAL_SURFACE();
        polygonal_surface->Read(input);
    //    polygonal_surface->Print();

        triangulated_surface = RIGID_BODY_BUILDER<T>::Create_Triangulated_Surface(polygonal_surface);
    }
    else if (Is_Tri_File(input_filename))
    {
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_filename);
        triangulated_surface=TRIANGULATED_SURFACE<T>::Create();
        triangulated_surface->template Read<T>(*input);delete input;
    }

    if (parse_args.Is_Value_Set("-i")) 
    {
        std::string rgd_input_filename;
        rgd_input_filename=parse_args.Get_String_Value("-i");

        RIGID_BODY_3D<T> rigid_body;
        FILE_UTILITIES::Read_From_File<T>(rgd_input_filename,rigid_body);

        cout << "** Using input rgd file " << rgd_input_filename << endl;
        cout << "** Got mass " << rigid_body.mass << endl;
        cout << "** Got position " << rigid_body.position << endl;
        cout << "** Got orientation " << rigid_body.orientation << endl;

        builder.Set_Mass(rigid_body.mass);

        for (int i = 1; i <= triangulated_surface->particles.array_size; i++)
        {
            triangulated_surface->particles.X(i) = rigid_body.position + 
                    rigid_body.orientation.Rotate(triangulated_surface->particles.X(i));
        }
    }

    triangulated_surface->Update_Bounding_Box();

    RIGID_BODY_3D<T> *rigid_body = builder.Build_Rigid_Body(triangulated_surface, thin_shell);

    if (parse_args.Get_Option_Value("-no"))
    {
        cout << "RIGID BODY parameters:" << endl;
        cout << "mass = " << rigid_body->mass << endl;
        cout << "inertia_tensor = " << rigid_body->inertia_tensor << endl;
        cout << "surface_roughness = " << rigid_body->surface_roughness << endl;
        cout << "position = " << rigid_body->position << endl;
        cout << "orientation = " << rigid_body->orientation << endl;
    }
    else
    {
        string basename;
        char tri_filename[256], rgd_filename[256];

        if (parse_args.Is_Value_Set("-o"))
        {
            basename=Get_Basename(parse_args.Get_String_Value("-o"));
        }
        else
        {
            basename=Get_Basename(input_filename);
        }
        sprintf(tri_filename, "%s.tri", basename.c_str());
        sprintf(rgd_filename, "%s.rgd", basename.c_str());

//        if (Is_Ply_File(input_filename))
        {
            cout << "Writing triangulated surface to " << tri_filename << endl;
            ofstream tri_output(tri_filename, ios::binary);
            rigid_body->triangulated_surface->template Write<T>(tri_output);
        }

        cout << "Writing rigid body to " << rgd_filename << endl;
        ofstream rgd_output(rgd_filename, ios::binary);
        rigid_body->template Write<T>(rgd_output);
    }
}

int main(int argc, const char *argv[])
{
    PARSE_ARGS parse_args;
    bool type_double = false;
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-no", "print results, don't output to file");
    parse_args.Add_Option_Argument("-thin_shell");
    parse_args.Add_String_Argument("-o", "", "output filename");
    parse_args.Add_Double_Argument("-m", 1, "mass");
    parse_args.Add_Double_Argument("-d", 1, "density");
    parse_args.Add_Double_Argument("-r", 1e-8, "surface_roughness");
    parse_args.Add_String_Argument("-i", "", "use this rgd file's mass, position, orientation");
    parse_args.Set_Extra_Arguments(1, "<ply or tri file>");

    int extraarg = parse_args.Parse(argc, argv);

    char input_filename[256];
    if (extraarg < argc)
        strcpy(input_filename, argv[extraarg]);
    else
        return -1;

    if (parse_args.Is_Value_Set("-float")) type_double = false;
    if (parse_args.Is_Value_Set("-double")) type_double = true;

    if (type_double) Build_It<double>(input_filename, parse_args); 
    else Build_It<float>(input_filename, parse_args);

}
