#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T> void Print_Info(const char *filename, bool verbose)
{
    SEGMENTED_CURVE_2D<T> *segmented_curve=0;
    FILE_UTILITIES::Create_From_File<T>(filename,segmented_curve);
    segmented_curve->Update_Bounding_Box();

    SEGMENT_MESH &mesh = segmented_curve->segment_mesh;
    cout << "segmented_curve.segment_mesh" << endl;
    cout << "number_nodes = " << mesh.number_nodes << endl;
    cout << "segments.m = " << mesh.segments.m << endl;
    cout << endl;

    BOX_2D<T> &box = *(segmented_curve->bounding_box);
    cout << "segmented_curve.bounding_box" << endl;
    cout << "(xmin,xmax) = (" << box.xmin << ", " << box.xmax << ")" << endl;
    cout << "(ymin,ymax) = (" << box.ymin << ", " << box.ymax << ")" << endl;
    cout << endl;

    if (verbose)
    {
        int i;

        cout << endl;
        cout << "VERTICES" << endl;
        SOLIDS_PARTICLES<T,VECTOR<T,2> > &particles = segmented_curve->particles;
        for (i = 1; i <= particles.number; i++)
            cout << i << ": " << particles.X(i) << endl;
        cout << endl << endl;

        cout << "SEGMENTS" << endl;
        for (i = 1; i <= mesh.segments.m; i++)
        {
            cout << i << ":";
            for (int j = 1; j <= 2; j++)
                cout << " " << mesh.segments(j, i);
            cout << endl;
        }
    }
}

int main(int argc, char *argv[])
{
    bool type_double = false;   // float by default
    bool verbose = false;
    char filename[256];

    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-v", "display verbose information");
    parse_args.Set_Extra_Arguments(1, "<filename>");

    int extraarg = parse_args.Parse();

    if (extraarg < argc)
        strcpy(filename, argv[extraarg]);
    else
        return -1;

    if (!Is_Curve2D_File(filename))
    {
        cerr << "Not a curve2d file: " << filename << endl;
        return -1;
    }

    if (parse_args.Get_Option_Value("-double")) type_double = true;
    if (parse_args.Get_Option_Value("-float")) type_double = false;

    verbose = parse_args.Get_Option_Value("-v");

    cout << "Filename: " << filename << endl;

    if(!type_double) Print_Info<float>(filename, verbose);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Print_Info<double>(filename, verbose);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
