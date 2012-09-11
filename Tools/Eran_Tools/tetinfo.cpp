#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

bool verbose, verify;

template<class T> void Print_Info(const char *filename)
{
    TETRAHEDRON_MESH tet_mesh;
    SOLIDS_PARTICLES<T,VECTOR<T,3> > particles;
    TETRAHEDRALIZED_VOLUME<T> tet_volume(tet_mesh,particles);
    EMBEDDED_TETRAHEDRALIZED_VOLUME<T> etv(tet_volume);

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    if(!input) exit(1);
    if (Get_File_Extension(filename) == "etv") etv.template Read<T>(*input);
    else tet_volume.template Read<T>(*input);
    delete input;

    tet_volume.Update_Bounding_Box();

    cout << "tetrahedralized_volume.tetrahedron_mesh" << endl;
    cout << "number_nodes = " << tet_mesh.number_nodes << endl;
    cout << "tetrahedrons.m = " << tet_mesh.tetrahedrons.m << endl;
    cout << endl;

    tet_volume.Initialize_Triangulated_Surface();
    TRIANGULATED_SURFACE<T> &tri_surface=*tet_volume.triangulated_surface;
    TRIANGLE_MESH &triangle_mesh=tri_surface.triangle_mesh;
    cout << "tetrahedralized_volume.triangulated_surface.triangle_mesh" << endl;
    cout << "number_nodes = " << triangle_mesh.number_nodes << endl;
    cout << "triangles.m = " << triangle_mesh.triangles.m << endl;
    cout << endl;

    BOX_3D<T> &box = *tet_volume.bounding_box;
    cout << "tetrahedralized_volume.bounding_box" << endl;
    cout << "(xmin,xmax) = (" << box.xmin << ", " << box.xmax << ")" << endl;
    cout << "(ymin,ymax) = (" << box.ymin << ", " << box.ymax << ")" << endl;
    cout << "(zmin,zmax) = (" << box.zmin << ", " << box.zmax << ")" << endl;
    cout << endl;

    if (verbose)
    {
        int i;

        cout << endl;
        cout << "VERTICES" << endl;
        for (i = 1; i <= particles.array_size; i++)
            cout << i << ": " << particles.X(i) << endl;
        cout << endl << endl;

        cout << "TETRAHEDRA" << endl;
        for (i = 1; i <= tet_mesh.tetrahedrons.m; i++)
        {
            cout << i << ":";
            for (int j = 1; j <= 4; j++)
                cout << " " << tet_mesh.tetrahedrons(j, i);
            cout << endl;
        }

        cout << "TRIANGLES" << endl;
        for (i = 1; i <= triangle_mesh.triangles.m; i++)
        {
            cout << i << ":";
            for (int j = 1; j <= 3; j++)
                cout << " " << triangle_mesh.triangles(j, i);
            cout << endl;
        }
    }

    if (verify)
    {
        ARRAYS<VECTOR<bool,1> > node_is_referenced(1,tet_mesh.number_nodes);
        for(int i=0;i<tet_mesh.tetrahedrons.m;i++){
            int node1,node2,node3,node4;tet_mesh.tetrahedrons.Get(i,node1,node2,node3,node4);
            node_is_referenced(node1)=node_is_referenced(node2)=node_is_referenced(node3)=node_is_referenced(node4)=true;
        }

        for (int i=1;i<=particles.array_size;i++){
            if (i <= tet_mesh.number_nodes && !node_is_referenced(i))
                std::cout << "Node " << i << " is active but not referenced" << std::endl;
            if (i > tet_mesh.number_nodes)
                std::cout << "Node " << i << " is active but out of bounds" << std::endl;
        }
    }
}

int main(int argc, char *argv[])
{
    bool type_double = false;   // float by default
    char filename[256];

    PARSE_ARGS parse_args;

    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add_Option_Argument("-v", "display verbose information");
    parse_args.Add_Option_Argument("-verify", "verify triangles");
    parse_args.Set_Extra_Arguments(1, "<filename>");

    int extraarg = parse_args.Parse();

    if (extraarg < argc)
        strcpy(filename, argv[extraarg]);
    else
        return -1;

    verbose = parse_args.Get_Option_Value("-v");
    verify = parse_args.Get_Option_Value("-verify");

    cout << "Filename: " << filename << " [" << ((type_double)?"double":"float") << "]" << endl;

    if(!type_double) Print_Info<float>(filename);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Do_It<double>(filename);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
