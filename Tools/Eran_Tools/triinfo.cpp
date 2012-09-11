#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

bool verbose;
bool verify;

template<class T> void Print_Info(const char *filename)
{
    TRIANGULATED_SURFACE<T> *triangulated_surface=0;
    FILE_UTILITIES::Create_From_File<T>(filename,triangulated_surface);
    triangulated_surface->Update_Bounding_Box();

    TRIANGLE_MESH &mesh = triangulated_surface->triangle_mesh;
    cout << "triangulated_surface.triangle_mesh" << endl;
    cout << "number_nodes = " << mesh.number_nodes << endl;
    cout << "triangles.m = " << mesh.triangles.m << endl;
    cout << endl;

    SOLIDS_PARTICLES<T,VECTOR<T,3> > &particles = triangulated_surface->particles;
    cout << "triangulated_surface.particles" << endl;
    cout << "number = " << particles.number << endl;
    cout << "array_size = " << particles.array_size << endl;
    ARRAY<bool> particle_referenced(particles.array_size);
    for(int i=0;i<mesh.triangles.m;i++) {
        int node1,node2,node3;
        mesh.triangles.Get(i,node1,node2,node3);
        particle_referenced(node1)=particle_referenced(node2)=particle_referenced(node3)=true;
    }
    cout << "particles referenced = " << particle_referenced.Number_True() << endl;
    cout << endl;

    BOX_3D<T> &box = *(triangulated_surface->bounding_box);
    cout << "triangulated_surface.bounding_box" << endl;
    cout << "(xmin,xmax) = (" << box.xmin << ", " << box.xmax << ")" << endl;
    cout << "(ymin,ymax) = (" << box.ymin << ", " << box.ymax << ")" << endl;
    cout << "(zmin,zmax) = (" << box.zmin << ", " << box.zmax << ")" << endl;
    cout << "size = " << box.Size() << endl;
    cout << endl;

    if (verbose)
    {
        cout << endl;
        cout << "VERTICES" << endl;
        SOLIDS_PARTICLES<T,VECTOR<T,3> > &particles = triangulated_surface->particles;
        for (int i = 1; i <= particles.array_size; i++)
            cout << i << ": " << particles.X(i) << endl;
        cout << endl << endl;

        cout << "TRIANGLES" << endl;
        for (int i = 1; i <= mesh.triangles.m; i++)
        {
            cout << i << ":";
            for (int j = 1; j <= 3; j++)
                cout << " " << mesh.triangles(j, i);
            cout << endl;
        }
    }

    if (verify)
    {
        cout << endl;
        cout << "VERIFYING..." << endl;
        SOLIDS_PARTICLES<T,VECTOR<T,3> > &particles = triangulated_surface->particles;
        T min_area=FLT_MAX;
        for (int t = 1; t <= mesh.triangles.m; t++)
        {
            int node1,node2,node3;
            mesh.triangles.Get(t,node1,node2,node3);
            VECTOR<T,3> x1 = particles.X(node1), x2 = particles.X(node2), x3 = particles.X(node3);
            T area = TRIANGLE_3D<T>::Area(x1,x2,x3);
            min_area = PhysBAM::min(min_area, area);
            if (area < 1e-6) {
                std::cout << "Triangle " << t << ": [nodes: " << node1 << ", " << node2 << ", " << node3 << "] area=" << area << std::endl;
            }
        }
        std::cout << "Minimum area: " << min_area << std::endl;
    }
}

int main(int argc, char *argv[])
{
    bool type_double=false;   // float by default
    char filename[256];

    PARSE_ARGS parse_args;

    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-v",&verbose, "display verbose information");
    parse_args.Add("-verify",&verify, "verify triangles");
    parse_args.Set_Extra_Arguments(1, "<filename>");

    int extraarg = parse_args.Parse();

    if (extraarg < argc)
        strcpy(filename, argv[extraarg]);
    else
        return -1;

    cout << "Filename: " << filename << " [" << ((type_double)?"double":"float") << "]" << endl;

    if(!type_double) Print_Info<float>(filename);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Do_It<double>(filename);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
