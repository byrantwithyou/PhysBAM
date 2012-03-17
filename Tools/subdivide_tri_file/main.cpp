#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include <string.h> // for gcc 2.96 for str* functions

using namespace std;
using namespace PhysBAM;

void get_basename(const char *name, char *basename)
{
    if (!strcmp(&name[strlen(name)-4],".tri")){
        strncpy(basename,name,strlen(name)-4);
        basename[strlen(name)-4]=0;}
    else strcpy(basename, name);
}

void print_usage_and_exit(const char *progname)
{
    cerr << "Usage: " << progname << " [-l <num_levels>] <filename>" << endl;
    exit(-1);
}

int main(int argc,char *argv[])
{
    char input_name[1024]="../Public_Library/Data/Rigid_Bodies/Bones/Cranium.tri";
    int levels=1;

    // Parse arguments
    if(argc<=1) cout << "Using default values" << endl;
    else if(argc==2) strcpy(input_name, argv[1]);
    else if (argc==4){
            if(strcmp(argv[1],"-l")) print_usage_and_exit(argv[0]);
            levels=atoi(argv[2]);
        if(levels<1 || levels>20) print_usage_and_exit(argv[0]);
        strcpy(input_name,argv[3]);}
    else print_usage_and_exit(argv[0]);

    char basename[1024],input_filename[1024],output_filename[1024];
    get_basename(input_name,basename);
    sprintf(input_filename,"%s.tri", basename);
    sprintf(output_filename,"%s_%d.tri", basename, levels);

    cout << "Input filename: " << input_filename << endl;
    cout << "Output filename: " << output_filename << endl;
    cout << "Levels: " << levels << endl;
    
    TRIANGULATED_SURFACE<float> &surface=*TRIANGULATED_SURFACE<float>::Create();
        FILE_UTILITIES::Read_From_File<float>(input_filename,surface);

    for(int i=0;i<levels;i++) surface.Loop_Subdivide();

        FILE_UTILITIES::Write_To_File<float>(output_filename,surface);
}
