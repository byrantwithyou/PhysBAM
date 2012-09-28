#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_SURFACE.h"
#include <stdio.h>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T> void Phi_To_Tri(const char *input_filename,
                                  const char *output_filename,
                                  T resample_factor=1)
{
    LEVELSET_IMPLICIT_SURFACE<T> *implicit_surface=0;
    FILE_UTILITIES::Create_From_File<T>(input_filename,implicit_surface);

    TRIANGULATED_SURFACE<T> *triangulated_surface=TRIANGULATED_SURFACE<T>::Create();

    if(resample_factor!=1)
    {
        GRID_3D<T> grid=implicit_surface->levelset.grid;
        T MAC_offset=grid.MAC_offset;
        grid.Initialize((int)(resample_factor*grid.m),(int)(resample_factor*grid.n),(int)(resample_factor*grid.mn),grid.xmin,grid.xmax,grid.ymin,grid.ymax,grid.zmin,grid.zmax);
        if(MAC_offset==0.5) grid.Set_MAC_Grid();
        implicit_surface->levelset.Calculate_Triangulated_Surface_From_Marching_Tetrahedra(grid,*triangulated_surface);
    }
    else
    {
        implicit_surface->levelset.Calculate_Triangulated_Surface_From_Marching_Tetrahedra(*triangulated_surface);
    }

    ofstream out(output_filename,ios::binary);
    triangulated_surface->template Write<T>(out);

    delete implicit_surface;
    delete triangulated_surface;
}

int main(int argc,char *argv[])
{
    bool type_double=false;
    char output_filename[256];
    std::string filename;

    double resample=1;
    std::string output_filename;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-resample",&resample,"value","resample grid");
    parse_args.Add("-o",&output_filename,"file","output","output file name");
    parse_args.Extra(&filename,"filename","filename");
    parse_args.Parse();

    if(!Is_Phi_File(filename))
    {
        cerr<<"Not a phi file: "<<filename<<endl;
        return -1;
    }

    if(!output_filename.size())
        output_filename=FILE_UTILITIES::Get_Basename(filename)+".tri";

    cout<<"Input filename: "<<filename<<endl;
    cout<<"Output filename: "<<output_filename<<endl;

    if(!type_double) Phi_To_Tri<float>(filename,output_filename,resample);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Phi_To_Tri<double>(filename,output_filename,resample);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
