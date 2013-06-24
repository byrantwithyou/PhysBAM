#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>

using namespace PhysBAM;

template<class T>
void Convert(PARSE_ARGS &parse_args)
{
    T scale=1;
    bool opt_fmm=false;
    std::string output_filename,input_filename;
    parse_args.Add("-scale",&scale,"scale","scale");
    parse_args.Add("-o",&output_filename,"file","output filename");
    parse_args.Add("-fmm",&opt_fmm,"use fast marching method");
    parse_args.Extra(&input_filename, "filename", "filename");
    parse_args.Parse();

    ARRAYS<VECTOR<VECTOR<T,3> ,2> > image;
    IMAGE<T>::Read(input_filename,image);
    GRID_2D<T> image_grid(image.m,image.n,0,image.m,0,image.n);image_grid.Set_MAC_Grid();

    if(!output_filename.size()) output_filename=FILE_UTILITIES::Get_Basename(input_filename)+".phi";

    GRID_2D<T> grid(image.m,image.n,-(T).5*scale*(image.m-1),(T).5*scale*(image.m-1),-(T).5*scale*(image.n-1),(T).5*scale*(image.n-1));
    ARRAYS<VECTOR<T,2> > phi(grid);LEVELSET_2D<T> levelset(grid,phi);

    LINEAR_INTERPOLATION<T,VECTOR<T,3> > color_interpolation;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){phi(i,j)=2*image(i,j).Magnitude()/root_three-(T)1;}

    if(opt_fmm) levelset.Fast_Marching_Method(0,10*grid.max_dx_dy);

    FILE_UTILITIES::Write_To_File<T>(output_filename,levelset);
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args(argc,argv);

    bool type_double=false;
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Parse(true);

    if(!type_double) Convert<float>(parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double>(parse_args);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
