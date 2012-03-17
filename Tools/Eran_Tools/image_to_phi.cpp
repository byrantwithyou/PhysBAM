#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T>
void Convert(PARSE_ARGS &parse_args)
{
    std::string input_filename=parse_args.Extra_Arg(1);
    std::string output_filename;

    ARRAYS<VECTOR<VECTOR<T,3> ,2> > image;
    IMAGE<T>::Read(input_filename,image);
    GRID_2D<T> image_grid(image.m,image.n,0,image.m,0,image.n);image_grid.Set_MAC_Grid();

    if(parse_args.Is_Value_Set("-o")) output_filename=parse_args.Get_String_Value("-o");
    else output_filename=FILE_UTILITIES::Get_Basename(input_filename)+".phi";

    T scale=(T)parse_args.Get_Double_Value("-scale");
    GRID_2D<T> grid(image.m,image.n,-(T).5*scale*(image.m-1),(T).5*scale*(image.m-1),-(T).5*scale*(image.n-1),(T).5*scale*(image.n-1));
    ARRAYS<VECTOR<T,2> > phi(grid);LEVELSET_2D<T> levelset(grid,phi);

    LINEAR_INTERPOLATION<T,VECTOR<T,3> > color_interpolation;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){phi(i,j)=2*image(i,j).Magnitude()/root_three-(T)1;}

    if(parse_args.Get_Option_Value("-fmm")) levelset.Fast_Marching_Method(0,10*grid.max_dx_dy);

    FILE_UTILITIES::Write_To_File<T>(output_filename,levelset);
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Double_Argument("-scale",1);
    parse_args.Add_Double_Argument("-resolution_scale",1);
    parse_args.Add_String_Argument("-o", "");
    parse_args.Add_Option_Argument("-fmm");
    parse_args.Set_Extra_Arguments(-1, "<filename>");

    parse_args.Parse(argc, argv);

    if(parse_args.Num_Extra_Args() < 1) return 1;

    if(!parse_args.Get_Option_Value("-double")) Convert<float>(parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double>(parse_args);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
