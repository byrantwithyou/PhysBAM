#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T>
void Convert(PARSE_ARGS &parse_args)
{
    std::string input_filename=parse_args.Extra_Arg(1);
    std::string output_filename;

    if (parse_args.Is_Value_Set("-o")) output_filename=parse_args.Get_String_Value("-o");
    else output_filename=FILE_UTILITIES::Get_Basename(input_filename)+".ppm";
    if(!IMAGE<T>::Is_Supported(output_filename)){std::cerr << "Image format for '" << output_filename << "' not supported" << std::endl;exit(1);}

    GRID_2D<T> grid;ARRAYS<VECTOR<T,2> > phi;LEVELSET_2D<T> levelset(grid,phi);
    FILE_UTILITIES::Read_From_File<T>(input_filename,levelset);

    int scale=(int)parse_args.Get_Integer_Value("-scale");
    GRID_2D<T> image_grid(grid.m*scale,grid.n*scale,grid.xmin-grid.dx/2,grid.xmax+grid.dx/2,grid.ymin-grid.dy/2,grid.ymax+grid.dy/2);
    image_grid.Set_MAC_Grid();
    ARRAYS<VECTOR<VECTOR<T,3> ,2> > image(image_grid);

    int samples=parse_args.Get_Integer_Value("-samples");
    GRID_2D<T> sample_offset_grid(samples,samples,-image_grid.dx/2,image_grid.dx/2,-image_grid.dy/2,image_grid.dy/2);
    int subsamples=sample_offset_grid.m*sample_offset_grid.n;

    VECTOR<T,3> negative_color(parse_args.Get_Vector_3D_Value("-negative_color"));
    VECTOR<T,3> positive_color(parse_args.Get_Vector_3D_Value("-positive_color"));
    for(int i=0;i<image_grid.m;i++) for(int j=0;j<image_grid.n;j++){
        int in_count=0;
        for(int subi=0;subi<sample_offset_grid.m;subi++) for(int subj=0;subj<sample_offset_grid.n;subj++)
            if(levelset.Phi(image_grid.X(i,j)+sample_offset_grid.X(subi,subj))<=0) in_count++;
        T alpha=(T)in_count/subsamples;
        image(i,j)=alpha*negative_color+(1-alpha)*positive_color;}

    IMAGE<T>::Write(output_filename,image);
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Integer_Argument("-scale",1);
    parse_args.Add_Integer_Argument("-samples",2);
    parse_args.Add_Vector_3D_Argument("-negative_color",VECTOR<double,3>(0,0,1));
    parse_args.Add_Vector_3D_Argument("-positive_color",VECTOR<double,3>(1,1,1));
    parse_args.Add_String_Argument("-o", "");
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
