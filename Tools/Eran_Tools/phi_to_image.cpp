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
    std::string output_filename;
    int scale=1,samples=2;
    VECTOR<T,3> negative_color(0,0,1),positive_color(1,1,1);
    parse_args.Add("-scale",&scale,"scale","image scale");
    parse_args.Add("-samples",&samples,"samples","samples");
    parse_args.Add("-negative_color",&negative_color,"negative_color","negative_color");
    parse_args.Add("-positive_color",&positive_color,"positive_color","positive_color");
    parse_args.Add("-o",&output_filename,"file","output filename");
    parse_args.Set_Extra_Arguments(-1, "<filename>");
    parse_args.Parse();

    if(parse_args.Num_Extra_Args() < 1) return;
    std::string input_filename=parse_args.Extra_Arg(0);

    if(output_filename.empty()) output_filename=FILE_UTILITIES::Get_Basename(input_filename)+".ppm";
    if(!IMAGE<T>::Is_Supported(output_filename)){std::cerr << "Image format for '" << output_filename << "' not supported" << std::endl;exit(1);}

    GRID_2D<T> grid;ARRAYS<VECTOR<T,2> > phi;LEVELSET_2D<T> levelset(grid,phi);
    FILE_UTILITIES::Read_From_File<T>(input_filename,levelset);

    GRID_2D<T> image_grid(grid.m*scale,grid.n*scale,grid.xmin-grid.dx/2,grid.xmax+grid.dx/2,grid.ymin-grid.dy/2,grid.ymax+grid.dy/2);
    image_grid.Set_MAC_Grid();
    ARRAYS<VECTOR<VECTOR<T,3> ,2> > image(image_grid);

    GRID_2D<T> sample_offset_grid(samples,samples,-image_grid.dx/2,image_grid.dx/2,-image_grid.dy/2,image_grid.dy/2);
    int subsamples=sample_offset_grid.m*sample_offset_grid.n;

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
    bool use_double=false,use_float=true;
    parse_args.Add_Option_Argument("-double",&use_double,"use doubles");
    parse_args.Add_Option_Argument("-float",&use_float,"use floats");
    parse_args.Parse(true);

    if(!use_doubles) Convert<float>(parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double>(parse_args);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
