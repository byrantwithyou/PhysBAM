//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/EXR_FILE.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class T2>
void Resample(PARSE_ARGS &parse_args)
{
    std::string input_filename_pattern=parse_args.Extra_Arg(1);
    std::string output_filename_pattern;
    if (parse_args.Is_Value_Set("-o")) output_filename_pattern=parse_args.Get_String_Value("-o");
    else output_filename_pattern="resampled_"+input_filename_pattern;

    GRID_2D<T> original_grid;
    FILE_UTILITIES::Read_From_File<T>("grid",original_grid);

    std::cout << "Got original grid " << original_grid << std::endl;

    int m=parse_args.Get_Integer_Value("-m"),n=parse_args.Get_Integer_Value("-n");
    T initial_xmin=parse_args.Get_Double_Value("-xmin"),initial_xmax=parse_args.Get_Double_Value("-xmax");
    T initial_ymin=parse_args.Get_Double_Value("-ymin"),initial_ymax=parse_args.Get_Double_Value("-ymax");
    VECTOR<T,2> initial_minimum_corner(initial_xmin,initial_ymin),initial_maximum_corner(initial_xmax,initial_ymax);
    int fps=parse_args.Get_Integer_Value("-fps");
    VECTOR<T,2> velocity=parse_args.Get_Vector_2D_Value("-velocity");

    if(m==0 || n==0 || (initial_xmax-initial_xmin)<=0 || (initial_ymax-initial_ymin)<=0) {
        std::cerr << "Invalid resampled grid dimensions" << std::endl;
        exit(1);
    }

    GRID_2D<T> resampled_grid;
    ARRAYS<VECTOR<T2,2> > original_field;
    ARRAYS<VECTOR<T2,2> > resampled_field(1,m,1,n);

    LINEAR_INTERPOLATION<T,T2> interpolation;

    int start_frame=parse_args.Get_Integer_Value("-start_frame");
    int end_frame=parse_args.Get_Integer_Value("-end_frame");
    for(int frame=start_frame;frame<=end_frame;frame++)
    {
        T time=(T)frame/fps;

        std::string input_filename=FILE_UTILITIES::Get_Frame_Filename(input_filename_pattern,frame);
        if(!FILE_UTILITIES::File_Exists(input_filename)) break;

        std::cout << "Frame " << frame << ", Time " << time << ", Grid " << resampled_grid << std::endl;
        VECTOR<T,2> minimum_corner=initial_minimum_corner+time*velocity,maximum_corner=initial_maximum_corner+time*velocity;
        resampled_grid=GRID_2D<T>(m,n,minimum_corner.x,maximum_corner.x,minimum_corner.y,maximum_corner.y);

        std::cout << "Reading..." << std::flush;
        FILE_UTILITIES::Read_From_File<T>(input_filename,original_field);

        std::cout << "Resampling..." << std::flush;
        for(int i=1;i<=resampled_grid.m;i++) for(int j=1;j<=resampled_grid.n;j++)
            resampled_field(i,j)=interpolation.Clamped_To_Array(original_grid,original_field,resampled_grid.X(i,j));

        std::cout << "Writing" << std::endl;
        FILE_UTILITIES::Write_To_File<T>(STRING_UTILITIES::string_sprintf("resampled_grid.%d",frame),resampled_grid);

        std::string output_filename=FILE_UTILITIES::Get_Frame_Filename(output_filename_pattern,frame);
        FILE_UTILITIES::Write_To_File<T>(output_filename,resampled_field);
    }
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_String_Argument("-o", "", "output filename pattern");
    parse_args.Add_Integer_Argument("-m",0);
    parse_args.Add_Integer_Argument("-n",0);
    parse_args.Add_Double_Argument("-xmin",0);
    parse_args.Add_Double_Argument("-xmax",0);
    parse_args.Add_Double_Argument("-ymin",0);
    parse_args.Add_Double_Argument("-ymax",0);
    parse_args.Add_Vector_2D_Argument("-velocity",VECTOR<double,2>(0,0));
    parse_args.Add_Integer_Argument("-fps",24);
    parse_args.Add_Integer_Argument("-start_frame",0);
    parse_args.Add_Integer_Argument("-end_frame",1000000);
    parse_args.Set_Extra_Arguments(-1, "<filename_pattern>");

    parse_args.Parse(argc, argv);

    if (parse_args.Num_Extra_Args() < 1) return 1;

    if(!parse_args.Get_Option_Value("-double")) Resample<float,float>(parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Resample<double,double>(parse_args);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
