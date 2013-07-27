//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Read_Write/GNUPLOT_OUTPUT.h>
using namespace PhysBAM;
//#####################################################################
// Class PHYSBAM_TO_GNUPLOT_CONVERTER 
//#####################################################################
template<class TV,class RW>
class PHYSBAM_TO_GNUPLOT_CONVERTER
{
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::SCALAR T;

    const std::string input_directory,output_directory;
    GNUPLOT_OUTPUT gnuplot_output;
public:
    bool convert_density,convert_momentum,convert_energy,convert_velocity,convert_pressure,convert_internal_energy,convert_entropy,convert_machnumber,convert_log;

    PHYSBAM_TO_GNUPLOT_CONVERTER(const std::string& input_directory_input,const std::string& output_directory_input)
        :input_directory(input_directory_input),output_directory(output_directory_input)
    {}

    void Convert(const int frame,const std::string variable_name="")
    {GRID<TV> grid;FILE_UTILITIES::Read_From_File<RW>(input_directory+"/common/grid",grid);
    if(variable_name.size()) Convert_Data<ARRAY<T,TV_INT>>(variable_name,grid,frame);
    if(convert_density) Convert_Data<ARRAY<T,TV_INT>>("density",grid,frame);
    if(convert_momentum) Convert_Data<ARRAY<T,TV_INT>>("momentum",grid,frame);
    if(convert_energy) Convert_Data<ARRAY<T,TV_INT>>("energy",grid,frame);
    if(convert_velocity) Convert_Data<ARRAY<T,TV_INT>>("centered_velocities",grid,frame);
    if(convert_pressure) Convert_Data<ARRAY<T,TV_INT>>("pressure",grid,frame);
    if(convert_internal_energy) Convert_Data<ARRAY<T,TV_INT>>("internal_energy",grid,frame);
    if(convert_entropy) Convert_Data<ARRAY<T,TV_INT>>("entropy",grid,frame);
    if(convert_machnumber) Convert_Data<ARRAY<T,TV_INT>>("machnumber",grid,frame);
    }

    void Convert_All_Frames(const int first_frame,const int last_frame,const std::string variable_name="")
    {for(int frame=first_frame;frame<=last_frame;frame++){Convert(frame,variable_name);}}

private:
    template <class T_TYPE> void Convert_Data(const std::string& file_name_prefix,const GRID<TV>& grid,const int frame)
    {std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    std::string input_file=input_directory+"/"+f+"/"+file_name_prefix;
    std::string output_file=output_directory+"/"+file_name_prefix;
    T_TYPE data;FILE_UTILITIES::Read_From_File<RW>(input_file,data);
    if(convert_log){
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){T tmp=data(iterator.Cell_Index());
            data(iterator.Cell_Index()) = log(tmp)/log((T)10);}}
    gnuplot_output.Write_Output_File(output_file,grid,data,frame);}
};
template<class TV,class RW> void PhysBAM_To_Gnuplot(PARSE_ARGS& parse_args,const int verbosity)
{
    bool convert_density=false,convert_momentum=false,convert_energy=false,convert_velocity=false;
    bool convert_pressure=false,convert_internal_energy=false,convert_entropy=false,convert_machnumber=false,convert_log=false;
    int first_frame=0,last_frame=0;
    std::string output_directory,input_directory,variable_name;
    parse_args.Add("-density",&convert_density,"convert density");
    parse_args.Add("-log",&convert_log,"output log (base 10) of the data");
    parse_args.Add("-momentum",&convert_momentum,"convert momentum");
    parse_args.Add("-machnumber",&convert_machnumber,"convert machnumber");
    parse_args.Add("-energy",&convert_energy,"convert energy");
    parse_args.Add("-entropy",&convert_entropy,"convert entropy");
    parse_args.Add("-velocity",&convert_velocity,"convert velocity");
    parse_args.Add("-pressure",&convert_pressure,"convert pressure");
    parse_args.Add("-internal_energy",&convert_internal_energy,"convert internal_energy");
    parse_args.Add("-start_frame",&first_frame,"frame","start frame number");
    parse_args.Add("-last_frame",&last_frame,"frame","last frame number");
    parse_args.Add("-o",&output_directory,"file","output directory");
    parse_args.Add("-v",&variable_name,"var","variable to read");
    parse_args.Extra(&input_directory,"input_directory","input_directory");
    parse_args.Parse();

    if(!output_directory.size()) output_directory=input_directory;
    FILE_UTILITIES::Create_Directory(output_directory);

    FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/first_frame",first_frame);
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);

    if(verbosity>0){
        std::cout<<"input_directory="<<input_directory<<", output_directory="<<output_directory<<std::endl;
        std::cout<<"first_frame="<<first_frame<<std::endl<<"last_frame="<<last_frame<<std::endl;
        if(variable_name.size()) std::cout<<"variable_name="<<variable_name<<std::endl;
        else std::cout<<"no variable_name specified"<<std::endl;}

    PHYSBAM_TO_GNUPLOT_CONVERTER<TV,typename TV::SCALAR> physbam_to_matlab_converter(input_directory,output_directory);

    physbam_to_matlab_converter.convert_density=convert_density;
    physbam_to_matlab_converter.convert_momentum=convert_momentum;
    physbam_to_matlab_converter.convert_energy=convert_energy;
    physbam_to_matlab_converter.convert_velocity=convert_velocity;
    physbam_to_matlab_converter.convert_pressure=convert_pressure;
    physbam_to_matlab_converter.convert_internal_energy=convert_internal_energy;
    physbam_to_matlab_converter.convert_entropy=convert_entropy;
    physbam_to_matlab_converter.convert_machnumber=convert_machnumber;
    physbam_to_matlab_converter.convert_log=convert_log;

    physbam_to_matlab_converter.Convert_All_Frames(first_frame,last_frame,variable_name);
}
int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    bool use_double=false;
    int verbosity=0;
    int dimension=1;
    parse_args.Add("-verbosity",&verbosity,"value","Verbosity level");
    parse_args.Add("-dimension",&dimension,"value","Grid dimension");
    parse_args.Add("-double",&use_double,"Read in file in double format");
    parse_args.Parse(true);

#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(use_double) PHYSBAM_FATAL_ERROR("No double support");
#endif

    if(verbosity>0) std::cout<<"dimension="<<dimension<<std::endl;
    if(use_double){
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        switch(dimension){
            case 1:
                PhysBAM_To_Gnuplot<VECTOR<double,1>,double>(parse_args,verbosity);
                break;
            case 2:
                PhysBAM_To_Gnuplot<VECTOR<double,2>,double>(parse_args,verbosity);
                break;
            case 3:
                PhysBAM_To_Gnuplot<VECTOR<double,3>,double>(parse_args,verbosity);
                break;}
#endif
    }
    else{
        switch(dimension){
            case 1:
                PhysBAM_To_Gnuplot<VECTOR<float,1>,float>(parse_args,verbosity);
                break;
            case 2:
                PhysBAM_To_Gnuplot<VECTOR<float,2>,float>(parse_args,verbosity);
                break;
            case 3:
                PhysBAM_To_Gnuplot<VECTOR<float,3>,float>(parse_args,verbosity);
                break;}}
}
