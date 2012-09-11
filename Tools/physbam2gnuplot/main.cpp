//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/GNUPLOT_OUTPUT.h>
using namespace PhysBAM;
//#####################################################################
// Class PHYSBAM_TO_GNUPLOT_CONVERTER 
//#####################################################################
template<class T_GRID,class RW>
class PHYSBAM_TO_GNUPLOT_CONVERTER
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;

    const std::string input_directory,output_directory;
    GNUPLOT_OUTPUT gnuplot_output;
public:
    bool convert_density,convert_momentum,convert_energy,convert_velocity,convert_pressure,convert_internal_energy,convert_entropy,convert_machnumber,convert_log;

    PHYSBAM_TO_GNUPLOT_CONVERTER(const std::string& input_directory_input,const std::string& output_directory_input)
        :input_directory(input_directory_input),output_directory(output_directory_input)
    {}

    void Convert(const int frame,const std::string* variable_name=0)
    {T_GRID grid;FILE_UTILITIES::Read_From_File<RW>(input_directory+"/common/grid",grid);
    if(variable_name) Convert_Data<T_ARRAYS_SCALAR>(*variable_name,grid,frame);
    if(convert_density) Convert_Data<T_ARRAYS_SCALAR>("density",grid,frame);
    if(convert_momentum) Convert_Data<T_ARRAYS_SCALAR>("momentum",grid,frame);
    if(convert_energy) Convert_Data<T_ARRAYS_SCALAR>("energy",grid,frame);
    if(convert_velocity) Convert_Data<T_ARRAYS_SCALAR>("centered_velocities",grid,frame);
    if(convert_pressure) Convert_Data<T_ARRAYS_SCALAR>("pressure",grid,frame);
    if(convert_internal_energy) Convert_Data<T_ARRAYS_SCALAR>("internal_energy",grid,frame);
    if(convert_entropy) Convert_Data<T_ARRAYS_SCALAR>("entropy",grid,frame);
    if(convert_machnumber) Convert_Data<T_ARRAYS_SCALAR>("machnumber",grid,frame);
    }

    void Convert_All_Frames(const int first_frame,const int last_frame,const std::string* variable_name=0)
    {for(int frame=first_frame;frame<=last_frame;frame++){Convert(frame,variable_name);}}

private:
    template <class T_TYPE> void Convert_Data(const std::string& file_name_prefix,const T_GRID& grid,const int frame)
    {std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    std::string input_file=input_directory+"/"+f+"/"+file_name_prefix;
    std::string output_file=output_directory+"/"+file_name_prefix;
    T_TYPE data;FILE_UTILITIES::Read_From_File<RW>(input_file,data);
    if(convert_log){
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){T tmp=data(iterator.Cell_Index());
            data(iterator.Cell_Index()) = log(tmp)/log((T)10);}}
    gnuplot_output.Write_Output_File(output_file,grid,data,frame);}
};
template<class T_GRID,class RW> void PhysBAM_To_Gnuplot(const PARSE_ARGS& parse_args,const int verbosity)
{
    bool convert_density=parse_args.Is_Value_Set("-density");
    bool convert_momentum=parse_args.Is_Value_Set("-momentum");
    bool convert_energy=parse_args.Is_Value_Set("-energy"); 
    bool convert_velocity=parse_args.Is_Value_Set("-velocity"); 
    bool convert_pressure=parse_args.Is_Value_Set("-pressure"); 
    bool convert_internal_energy=parse_args.Is_Value_Set("-internal_energy"); 
    bool convert_entropy=parse_args.Is_Value_Set("-entropy"); 
    bool convert_machnumber=parse_args.Is_Value_Set("-machnumber"); 
    bool convert_log=parse_args.Is_Value_Set("-log"); 

    std::string input_directory=parse_args.Extra_Arg(0),output_directory=input_directory;
    if(parse_args.Is_Value_Set("-o")) output_directory=parse_args.Get_String_Value("-o");
    FILE_UTILITIES::Create_Directory(output_directory);

    const std::string* variable_name=(parse_args.Is_Value_Set("-v"))?&parse_args.Get_String_Value("-v"):0;

    int first_frame,last_frame;
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/first_frame",first_frame);
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);
    if(parse_args.Is_Value_Set("-start_frame")) first_frame=parse_args.Get_Integer_Value("-start_frame");
    if(parse_args.Is_Value_Set("-last_frame")) last_frame=parse_args.Get_Integer_Value("-last_frame");

    if(verbosity>0){
        std::cout<<"input_directory="<<input_directory<<", output_directory="<<output_directory<<std::endl;
        std::cout<<"first_frame="<<first_frame<<std::endl<<"last_frame="<<last_frame<<std::endl;
        if(variable_name) std::cout<<"variable_name="<<*variable_name<<std::endl;
        else std::cout<<"no variable_name specified"<<std::endl;}


    PHYSBAM_TO_GNUPLOT_CONVERTER<T_GRID,typename T_GRID::SCALAR> physbam_to_matlab_converter(input_directory,output_directory);

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
    parse_args.Add_Integer_Argument("-start_frame",0,"start frame number");
    parse_args.Add_Integer_Argument("-last_frame",0,"last frame number");
    parse_args.Add_Integer_Argument("-verbosity",0,"Verbosity level");
    parse_args.Add_String_Argument("-o","","output directory");
    parse_args.Add_Integer_Argument("-dimension",1,"Grid dimension");
    parse_args.Add_String_Argument("-v","","variable to read");
    parse_args.Add("-double",&use_double,"Read in file in double format");
    parse_args.Add_Option_Argument("-density","convert density");
    parse_args.Add_Option_Argument("-log","output log (base 10) of the data");
    parse_args.Add_Option_Argument("-momentum","convert momentum");
    parse_args.Add_Option_Argument("-machnumber","convert machnumber");
    parse_args.Add_Option_Argument("-energy","convert energy");
    parse_args.Add_Option_Argument("-entropy","convert entropy");
    parse_args.Add_Option_Argument("-velocity","convert velocity");
    parse_args.Add_Option_Argument("-pressure","convert pressure");
    parse_args.Add_Option_Argument("-internal_energy","convert internal_energy");
    parse_args.Set_Extra_Arguments(1,"<input_directory>");
    parse_args.Parse();

    int verbosity=parse_args.Get_Integer_Value("-verbosity");
    int dimension=parse_args.Get_Integer_Value("-dimension");
    
#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(use_double) PHYSBAM_FATAL_ERROR("No double support");
#endif

    if(verbosity>0) std::cout<<"dimension="<<dimension<<std::endl;
    if(use_double){
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        switch(dimension){
            case 1:
                PhysBAM_To_Gnuplot<GRID<VECTOR<double,1> >,double>(parse_args,verbosity);
                break;
            case 2:
                PhysBAM_To_Gnuplot<GRID<VECTOR<double,2> >,double>(parse_args,verbosity);
                break;
            case 3:
                PhysBAM_To_Gnuplot<GRID<VECTOR<double,3> >,double>(parse_args,verbosity);
                break;}
#endif
    }
    else{
        switch(dimension){
            case 1:
                PhysBAM_To_Gnuplot<GRID<VECTOR<float,1> >,float>(parse_args,verbosity);
                break;
            case 2:
                PhysBAM_To_Gnuplot<GRID<VECTOR<float,2> >,float>(parse_args,verbosity);
                break;
            case 3:
                PhysBAM_To_Gnuplot<GRID<VECTOR<float,3> >,float>(parse_args,verbosity);
                break;}}
}
