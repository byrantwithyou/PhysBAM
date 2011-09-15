//#####################################################################
// Copyright 2008, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Matlab/MATLAB_OUTPUT.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Class PHYSBAM_TO_MATLAB_CONVERTER 
//#####################################################################
template<class T_GRID,class RW>
class PHYSBAM_TO_MATLAB_CONVERTER
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;

    const std::string input_directory,output_directory;
    MATLAB_OUTPUT matlab_output;
public:
    bool convert_density,convert_momentum,convert_energy,convert_velocity,convert_pressure,convert_internal_energy,convert_entropy,convert_machnumber;

    PHYSBAM_TO_MATLAB_CONVERTER(const std::string& input_directory_input,const std::string& output_directory_input)
        :input_directory(input_directory_input),output_directory(output_directory_input),convert_density(false),convert_momentum(false),convert_energy(false),convert_velocity(false),
        convert_pressure(false),convert_internal_energy(false)
    {}

    void Convert(const int frame)
    {T_GRID grid;FILE_UTILITIES::Read_From_File<RW>(input_directory+"/common/grid",grid);
    matlab_output.Write_Header_File(output_directory+"/header",grid,frame);
    if(convert_density) Convert_Data<T_ARRAYS_SCALAR>("density",frame);
    if(convert_momentum) Convert_Data<T_ARRAYS_SCALAR>("momentum",frame);
    if(convert_energy) Convert_Data<T_ARRAYS_SCALAR>("energy",frame);
    if(convert_velocity) Convert_Data<T_ARRAYS_SCALAR>("centered_velocities",frame);
    if(convert_pressure) Convert_Data<T_ARRAYS_SCALAR>("pressure",frame);
    if(convert_internal_energy) Convert_Data<T_ARRAYS_SCALAR>("internal_energy",frame);
    if(convert_entropy) Convert_Data<T_ARRAYS_SCALAR>("entropy",frame);
    if(convert_machnumber) Convert_Data<T_ARRAYS_SCALAR>("machnumber",frame);
    }

    void Convert_All_Frames(const int first_frame,const int last_frame)
    {for(int frame=first_frame;frame<=last_frame;frame++){std::cout<<"frame="<<frame<<std::endl;Convert(frame);}}

private:
    template <class T_TYPE> void Convert_Data(const std::string& file_name_prefix,const int frame)
    {std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    std::string input_file=input_directory+"/"+f+"/"+file_name_prefix;
    std::string output_file=output_directory+"/"+file_name_prefix;
    T_TYPE data;FILE_UTILITIES::Read_From_File<RW>(input_file,data);
    matlab_output.Write_Output_File(output_file,data,frame);}
};
int main(int argc,char* argv[])
{
#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
    typedef float T;
#else
    typedef double T;
#endif
    typedef double RW;

    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-start_frame",0,"start frame number");
    parse_args.Add_Integer_Argument("-last_frame",0,"last frame number");
    parse_args.Add_String_Argument("-o","","output directory");
    parse_args.Add_Option_Argument("-density","convert density");
    parse_args.Add_Option_Argument("-momentum","convert momentum");
    parse_args.Add_Option_Argument("-machnumber","convert machnumber");
    parse_args.Add_Option_Argument("-energy","convert energy");
    parse_args.Add_Option_Argument("-entropy","convert entropy");
    parse_args.Add_Option_Argument("-velocity","convert velocity");
    parse_args.Add_Option_Argument("-pressure","convert pressure");
    parse_args.Add_Option_Argument("-internal_energy","convert internal_energy");
    parse_args.Set_Extra_Arguments(1,"<input_directory>");
    parse_args.Parse(argc,argv);

    std::string input_directory=parse_args.Extra_Arg(1),output_directory=input_directory;
    if(parse_args.Is_Value_Set("-o")) output_directory=parse_args.Get_String_Value("-o");
    FILE_UTILITIES::Create_Directory(output_directory);

    int first_frame,last_frame;
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/first_frame",first_frame);
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);
    if(parse_args.Is_Value_Set("-start_frame")) first_frame=parse_args.Get_Integer_Value("-start_frame");
    if(parse_args.Is_Value_Set("-last_frame")) last_frame=parse_args.Get_Integer_Value("-last_frame");

#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
    PHYSBAM_FATAL_ERROR("No double support");
#else
    bool convert_density=parse_args.Is_Value_Set("-density");
    bool convert_momentum=parse_args.Is_Value_Set("-momentum");
    bool convert_energy=parse_args.Is_Value_Set("-energy"); 
    bool convert_velocity=parse_args.Is_Value_Set("-velocity"); 
    bool convert_pressure=parse_args.Is_Value_Set("-pressure"); 
    bool convert_internal_energy=parse_args.Is_Value_Set("-internal_energy"); 
    bool convert_entropy=parse_args.Is_Value_Set("-entropy"); 
    bool convert_machnumber=parse_args.Is_Value_Set("-machnumber"); 

    std::cout<<"input_directory="<<input_directory<<"output_directory="<<output_directory<<std::endl;
    std::cout<<"first_frame="<<first_frame<<std::endl<<"last_frame="<<last_frame<<std::endl;
    PHYSBAM_TO_MATLAB_CONVERTER<GRID<VECTOR<T,1> >,RW> physbam_to_matlab_converter(input_directory,output_directory);
    physbam_to_matlab_converter.convert_density=convert_density;
    physbam_to_matlab_converter.convert_momentum=convert_momentum;
    physbam_to_matlab_converter.convert_energy=convert_energy;
    physbam_to_matlab_converter.convert_velocity=convert_velocity;
    physbam_to_matlab_converter.convert_pressure=convert_pressure;
    physbam_to_matlab_converter.convert_internal_energy=convert_internal_energy;
    physbam_to_matlab_converter.convert_entropy=convert_entropy;
    physbam_to_matlab_converter.convert_machnumber=convert_machnumber;

    physbam_to_matlab_converter.Convert_All_Frames(first_frame,last_frame);
#endif
}
