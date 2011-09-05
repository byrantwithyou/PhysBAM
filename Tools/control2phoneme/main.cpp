//#####################################################################
// Copyright 2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Dynamics/Motion/PHONEME.h>
using namespace PhysBAM;
//#####################################################################
// Utility Functions
//#####################################################################
template <class T> int
main_templatized(int argc, char**argv)
{
    LOG::Initialize_Logging();

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-i","","input control prefix","input control prefix");
    parse_args.Add_String_Argument("-o","","output filename","output filename");
    parse_args.Add_String_Argument("-phoneme","","phoneme name","phoneme name");
    parse_args.Add_String_Argument("-previous_phoneme","","previous phoneme name","previous phoneme name");
    parse_args.Add_String_Argument("-next_phoneme","","next phoneme name","next phoneme name");
    parse_args.Add_Integer_Argument("-sample_start",0,"first control sample","first control sample");
    parse_args.Add_Integer_Argument("-sample_end",0,"last control sample","last control sample");
    parse_args.Add_Integer_Argument("-phoneme_start",0,"first sample of phoneme","first sample of phoneme");
    parse_args.Add_Integer_Argument("-phoneme_end",0);
    parse_args.Add_Double_Argument("-frame_rate",120.0);
    parse_args.Parse(argc,argv);

    std::string input_prefix;if(parse_args.Is_Value_Set("-i"))input_prefix=parse_args.Get_String_Value("-i");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    std::string output_filename;if(parse_args.Is_Value_Set("-o"))output_filename=parse_args.Get_String_Value("-o");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    std::string phoneme_name;if(parse_args.Is_Value_Set("-phoneme"))phoneme_name=parse_args.Get_String_Value("-phoneme");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    std::string previous_phoneme_name;if(parse_args.Is_Value_Set("-previous_phoneme"))previous_phoneme_name=parse_args.Get_String_Value("-previous_phoneme");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    std::string next_phoneme_name;if(parse_args.Is_Value_Set("-next_phoneme"))next_phoneme_name=parse_args.Get_String_Value("-next_phoneme");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    int sample_start;if(parse_args.Is_Value_Set("-sample_start"))sample_start=parse_args.Get_Integer_Value("-sample_start");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    int sample_end;if(parse_args.Is_Value_Set("-sample_end"))sample_end=parse_args.Get_Integer_Value("-sample_end");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    int phoneme_start;if(parse_args.Is_Value_Set("-phoneme_start"))phoneme_start=parse_args.Get_Integer_Value("-phoneme_start");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    int phoneme_end;if(parse_args.Is_Value_Set("-phoneme_end"))phoneme_end=parse_args.Get_Integer_Value("-phoneme_end");else {LOG::cerr<<"Incorrect arguments"<<std::endl;exit(0);}
    T frame_rate=(T)parse_args.Get_Double_Value("-frame_rate");

    LOG::Time("Initializing");
    PHONEME<T> phoneme;
    phoneme.Set_Names(phoneme_name,previous_phoneme_name,next_phoneme_name);
    phoneme.Set_Frame_Length(phoneme_end-phoneme_start+1);
    phoneme.Set_Time_Length(((T)phoneme_end-phoneme_start)/frame_rate);
    phoneme.Resize_Controls(phoneme_start-sample_start,sample_end-phoneme_end);
    LOG::cout<<"Phoneme name          : "<<phoneme.name<<std::endl;
    LOG::cout<<"Previous phoneme name : "<<phoneme.previous_name<<std::endl;
    LOG::cout<<"Next phoneme name     : "<<phoneme.next_name<<std::endl;
    LOG::cout<<"Frame length          : "<<phoneme.frame_length<<std::endl;
    LOG::cout<<"Time length           : "<<phoneme.time_length<<std::endl;
    LOG::Stop_Time();

    LOG::Time("Reading control frames");
    for(int frame=sample_start;frame<=sample_end;frame++)FILE_UTILITIES::Read_From_File<T>(input_prefix+STRING_UTILITIES::string_sprintf(".%d",frame),phoneme.Frame_Controls(frame-phoneme_start+1));
    LOG::Stop_Time();

    LOG::Time("Writing phoneme");
    FILE_UTILITIES::Write_To_File<T>(output_filename,phoneme);
    LOG::Stop_Time();

    LOG::Finish_Logging();
    return 0;
}
int main(int argc,char **argv)
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    return main_templatized<float>(argc,argv);
}
//#######################################################################
