//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARSE_ARGS
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <cstring>
using namespace PhysBAM;
//#####################################################################
// Function Use_Help_Option
//#####################################################################
PARSE_ARGS::
PARSE_ARGS(int argc_input,char** argv_input)
    :num_expected_extra_args(0),use_help_option(true),extra_usage_callback(0),argc(argc_input),argv(argv_input)
{
}
PARSE_ARGS::
~PARSE_ARGS()
{
}
void PARSE_ARGS::Use_Help_Option(bool use_it)
{
    use_help_option=use_it;
}
//#####################################################################
// Function Add_Option_Argument
//#####################################################################
void PARSE_ARGS::
Add_Option_Argument(const std::string& arg_str,const std::string& desc)
{
    arg_data_list.Append(ARG_DATA(arg_str,desc));
}
//#####################################################################
// Function Add_Integer_Argument
//#####################################################################
void PARSE_ARGS::
Add_Integer_Argument(const std::string& arg_str,int default_value,const std::string& val_name,const std::string& desc)
{
    arg_data_list.Append(ARG_DATA(arg_str,val_name,desc,default_value));
}
//#####################################################################
// Function Add_Double_Argument
//#####################################################################
void PARSE_ARGS::
Add_Double_Argument(const std::string& arg_str,double default_value,const std::string& val_name,const std::string& desc)
{
    arg_data_list.Append(ARG_DATA(arg_str,val_name,desc,default_value));
}
//#####################################################################
// Function Add_Vector_2D_Argument
//#####################################################################
void PARSE_ARGS::
Add_Vector_2D_Argument(const std::string& arg_str,const VECTOR<double,2> &default_value,const std::string& val_name,const std::string& desc)
{
    arg_data_list.Append(ARG_DATA(arg_str,val_name,desc,default_value));
}
//#####################################################################
// Function Add_Vector_3D_Argument
//#####################################################################
void PARSE_ARGS::
Add_Vector_3D_Argument(const std::string& arg_str,const VECTOR<double,3> &default_value,const std::string& val_name,const std::string& desc)
{
    arg_data_list.Append(ARG_DATA(arg_str,val_name,desc,default_value));
}
//#####################################################################
// Function Add_String_Argument
//#####################################################################
void PARSE_ARGS::
Add_String_Argument(const std::string& arg_str,const std::string& default_value,const std::string& val_name,const std::string& desc)
{
    arg_data_list.Append(ARG_DATA(arg_str,val_name,desc,default_value));
}
//#####################################################################
// Function Set_Extra_Arguments
//#####################################################################
void PARSE_ARGS::
Set_Extra_Arguments(int num,const std::string& synopsis,const std::string& desc) // num=-1 for arbitrary extra rguments
{
    num_expected_extra_args=num;if(synopsis.length())extra_args_synopsis=synopsis;if(desc.length())extra_args_desc=desc;
}
//#####################################################################
// Function Set_Extra_Usage_Callback
//#####################################################################
void PARSE_ARGS::
Set_Extra_Usage_Callback(void (*extra_usage_callback_input)())
{
    extra_usage_callback=extra_usage_callback_input;
}
//#####################################################################
// Function Parse
//#####################################################################
void PARSE_ARGS::
Parse(bool partial)
{
    program_name=argv[0];
    extra_arg_list.Remove_All();
    int kept=1;
    for(int i=1;i<argc;i++){
        if(OPTION* o=options.Get_Pointer(argv[i])){
            if(o->found) *o->found=true;
            if(o->store)
                if(!argv[++i] || !o->store_func(o->store,argv[i]))
                    Print_Usage(true);}
        else{
            int match=Find_Match(argv[i]);
            if(match>=0){
                if(!arg_data_list(match).Parse_Value(argc,argv,i))
                    Print_Usage(true);
                else arg_data_list(match).value_set=true;
                i--;}
            else argv[kept++]=argv[i];}}
    argc=kept;
    argv[argc]=0;

    if(!partial){
        for(int i=1;i<argc;i++){
            if(argv[i][0]=='-' && isalpha(argv[i][1])) Print_Usage(true);
            extra_arg_list.Append(argv[i]);}
        if(num_expected_extra_args!=-1 && extra_arg_list.m<num_expected_extra_args) Print_Usage(true);} // didn't get the expected number of extra args
}
//#####################################################################
// Function Get_Option_Value
//#####################################################################
bool PARSE_ARGS::
Get_Option_Value(const std::string& arg_str) const
{
    return arg_data_list(Find_Match(arg_str,ARG_DATA::OPTION)).boolean_value;
}
//#####################################################################
// Function Get_Integer_Value
//#####################################################################
int PARSE_ARGS::
Get_Integer_Value(const std::string& arg_str) const
{
    return arg_data_list(Find_Match(arg_str,ARG_DATA::INTEGER)).integer_value;
}
//#####################################################################
// Function Get_Double_Value
//#####################################################################
double PARSE_ARGS::
Get_Double_Value(const std::string& arg_str) const
{
    return arg_data_list(Find_Match(arg_str,ARG_DATA::DOUBLE)).double_value;
}
//#####################################################################
// Function Get_Vector_2D_Value
//#####################################################################
VECTOR<double,2> PARSE_ARGS::
Get_Vector_2D_Value(const std::string& arg_str) const
{
    return arg_data_list(Find_Match(arg_str,ARG_DATA::VECTOR2)).vector_2d_value;
}
//#####################################################################
// Function Get_Vector_3D_Value
//#####################################################################
VECTOR<double,3> PARSE_ARGS::
Get_Vector_3D_Value(const std::string& arg_str) const
{
    return arg_data_list(Find_Match(arg_str,ARG_DATA::VECTOR3)).vector_3d_value;
}
//#####################################################################
// Function Get_String_Value
//#####################################################################
const std::string& PARSE_ARGS::
Get_String_Value(const std::string& arg_str) const
{
    return arg_data_list(Find_Match(arg_str,ARG_DATA::STRING)).string_value;
}
//#####################################################################
// Function Is_Value_Set
//#####################################################################
bool PARSE_ARGS::
Is_Value_Set(const std::string& arg_str) const
{
    int match=Find_Match(arg_str);
    if(match==-1){LOG::cout<<"Argument "<<arg_str<<" undeclared"<<std::endl;PHYSBAM_FATAL_ERROR();}
    return arg_data_list(match).value_set;
}
//#####################################################################
// Function Override_String_Value
//#####################################################################
void PARSE_ARGS::
Override_String_Value(const std::string& arg_str,const std::string& value)
{
    arg_data_list(Find_Match(arg_str,ARG_DATA::STRING)).string_value=value;
}
//#####################################################################
// Function Find_Match
//#####################################################################
int PARSE_ARGS::
Find_Match(const std::string& str) const
{
    for(int i=0;i<arg_data_list.m;i++) if(arg_data_list(i).str==str) return i;
    return -1;
}
//#####################################################################
// Function Find_Match
//#####################################################################
int PARSE_ARGS::
Find_Match(const std::string& str,const ARG_DATA::TYPE& type) const
{
    int match=Find_Match(str);
    if(match==-1){LOG::cout<<"Argument "<<str<<" undeclared"<<std::endl;PHYSBAM_FATAL_ERROR();}
    if(arg_data_list(match).type!=type){LOG::cout<<"Type mismatch in Find_Match("<<str<<")"<<std::endl;PHYSBAM_FATAL_ERROR();}
    return match;
}
//#####################################################################
// Function Num_Extra_Args
//#####################################################################
int PARSE_ARGS::
Num_Extra_Args() const
{
    return extra_arg_list.m;
}
//#####################################################################
// Function Extra_Arg
//#####################################################################
const std::string& PARSE_ARGS::
Extra_Arg(int i) const
{
    PHYSBAM_ASSERT((unsigned)i<(unsigned)extra_arg_list.m);
    return extra_arg_list(i);
}
//#####################################################################
// Function Get_Program_Name
//#####################################################################
const std::string& PARSE_ARGS::
Get_Program_Name() const
{
    return program_name;
}
//#####################################################################
// Function Find_And_Remove
//#####################################################################
bool PARSE_ARGS::
Find_And_Remove(const char *str,int& argc,char** argv)
{
    int i;for(i=0;i<argc;i++)if(!strcmp(str,argv[i]))break;
    if(i<argc){for(;i<argc-1;i++)argv[i]=argv[i+1];argc--;argv[argc]=0;return true;}
    return false;
}
//#####################################################################
// Function Find_And_Remove_Double
//#####################################################################
double PARSE_ARGS::
Find_And_Remove_Double(const char *str,int& argc,char** argv)
{
    double value;
    int i;for(i=0;i<argc;i++)if(!strcmp(str,argv[i]))break;
    if(i+1<argc){value=atof(argv[i+1]);for(;i<argc-2;i++)argv[i]=argv[i+2];argc--;argv[argc]=0;argc--;argv[argc]=0;return value;}
    return 0.0;
}
//#####################################################################
// Function Find_And_Remove_Integer
//#####################################################################
int PARSE_ARGS::
Find_And_Remove_Integer(const char *str,int& argc,char** argv)
{
    int value;
    int i;for(i=0;i<argc;i++)if(!strcmp(str,argv[i]))break;
    if(i+1<argc){value=atoi(argv[i+1]);for(;i<argc-2;i++)argv[i]=argv[i+2];argc--;argv[argc]=0;argc--;argv[argc]=0;return value;}
    return 0;
}
//#####################################################################
// Function Print_Usage
//#####################################################################
void PARSE_ARGS::
Print_Usage(bool do_exit) const
{
    ARRAY<std::string> args;
    options.Get_Keys(args);
    args.Sort();

    LOG::cerr<<"Usage: "<<program_name;
    for(int i=0;i<arg_data_list.m;i++){LOG::cerr<<" ";arg_data_list(i).Print_Synopsis();}
    for(int i=0;i<args.m;i++){
        const OPTION& o=options.Get(args(i));
        LOG::cerr<<" ["<<o.opt;
        if(o.store) LOG::cerr<<" <"<<(o.name.size()?o.name:"arg")<<">";
        LOG::cerr<<"]";}
    LOG::cerr<<extra_args_synopsis<<std::endl;

    int width=0;
    for(int i=0;i<arg_data_list.m;i++){int len=(int)arg_data_list(i).str.length();if(len>width)width=len;}
    for(int i=0;i<args.m;i++) width=max((int)args(i).size(),width);

    for(int i=0;i<arg_data_list.m;i++){arg_data_list(i).Print_Description(width+2);LOG::cerr<<std::endl;}
    for(int i=0;i<args.m;i++){
        const OPTION& o=options.Get(args(i));
        LOG::cerr.flags(std::ios::left);
        LOG::cerr.width(width+2);
        LOG::cerr<<o.opt<<o.desc;
        if(o.store){
            LOG::cerr<<" (";
            o.print_default_func(o.store);
            LOG::cerr<<")";}
        LOG::cerr<<std::endl;}
    if(extra_args_desc.size()) LOG::cerr<<" "<<extra_args_desc<<std::endl;
    if(extra_usage_callback) extra_usage_callback();
    if(do_exit) exit(-1);
}
//#####################################################################
// Function Print_Arguments
//#####################################################################
std::string PARSE_ARGS::
Print_Arguments() const
{
    std::string s="command = ";
    for(int i=0;i<argc;i++){s+=argv[i];s+=' ';}
    s+="\nworking directory = "+FILE_UTILITIES::Get_Working_Directory()+"\n";
    return s;
}
//#####################################################################
