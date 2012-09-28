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
//#####################################################################
// Destructor
//#####################################################################
PARSE_ARGS::
~PARSE_ARGS()
{
}
//#####################################################################
// Function Use_Help_Option
//#####################################################################
void PARSE_ARGS::
Use_Help_Option(bool use_it)
{
    use_help_option=use_it;
}
//#####################################################################
// Function Set_Extra_Arguments
//#####################################################################
void PARSE_ARGS::
Set_Extra_Arguments(int num,const std::string& synopsis,const std::string& desc) // num=-1 for arbitrary extra rguments
{
    num_expected_extra_args=num;
    if(synopsis.length()) extra_args_synopsis=synopsis;
    if(desc.length()) extra_args_desc=desc;
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
            if(o->found) *o->found=o->found_value;
            if(o->store)
                if(!argv[++i] || !o->store_func(o->store,argv[i]))
                    Print_Usage(true);}
        else argv[kept++]=argv[i];}
    argc=kept;
    argv[argc]=0;

    if(!partial){
        for(int i=1;i<argc;i++){
            if(argv[i][0]=='-' && isalpha(argv[i][1])) Print_Usage(true);
            extra_arg_list.Append(argv[i]);}
        if(num_expected_extra_args!=-1 && extra_arg_list.m<num_expected_extra_args) Print_Usage(true);} // didn't get the expected number of extra args
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
// Function Print_Usage
//#####################################################################
void PARSE_ARGS::
Print_Usage(bool do_exit) const
{
    ARRAY<std::string> args;
    options.Get_Keys(args);
    args.Sort();

    LOG::cerr<<"Usage: "<<program_name;
    for(int i=0;i<args.m;i++){
        const OPTION& o=options.Get(args(i));
        LOG::cerr<<" ["<<o.opt;
        if(o.store) LOG::cerr<<" <"<<(o.name.size()?o.name:"arg")<<">";
        LOG::cerr<<"]";}
    LOG::cerr<<extra_args_synopsis<<std::endl;

    int width=0;
    for(int i=0;i<args.m;i++) width=max((int)args(i).size(),width);

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
