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
    :argc(argc_input),argv(argv_input),unclaimed_arguments(false)
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
// Function Parse
//#####################################################################
void PARSE_ARGS::
Parse(bool partial)
{
    program_name=argv[0];
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
        int k=0;
        for(int i=1;i<argc;i++){
            if(argv[i][0]=='-' && isalpha(argv[i][1])) Print_Usage(true);
            if(k<extras.m){
                if(extras(k).found) *extras(k).found=true;
                if(extras(k).store)
                    if(!extras(k).store_func(extras(k).store,argv[i]))
                        Print_Usage(true);
                if(!extras(k).exhaust) k++;}
            else unclaimed_arguments=true;}
        if(extras(k).exhaust) k++;
        for(;k<extras.m;k++)
            if(extras(k).required)
                Print_Usage(true);}
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

    for(int i=0;i<extras.m;i++) LOG::cerr<<(extras(i).required?" <":" [<")<<extras(i).name<<(extras(i).required?">":">]");

    int width=0;
    for(int i=0;i<args.m;i++) width=max((int)args(i).size(),width);
    for(int i=0;i<extras.m;i++) width=max((int)extras(i).name.size()+2,width);

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
    for(int i=0;i<extras.m;i++){
        LOG::cerr.flags(std::ios::left);
        LOG::cerr.width(width+2);
        LOG::cerr<<"<"<<extras(i).name<<">"<<extras(i).desc;
        if(!extras(i).required){
            LOG::cerr<<" (";
            extras(i).print_default_func(extras(i).store);
            LOG::cerr<<")";}
        LOG::cerr<<std::endl;}
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
