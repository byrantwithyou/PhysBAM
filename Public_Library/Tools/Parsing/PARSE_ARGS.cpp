//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARSE_ARGS
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
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
                if(!argv[++i] || !o->store_func(o->store,argv[i])){
                    LOG::cerr<<"Option '"<<argv[i]<<"' expects an argument."<<std::endl;
                    Print_Usage(true);}}
        else argv[kept++]=argv[i];}
    argc=kept;
    argv[argc]=0;

    if(!partial){
        int k=0;
        for(int i=1;i<argc;i++){
            if(argv[i][0]=='-' && isalpha(argv[i][1])){
                LOG::cerr<<"Failed to parse option '"<<argv[i]<<"'."<<std::endl;
                Print_Usage(true);}
            if(k<extras.m){
                if(extras(k).found) *extras(k).found=true;
                if(extras(k).store)
                    if(!extras(k).store_func(extras(k).store,argv[i])){
                        LOG::cerr<<"Failed to parse extra argument '"<<argv[i]<<"'."<<std::endl;
                        Print_Usage(true);}
                if(!extras(k).exhaust) k++;}
            else unclaimed_arguments=true;}
        if(k<extras.m && extras(k).exhaust) k++;
        for(;k<extras.m;k++)
            if(extras(k).required){
                LOG::cerr<<"Missing required extra argument '<"<<extras(k).name<<">'."<<std::endl;
                Print_Usage(true);}}
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
    LOG::cerr<<"\n"<<std::endl;

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
        LOG::cerr<<("<"+extras(i).name+">")<<extras(i).desc;
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
    s+="\nworking directory = "+Get_Working_Directory()+"\n";
    return s;
}
//#####################################################################
// Function Add
//#####################################################################
void PARSE_ARGS::
Add(const std::string& arg_str,bool* found,const std::string& desc)
{
    Add(arg_str,(int*)0,found,"",desc);
}
//#####################################################################
// Function Add_Not
//#####################################################################
void PARSE_ARGS::
Add_Not(const std::string& arg_str,bool* found,const std::string& desc)
{
    OPTION o={arg_str,"",desc,found,false,(int*)0,&store_impl<int>,&print_default_impl<int>};
    options.Set(arg_str,o);
}
//#####################################################################
// Destructor
//#####################################################################
PARSE_ARGS::OPTION::
~OPTION()
{
}
//#####################################################################
// Destructor
//#####################################################################
PARSE_ARGS::EXTRA::
~EXTRA()
{
}
//#####################################################################
namespace PhysBAM{
template void PARSE_ARGS::print_default_impl<int>(const void*);
template void PARSE_ARGS::print_default_impl<std::string>(const void*);
template void PARSE_ARGS::print_default_impl<float>(const void*);
template void PARSE_ARGS::print_default_impl<double>(const void*);
template bool PARSE_ARGS::store_impl<int>(void*, const std::string&);
template bool PARSE_ARGS::store_impl<std::string>(void*, const std::string&);
template bool PARSE_ARGS::store_impl<float>(void*, const std::string&);
template bool PARSE_ARGS::store_impl<double>(void*, const std::string&);
template void PARSE_ARGS::Add<int>(const std::string&, int*, bool*, const std::string&, const std::string&);
template void PARSE_ARGS::Add<std::string>(const std::string&, std::string*, bool*, const std::string&, const std::string&);
template void PARSE_ARGS::Add<float>(const std::string&, float*, bool*, const std::string&, const std::string&);
template void PARSE_ARGS::Add<double>(const std::string&, double*, bool*, const std::string&, const std::string&);

template ARRAY<HASHTABLE_ENTRY_TEMPLATE<std::string, PARSE_ARGS::OPTION>, int>::~ARRAY();
template bool HASHTABLE<std::string, PARSE_ARGS::OPTION>::Set(const std::string&, const PARSE_ARGS::OPTION&);
template void HASHTABLE<std::string, PARSE_ARGS::OPTION>::Resize_Table(int);
}
