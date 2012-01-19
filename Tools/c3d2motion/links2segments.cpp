//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Dynamics/Motion/MOTION_SEQUENCE.h>
#include "C3D_FILE.h"

using namespace PhysBAM;
//#####################################################################
// Function Split
//#####################################################################
void 
Split(const std::string& text,const std::string& delimiter,ARRAY<std::string>& output)
{
    unsigned int offset=0,delimiter_index=0;
    delimiter_index=text.find(delimiter,offset);
    while(delimiter_index!=std::string::npos){
        output.Append(text.substr(offset,delimiter_index-offset));
        offset+=delimiter_index-offset+delimiter.length();
        delimiter_index=text.find(delimiter,offset);}
    output.Append(text.substr(offset));
}
//#####################################################################
// Function Convert
//#####################################################################
template<class T> bool
Convert(const std::string& input_filename,const std::string& output_filename)
{
    ARRAY<VECTOR<int,2> > segments(2,false);
    std::istream* istream=FILE_UTILITIES::Safe_Open_Input(input_filename);
    std::istream& infile=*istream;
    std::string line;bool found_linkages=false;
    while(getline(infile,line,'\n')){
        if(line=="[Linkages]"){found_linkages=true;break;}}
    if(!found_linkages){std::cerr<<"Failed to find [Linkages]"<<std::endl;delete istream;return false;}
    getline(infile,line,'\n'); // skip comment
    // get linkage count
    getline(infile,line,'\n');
    std::string::size_type start=line.rfind("=");
    if(start==std::string::npos){std::cerr<<"Failed to get the linkage count"<<std::endl;delete istream;return false;}
    start+=1;
    int linkages=atoi(line.substr(start,line.length()-start).c_str());
    // get linkages
    for(int i=0;i<linkages;i++){
        getline(infile,line,'\n');
        ARRAY<std::string> items,numbers;
        Split(line,",",items);
        if(items.m<4){std::cerr<<"Expected more items on linkage line"<<std::endl;delete istream;return false;}
        Split(items(4)," ",numbers);
        if(numbers.m!=3){std::cerr<<"Expected 2 numbers on linkage line"<<std::endl;delete istream;return false;}
        int i=atoi(numbers(2).c_str()),j=atoi(numbers(3).c_str());
        segments.Append(VECTOR<int,2>(i,j));}
    FILE_UTILITIES::Write_To_File<T>(output_filename,segments);
    return true;
}
//#####################################################################
// Function main
//#####################################################################
int 
main(int argc,char *argv[])
{
    if(argc != 3){std::cerr<<"Usage is: "<<argv[0]<<" <input .prj> <output motion>"<<std::endl;return 1;}
    Convert<float>(argv[1],argv[2]);
    return 0;
}
//#####################################################################
