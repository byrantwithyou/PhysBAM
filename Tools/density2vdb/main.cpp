//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <openvdb/openvdb.h>
using namespace PhysBAM;

template<class T>
void WriteVDB(const std::string& output_filename,const ARRAY<T,VECTOR<int,3> >& density)
{
    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid=openvdb::FloatGrid::create(); 
    grid->setGridClass(openvdb::GRID_FOG_VOLUME);
    openvdb::FloatGrid::Accessor accessor=grid->getAccessor();
    openvdb::Coord xyz(0,0,0);
    VECTOR<int,3> start=density.domain.min_corner;
    VECTOR<int,3> counts(density.domain.Edge_Lengths());
    for(int i=0;i<counts.x;i++){for(int j=0;j<counts.y;j++){for(int k=0;k<counts.z;k++){
                xyz.reset(i,j,k);
                accessor.setValue(xyz,density(VECTOR<int,3>(start.x+i,start.y+j,start.z+k)));}}}
    openvdb::io::File file(output_filename);
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    file.write(grids);
    file.close();
}
template<class T>
void Convert(const std::string& input_directory,const std::string& output_filename_pattern)
{
    ARRAY<T,VECTOR<int,3> > density; 
    int last_frame;
    Read_From_Text_File(input_directory+"/common/last_frame",last_frame);
    for(int i=0;i<=last_frame;++i){
        LOG::cout<<"Reading frame "<<i<<" density."<<std::endl;
        Read_From_File<T>(LOG::sprintf("%s/%d/density.gz",input_directory,i),density);
        WriteVDB<T>(LOG::sprintf(output_filename_pattern.c_str(),i),density);}
}
int main(int argc,char *argv[])
{
    bool type_double=false;
    std::string input_directory,output_filename_pattern;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-double",&type_double,"Load doubles (for the density file)");
    parse_args.Extra(&input_directory,"sim-dir","simulation output directory");
    parse_args.Extra(&output_filename_pattern,"output-pattern","output vdb file name pattern, use %d");
    parse_args.Parse();
    if(type_double){
        LOG::cout<<"Loading using double.. (But still write as floats.)"<<std::endl;
        Convert<double>(input_directory,output_filename_pattern);}
    else Convert<float>(input_directory,output_filename_pattern);
    return 0;
}
