//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Parsing/STRING_UTILITIES.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Read_Write/READ_WRITE_FUNCTIONS.h>

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef VECTOR<T,1> TV;
    typedef float RW;
    STREAM_TYPE stream_type(RW());

    std::string input_directory="input_directory",data_file="density";
    int frame=-1;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-in",&input_directory,"dir","input directory");
    parse_args.Add("-data",&data_file,"file","Data to parse");
    parse_args.Add("-frame",&frame,"frame","frame");
    parse_args.Parse();

//##########################  INITIALIZATION  #########################
    GRID<TV> grid;
    ARRAY<T,VECTOR<int,1> > data;
    FILE_UTILITIES::Read_From_File<RW>(input_directory+"/common/grid",grid);
    std::string f=STRING_UTILITIES::string_sprintf("/%d/",frame);
    if(frame >= 0){
        FILE_UTILITIES::Read_From_File<RW>(input_directory+f+"/"+data_file,data);}
//#####################################################################

    LOG::cerr<<"Grid = "<<grid<<std::endl;
    LOG::cerr<<"Grid.DX = "<<grid.DX()<<std::endl;

    if(frame >= 0){
        LOG::cerr<<"Array = "<<data.Domain_Indices()<<std::endl;
        LOG::cerr<<"min = "<<data.Min()<<", max = "<<data.Max()<<std::endl;

        LOG::cout.precision(std::numeric_limits< double >::digits10);
        for(CELL_ITERATOR<VECTOR<T,1> > iter(grid,3);iter.Valid();iter.Next())
            LOG::cout<<iter.Location().x<<"\t"<<data(iter.Cell_Index())<<std::endl;
    }

    return 0;
}
//#####################################################################
