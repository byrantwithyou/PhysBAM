//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef VECTOR<T,1> TV;
    typedef float RW;
    STREAM_TYPE stream_type(RW());

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-in","input_directory");
    parse_args.Add_String_Argument("-data","density","Data to parse");
    parse_args.Add_Integer_Argument("-frame",-1);
    parse_args.Parse(argc,argv);
    std::string input_directory=parse_args.Get_String_Value("-in"),
        data_file=parse_args.Get_String_Value("-data");
    int frame=parse_args.Get_Integer_Value("-frame");

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
        for(GRID<VECTOR<T,1> >::CELL_ITERATOR iter(grid,3);iter.Valid();iter.Next())
            LOG::cout<<iter.Location().x<<"\t"<<data(iter.Cell_Index())<<std::endl;
    }

    return 0;
}
//#####################################################################
