//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FUNCTIONS.h>

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef VECTOR<T,2> TV;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_String_Argument("-in","input_directory");
    parse_args.Add_String_Argument("-data","density","Data to parse");
    parse_args.Add_Integer_Argument("-frame",-1);
    parse_args.Add_Integer_Argument("-xcut",-1);
    parse_args.Add_String_Argument("-surface","");

    parse_args.Add_Double_Argument("-xstart",0);
    parse_args.Add_Double_Argument("-transverse_velocity",-1);
    parse_args.Add_Option_Argument("-follow_transverse_velocity");

    parse_args.Parse();
    std::string input_directory=parse_args.Get_String_Value("-in"),
        data_file=parse_args.Get_String_Value("-data"),
        surface_file=parse_args.Get_String_Value("-surface");
    int frame=parse_args.Get_Integer_Value("-frame"),
        xcut=parse_args.Get_Integer_Value("-xcut");
    T xstart=(T)parse_args.Get_Double_Value("-xstart");

//##########################  INITIALIZATION  #########################
    GRID<TV> grid;
    ARRAY<T,VECTOR<int,2> > data;
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
        if(xcut >= 0){
            GRID<VECTOR<T,1> > cut_grid(grid.counts.Remove_Index(1),grid.domain.Remove_Dimension(1),grid.Is_MAC_Grid());
            for(GRID<VECTOR<T,1> >::CELL_ITERATOR iter(cut_grid);iter.Valid();iter.Next())
                LOG::cout<<iter.Location().x<<"\t"<<data(iter.Cell_Index().Insert(xcut,1))<<std::endl;}
        else if(parse_args.Get_Option_Value("-follow_transverse_velocity")){
            T time;FILE_UTILITIES::Read_From_File<RW>(f+"time",time);
            T transverse_velocity;FILE_UTILITIES::Read_From_File<RW>(input_directory+"/common/transverse_velocity",transverse_velocity);
            T x_position=xstart+time*transverse_velocity;
            LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
            GRID<VECTOR<T,1> > cut_grid(grid.counts.Remove_Index(1),grid.domain.Remove_Dimension(1),grid.Is_MAC_Grid());
            for(GRID<VECTOR<T,1> >::CELL_ITERATOR iter(cut_grid,2);iter.Valid();iter.Next())
                LOG::cout<<iter.Location().x<<"\t"<<interpolation.Clamped_To_Array(grid, data,iter.Location().Insert(x_position,1))<<std::endl;}
        else if(surface_file != ""){
            OCTAVE_OUTPUT<T> octave_output(surface_file.c_str());
            octave_output.Write(data_file.c_str(),data);}}

    return 0;
}
//#####################################################################
