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

    std::string input_directory="input_directory",data_file="density",surface_file;
    int frame=-1,xcut=-1;
    T xstart=0;
    bool opt_follow_transverse_velocity=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-in",&input_directory,"dir","input_directory");
    parse_args.Add("-data",&data_file,"file","Data to parse");
    parse_args.Add("-frame",&frame,"frame","frame");
    parse_args.Add("-xcut",&xcut,"xcut","xcut");
    parse_args.Add("-surface",&surface_file,"file","surface file");
    parse_args.Add("-follow_transverse_velocity",&opt_follow_transverse_velocity,"","");
    parse_args.Parse();

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
            for(UNIFORM_GRID_ITERATOR_CELL<VECTOR<T,1> > iter(cut_grid);iter.Valid();iter.Next())
                LOG::cout<<iter.Location().x<<"\t"<<data(iter.Cell_Index().Insert(xcut,1))<<std::endl;}
        else if(opt_follow_transverse_velocity){
            T time;FILE_UTILITIES::Read_From_File<RW>(f+"time",time);
            T transverse_velocity;FILE_UTILITIES::Read_From_File<RW>(input_directory+"/common/transverse_velocity",transverse_velocity);
            T x_position=xstart+time*transverse_velocity;
            LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
            GRID<VECTOR<T,1> > cut_grid(grid.counts.Remove_Index(1),grid.domain.Remove_Dimension(1),grid.Is_MAC_Grid());
            for(UNIFORM_GRID_ITERATOR_CELL<VECTOR<T,1> > iter(cut_grid,2);iter.Valid();iter.Next())
                LOG::cout<<iter.Location().x<<"\t"<<interpolation.Clamped_To_Array(grid, data,iter.Location().Insert(x_position,1))<<std::endl;}
        else if(surface_file != ""){
            OCTAVE_OUTPUT<T> octave_output(surface_file.c_str());
            octave_output.Write(data_file.c_str(),data);}}

    return 0;
}
//#####################################################################
