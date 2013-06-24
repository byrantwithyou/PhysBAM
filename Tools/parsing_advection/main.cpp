//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
using namespace PhysBAM;
//#####################################################################
// Function Print_Mass
//#####################################################################
template<class TV> void
Print_Mass(std::string& input_directory,GRID<TV>& grid,int frame)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;
    
    if(frame==-1){int last_frame=0;FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);for(int i=0;i<last_frame;i++) Print_Mass(input_directory,grid,i);return;}
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    T total_density=0;ARRAY<T,TV_INT> density;FILE_UTILITIES::Read_From_File<T>(input_directory+"/"+f+"/density",density);
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        total_density+=density(iterator.Cell_Index())*grid.dX.Min();}
    LOG::cout<<total_density<<std::endl;
}
//#####################################################################
// Function Print_Momentum
//#####################################################################
template<class TV> void
Print_Momentum(std::string& input_directory,GRID<TV>& grid,int frame)
{
    typedef typename TV::SCALAR T;
    
    if(frame==-1){int last_frame=0;FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);for(int i=0;i<last_frame;i++) Print_Momentum(input_directory,grid,i);return;}
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    T total_l=0;ARRAY<T,FACE_INDEX<TV::dimension> > u;FILE_UTILITIES::Read_From_File<T>(input_directory+"/"+f+"/mac_velocities",u);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        total_l+=u(iterator.Full_Index())*grid.dX(iterator.Axis());}
    LOG::cout<<total_l<<std::endl;
}
//#####################################################################
// Function Print_Energy
//#####################################################################
template<class TV> void
Print_Energy(std::string& input_directory,GRID<TV>& grid,int frame)
{
    typedef typename TV::SCALAR T;
    
    if(frame==-1){int last_frame=0;FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);for(int i=0;i<last_frame;i++) Print_Energy(input_directory,grid,i);return;}
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    std::string filename=input_directory+"/%d/kinetic_energy";
    T total_k=0;ARRAY<T,FACE_INDEX<TV::dimension> > u;bool kinetic=false;
    if(FILE_UTILITIES::Frame_File_Exists(filename,frame)){kinetic=true;FILE_UTILITIES::Read_From_File<T>(input_directory+"/"+f+"/kinetic_energy",u);}
    else FILE_UTILITIES::Read_From_File<T>(input_directory+"/"+f+"/mac_velocities",u);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        if(kinetic) total_k+=u(iterator.Full_Index())*grid.dX(iterator.Axis());
        else total_k+=.5*1000*u(iterator.Full_Index())*u(iterator.Full_Index())*grid.dX(iterator.Axis());}
    LOG::cout<<total_k<<std::endl;
}
//#####################################################################
// Function Print_Density
//#####################################################################
template<class TV> void
Print_Density(std::string& input_directory,GRID<TV>& grid,int frame,const RANGE<TV>& range)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;
    
    if(frame==-1){int last_frame=0;FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);for(int i=0;i<=last_frame;i++) Print_Density(input_directory,grid,i,range);return;}
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    ARRAY<T,TV_INT> density;FILE_UTILITIES::Read_From_File<T>(input_directory+"/"+f+"/density",density);
    RANGE<TV_INT> cell_range;cell_range.min_corner=grid.Index(range.min_corner);cell_range.max_corner=grid.Index(range.max_corner);
    for(CELL_ITERATOR<TV> iterator(grid,cell_range);iterator.Valid();iterator.Next()){
        LOG::cout<<iterator.Location().x<<"\t"<<density(iterator.Cell_Index())<<std::endl;}
}
//#####################################################################
// Function Write_Output
//#####################################################################
template<class TV> void
Write_Output(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;

    std::string input_directory;
    bool opt_m=false,opt_l=false,opt_e=false;
    int frame=-1;
    T start=0,end=0;
    parse_args.Add("-start",&start,"value","start range");
    parse_args.Add("-end",&end,"value","end range");
    parse_args.Add("-frame",&frame,"value","frame output");
    parse_args.Add("-m",&opt_m,"print mass");
    parse_args.Add("-l",&opt_l,"print monmentum");
    parse_args.Add("-e",&opt_e,"print energy");
    parse_args.Extra(&input_directory,"input_directory","input_directory");
    parse_args.Parse();

    GRID<TV> grid;FILE_UTILITIES::Read_From_File<T>(input_directory+"/common/grid",grid);
    if(opt_m) Print_Mass(input_directory,grid,frame);
    else if(opt_l) Print_Momentum(input_directory,grid,frame);
    else if(opt_e) Print_Energy(input_directory,grid,frame);
    else Print_Density(input_directory,grid,frame,RANGE<TV>(TV()+start,TV()+end));
}
//#####################################################################
// Function Find_Dimension
//#####################################################################
template<class T> void
Find_Dimension(PARSE_ARGS& parse_args)
{
    bool opt_2d=false,opt_3d=false;
    parse_args.Add("-2d",&opt_2d,"input data is 2-D");
    parse_args.Add("-3d",&opt_3d,"input data is 3-D");
    parse_args.Parse(true);

    if(opt_3d)
        Write_Output<VECTOR<T,3> >(parse_args);
    else if(opt_2d)
        Write_Output<VECTOR<T,2> >(parse_args);
    else
        Write_Output<VECTOR<T,1> >(parse_args);
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char* argv[])
{
    bool use_doubles=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-double",&use_doubles,"input data is in doubles");
    parse_args.Parse(true);

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(use_doubles) Find_Dimension<double>(parse_args); 
    else Find_Dimension<float>(parse_args);
#else
    Find_Dimension<float>(parse_args);
#endif
}
//#####################################################################
