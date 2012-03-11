//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
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
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        total_density+=density(iterator.Cell_Index())*grid.min_dX;}
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
    for(typename GRID<TV>::FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
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
    for(typename GRID<TV>::FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        if(kinetic) total_k+=u(iterator.Full_Index())*grid.dX(iterator.Axis());
        else total_k+=.5*1000*u(iterator.Full_Index())*u(iterator.Full_Index())*grid.dX(iterator.Axis());}
    LOG::cout<<total_k<<std::endl;
}
//#####################################################################
// Function Print_Density
//#####################################################################
template<class TV> void
Print_Density(std::string& input_directory,GRID<TV>& grid,int frame,RANGE<TV>& range)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;
    
    if(frame==-1){int last_frame=0;FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);for(int i=0;i<=last_frame;i++) Print_Density(input_directory,grid,i,range);return;}
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    ARRAY<T,TV_INT> density;FILE_UTILITIES::Read_From_File<T>(input_directory+"/"+f+"/density",density);
    RANGE<TV_INT> cell_range;cell_range.min_corner=grid.Index(range.min_corner);cell_range.max_corner=grid.Index(range.max_corner);
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid,cell_range);iterator.Valid();iterator.Next()){
        LOG::cout<<iterator.Location().x<<"\t"<<density(iterator.Cell_Index())<<std::endl;}
}
//#####################################################################
// Function Write_Output
//#####################################################################
template<class TV> void
Write_Output(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;

    std::string input_directory=parse_args.Extra_Arg(1);
    int frame=parse_args.Get_Integer_Value("-frame");
    GRID<TV> grid;FILE_UTILITIES::Read_From_File<T>(input_directory+"/common/grid",grid);
    if(parse_args.Is_Value_Set("-m")) Print_Mass(input_directory,grid,frame);
    else if(parse_args.Is_Value_Set("-l")) Print_Momentum(input_directory,grid,frame);
    else if(parse_args.Is_Value_Set("-e")) Print_Energy(input_directory,grid,frame);
    else{
        RANGE<TV> range;range.min_corner=TV::All_Ones_Vector()*parse_args.Get_Double_Value("-start");
        range.max_corner=TV::All_Ones_Vector()*parse_args.Get_Double_Value("-end");Print_Density(input_directory,grid,frame,range);}
}
//#####################################################################
// Function Find_Dimension
//#####################################################################
template<class T> void
Find_Dimension(PARSE_ARGS& parse_args)
{
    if(parse_args.Is_Value_Set("-3d")){
        Write_Output<VECTOR<T,3> >(parse_args);}
    if(parse_args.Is_Value_Set("-2d")){
        Write_Output<VECTOR<T,2> >(parse_args);}
    else{
        Write_Output<VECTOR<T,1> >(parse_args);}
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char* argv[])
{
    Initialize_Particles();
    Initialize_Read_Write_Structures();

    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-start",0,"start range");
    parse_args.Add_Double_Argument("-end",0,"end range");
    parse_args.Add_Integer_Argument("-frame",-1,"frame output");
    parse_args.Add_Option_Argument("-m","print mass");
    parse_args.Add_Option_Argument("-l","print monmentum");
    parse_args.Add_Option_Argument("-e","print energy");
    parse_args.Add_Option_Argument("-double","input data is in doubles");
    parse_args.Add_Option_Argument("-2d","input data is 2-D");
    parse_args.Add_Option_Argument("-3d","input data is 3-D");
    parse_args.Set_Extra_Arguments(1,"<input_directory>");
    parse_args.Parse(argc,argv);

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(parse_args.Is_Value_Set("-double")) Find_Dimension<double>(parse_args); 
    else Find_Dimension<float>(parse_args);
#else
    Find_Dimension<float>(parse_args);
#endif
}
//#####################################################################
