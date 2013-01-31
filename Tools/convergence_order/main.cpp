//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
using namespace PhysBAM;
//#####################################################################
// Function Write_Output
//#####################################################################
template<class TV> void
Write_Output(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;

    int frame=0,frame_rate=24;
    T start=0,end=0;
    std::string input_directory1,input_directory2;
    parse_args.Add("-start",&start,"value","start range");
    parse_args.Add("-end",&end,"value","end range");
    parse_args.Add("-frame",&frame,"value","frame output");
    parse_args.Add("-frame_rate",&frame_rate,"value","frame rate used");
    parse_args.Extra(&input_directory1,"low res input_directory","low res input_directory");
    parse_args.Extra(&input_directory2,"high res input_directory","high res input_directory");
    parse_args.Parse();

    if(frame==0) FILE_UTILITIES::Read_From_Text_File(input_directory1+"/common/last_frame",frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    GRID<TV> coarse_grid;FILE_UTILITIES::Read_From_File<T>(input_directory1+"/common/grid",coarse_grid);
    GRID<TV> grid;FILE_UTILITIES::Read_From_File<T>(input_directory2+"/common/grid",grid);
    int scale=grid.Counts().x/coarse_grid.Counts().x;
    ARRAY<T,TV_INT> density_coarse;FILE_UTILITIES::Read_From_File<T>(input_directory1+"/"+f+"/density",density_coarse);
    ARRAY<T,TV_INT> density;FILE_UTILITIES::Read_From_File<T>(input_directory2+"/"+f+"/density",density);
    RANGE<TV> range(TV()+start,TV()+end);
    RANGE<TV_INT> cell_range;cell_range.min_corner=coarse_grid.Index(range.min_corner);cell_range.max_corner=coarse_grid.Index(range.max_corner);
    if(end==0) cell_range=grid.Domain_Indices();
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(coarse_grid,cell_range);iterator.Valid();iterator.Next()){
        TV X_o=iterator.Location()-((T)frame)*1/(T)frame_rate;
        TV X_o_fine=grid.Center(scale*iterator.Cell_Index()-(scale-1)*TV_INT::All_Ones_Vector())-((T)frame)*1/(T)frame_rate;
        T analytic=0;if(X_o.x>=0.25&&X_o.x<=0.75) analytic=.5*(sin(2*pi/(.5)*(X_o.x-0.25)-pi/2.)+1);
//        T analytic_fine=0;if(X_o_fine.x>=0.25&&X_o_fine.x<=0.75) analytic_fine=.5*(sin(2*pi/(.5)*(X_o_fine.x-0.25)-pi/2.)+1);
        T coarse_error=analytic-(density_coarse(iterator.Cell_Index()));
        //T fine_error=analytic_fine-(density(scale*iterator.Cell_Index()-(scale-1)*TV_INT::All_Ones_Vector()));
        //PHYSBAM_ASSERT(scale*iterator.Cell_Index()-(scale-1)*TV_INT::All_Ones_Vector()==grid.Index(iterator.Location()));
        //T r=abs(fine_error)>1e-10?(coarse_error/fine_error):0;
        //T r=abs(fine_error)>1e-12?log(abs(coarse_error/fine_error))/log(2):0;
        //if(abs(r)>4) r=0;
        //LOG::cout<<iterator.Location().x<<" "<<r<<std::endl;}
        if(range.Lazy_Inside(iterator.Location())) LOG::cout<<iterator.Location().x<<" "<<coarse_error<<std::endl;}
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

    if(opt_3d){
        Write_Output<VECTOR<T,3> >(parse_args);}
    if(opt_2d){
        Write_Output<VECTOR<T,2> >(parse_args);}
    else{
        Write_Output<VECTOR<T,1> >(parse_args);}
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char* argv[])
{
    bool opt_double=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-double",&opt_double,"input data is in doubles");
    parse_args.Parse(true);

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(opt_double) Find_Dimension<double>(parse_args); 
    else Find_Dimension<float>(parse_args);
#else
    Find_Dimension<float>(parse_args);
#endif
}
//#####################################################################
