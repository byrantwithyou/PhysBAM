//#####################################################################
// Copyright 2004-2005, Jiayi Chong
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Math_Tools/minmag.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Heat_Flows/HEAT_UNIFORM.h>
#include <boost/bind.hpp>
#include <PhysBAM_Geometry/Grids_RLE_Computations/DUALCONTOUR_RLE_3D.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>

using namespace PhysBAM;

//#####################################################################
// Function Rasterize_RLE
//#####################################################################
template<class T,class TV> void Rasterize_RLE(const LEVELSET_IMPLICIT_OBJECT<TV>& implicit_surface,RLE_GRID_3D<T>& grid,LEVELSET_RLE<RLE_GRID_3D<T> >& levelset)
{
    typedef RLE_GRID_3D<T> T_GRID;
    grid.Initialize(boost::bind(&IMPLICIT_OBJECT<TV>::Extended_Phi,boost::cref(implicit_surface),_1),0);
    LOG::cout<<"initial grid: cell = "<<grid.number_of_cells<<std::endl;
    ARRAY<T>& phi=levelset.phi;
    phi.Resize(grid.number_of_cells);

    LOG::Time("rasterizing");
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,T> interpolation;
    for(typename T_GRID::CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){int c=cell.Cell();
        TV X=cell.Center();
        phi(c)=implicit_surface(X);
        if(cell.Long()) phi(c+1)=phi(c);}
}

//#####################################################################
// Function Convert
//#####################################################################
template<class T,class RW> void Convert(PARSE_ARGS& parse_args,const std::string& input_file,const std::string& output_file,const std::string& output_grid)
{
    typedef VECTOR<T,3> TV;
    LEVELSET_IMPLICIT_OBJECT<TV> implicit_surface(*(new GRID<VECTOR<T,3> >),*(new ARRAYS<VECTOR<T,3> >));
    FILE_UTILITIES::Read_From_File<RW>(input_file,implicit_surface); 
    const GRID<VECTOR<T,3> >& uniform_grid=implicit_surface.levelset.grid;
  
    RLE_GRID_3D<T> rle_grid;
    int negative_bandwidth=parse_args.Get_Integer_Value("-negative_bandwidth");
    int positive_bandwidth=parse_args.Get_Integer_Value("-positive_bandwidth");
    if(uniform_grid.Is_MAC_Grid()) rle_grid.Set_Uniform_Grid(uniform_grid.Get_Regular_Grid());
    else rle_grid.Set_Uniform_Grid(uniform_grid.Get_MAC_Grid_At_Regular_Positions().Get_Regular_Grid());
    rle_grid.Set_Positive_Bandwidth_In_Cells(positive_bandwidth);
    rle_grid.Set_Negative_Bandwidth_In_Cells(negative_bandwidth);
    LOG::cout<<"rle positive bandwidth = "<<rle_grid.positive_bandwidth<<", negative_bandwidth = "<<rle_grid.negative_bandwidth<<std::endl;
    rle_grid.Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile();
    ARRAY<T> phi;LEVELSET_RLE<RLE_GRID_3D<T> > levelset(rle_grid,phi);
    Rasterize_RLE(implicit_surface,rle_grid,levelset);
    FILE_UTILITIES::Write_To_File<RW>(output_file,phi);
    FILE_UTILITIES::Write_To_File<RW>(output_grid,rle_grid);
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char *argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-i","","input levelset");
    parse_args.Add_String_Argument("-o","","output levelset");
    parse_args.Add_String_Argument("-g","","output grid");
    parse_args.Add_Integer_Argument("-negative_bandwidth",3,"rle levelset negative_bandwidth");
    parse_args.Add_Integer_Argument("-positive_bandwidth",3,"rle levelset positive_bandwidth");
    parse_args.Parse(argc,argv);

    std::string input_file,output_file,output_grid;
    input_file=parse_args.Get_String_Value("-i");
    output_file=parse_args.Get_String_Value("-o");
    output_grid=parse_args.Get_String_Value("-g");
    if(!parse_args.Is_Value_Set("-i")){LOG::cerr<<"Must specify -i"<<std::endl;exit(1);}
    if(!parse_args.Is_Value_Set("-o")){LOG::cerr<<"Must specify -o"<<std::endl;exit(1);}
    if(!parse_args.Is_Value_Set("-g")){LOG::cerr<<"Must specify -g"<<std::endl;exit(1);}

    std::cout<<"Reading in: "<<input_file<<std::endl;
    Convert<float,float>(parse_args,input_file,output_file,output_grid);
    std::cout<<"Wrote out levelset: "<<output_file<<std::endl;
    std::cout<<"Wrote out grid: "<<output_grid<<std::endl;

    return 0;
}
//#####################################################################
