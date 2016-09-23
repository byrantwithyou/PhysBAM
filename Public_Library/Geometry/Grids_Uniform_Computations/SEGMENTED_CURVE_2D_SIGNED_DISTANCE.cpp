//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Grids_Uniform_Computations/SEGMENTED_CURVE_2D_SIGNED_DISTANCE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
namespace PhysBAM{

namespace SIGNED_DISTANCE{
//#####################################################################
// Function Calculate
//#####################################################################
template<class T> void Calculate(SEGMENTED_CURVE_2D<T>& curve,const GRID<VECTOR<T,2> >& grid,ARRAY<T,VECTOR<int,2> >& phi,bool print_progress)
{
    bool bounding_box_defined=curve.bounding_box!=0;if(!bounding_box_defined) curve.Update_Bounding_Box();

    phi.Resize(0,grid.counts.x,0,grid.counts.y);

    T epsilon=(T)1e-8*grid.dX.Min();
    int total_cells=grid.counts.x*grid.counts.y,cells_done=0,progress=-1;
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
        phi(i,j)=curve.Calculate_Signed_Distance(grid.X(VECTOR<int,2>(i,j)),epsilon);
        if(print_progress){
            cells_done++;int new_progress=(int)((T)100*cells_done/total_cells);
            if(new_progress > progress){
                LOG::cout<<new_progress<<"% "<<std::flush;
                progress=new_progress;}}}
    if(print_progress) LOG::cout<<std::endl;

    if(!bounding_box_defined){delete curve.bounding_box;curve.bounding_box=0;}
}
//####################################################################
template void Calculate(SEGMENTED_CURVE_2D<float>&,const GRID<VECTOR<float,2> >&,ARRAY<float,VECTOR<int,2> >&,bool);
template void Calculate(SEGMENTED_CURVE_2D<double>&,const GRID<VECTOR<double,2> >&,ARRAY<double,VECTOR<int,2> >&,bool);
};
};
