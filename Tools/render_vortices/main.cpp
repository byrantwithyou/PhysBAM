#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;
template<class T>
void Run(const STREAM_TYPE& stream_type,const std::string& file)
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    int ghost=6;
    GRID<TV> grid;
    FILE_UTILITIES::Read_From_File(stream_type,file+"/grid",grid);
    ARRAY<T,FACE_INDEX<2> > uo;
    FILE_UTILITIES::Read_From_File(stream_type,file+"/mac_velocities",uo);
    ARRAY<T,FACE_INDEX<2> > u(grid,ghost);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()) u(it.Full_Index())=uo(it.Full_Index());

    ARRAY<T,TV_INT> phi_array(grid.Domain_Indices(ghost));

    LEVELSET_2D<GRID<VECTOR<T,2> > > phi(grid,phi_array);

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,ghost);it.Valid();it.Next()){
        T val=-FLT_MAX;
        TV X=it.Location();
        val=std::max(val,(T).5-X.Magnitude());
        val=std::max(val,-X.x-(T)2.5);
        val=std::max(val,X.x-15);
        val=std::max(val,X.y-(T)3.5);
        val=std::max(val,-(T)3.5-X.y);
        phi_array(it.index)=val;}

    ARRAY<bool,FACE_INDEX<2> > inside(grid);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next())
        inside(it.Full_Index())=phi.Phi(it.Location());
    EXTRAPOLATION_HIGHER_ORDER<TV,T>::Extrapolate_Face(grid,phi,inside,ghost,u,100,3,ghost-1);

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,ghost);it.Valid();it.Next()){
        T du=u(FACE_INDEX<2>(1,it.index+TV_INT(1,0)))-u(FACE_INDEX<2>(1,it.index));
        T dv=u(FACE_INDEX<2>(2,it.index+TV_INT(0,1)))-u(FACE_INDEX<2>(2,it.index));
        T vor=du/grid.one_over_dX.x-dv/grid.one_over_dX.y;
        LOG::cout<<(it.index.x-1+ghost)<<" "<<(it.index.y-1+ghost)<<" "<<vor<<std::endl;
        if(it.index.y==grid.counts.y+ghost) LOG::cout<<std::endl;}
    
}

int main(int argc,char *argv[])
{
    //PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false;

    std::string filename;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Extra(&filename,"directory","directory to render vortices for");
    parse_args.Parse();

    STREAM_TYPE stream_type(type_double?STREAM_TYPE(0.0):STREAM_TYPE(0.0f));

    if(!type_double) Run<float>(stream_type,filename);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(type_double) Run<double>(stream_type,filename);
#endif
}
