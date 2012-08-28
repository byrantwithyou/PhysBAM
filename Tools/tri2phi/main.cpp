#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T> GRID<VECTOR<T,3> > Make_Cube_Grid(const GRID<VECTOR<T,3> >& input_grid,const int boundary_cells,const bool even_number_of_cells,const bool mac_grid)
{
    typedef VECTOR<T,3> TV;
    RANGE<TV> input_box=input_grid.Domain();

    // Want to make voxels into cubes
    T voxel_size=input_grid.Minimum_Edge_Length();
    int number_of_cells_m,number_of_cells_n,number_of_cells_mn;

    TV box_size=input_box.Edge_Lengths();
    TV box_center=input_box.Center();

    number_of_cells_m=(int)(box_size.x/voxel_size+(T).5)+2*boundary_cells;
    //if(fmod(box_size.x,voxel_size)!=0) adjusted_m+=2;
    if(even_number_of_cells && (number_of_cells_m&1)) number_of_cells_m++;

    number_of_cells_n=(int)(box_size.y/voxel_size+(T).5)+2*boundary_cells;
    //if(fmod(box_size.y,voxel_size)!=0) adjusted_n+=2;
    if(even_number_of_cells && (number_of_cells_n&1)) number_of_cells_n++;

    number_of_cells_mn=(int)(box_size.z/voxel_size+(T).5)+2*boundary_cells;
    //if(fmod(box_size.z,voxel_size)!=0) adjusted_mn+=2;
    if(even_number_of_cells && (number_of_cells_mn&1)) number_of_cells_mn++;

    TV new_size=voxel_size*TV((T)number_of_cells_m,(T)number_of_cells_n,(T)number_of_cells_mn);

    RANGE<TV> new_box;
    new_box.Reset_Bounds(box_center-(T).5*new_size);
    new_box.Enlarge_To_Include_Point(box_center+(T).5*new_size);

    if(mac_grid)
        return GRID<VECTOR<T,3> >(number_of_cells_m,number_of_cells_n,number_of_cells_mn,new_box,true);
    else
        return GRID<VECTOR<T,3> >(number_of_cells_m+1,number_of_cells_n+1,number_of_cells_mn+1,new_box,false);
}

template<class T> struct DEPTH_HELPER{ARRAY<PAIR<ORIENTED_BOX<VECTOR<T,3> >,int> > depths;int default_depth;};
template<class T>
static int Depth(int triangle_id,const RANGE<VECTOR<T,3> >& box,void* data)
{
    DEPTH_HELPER<T>& helper=*(DEPTH_HELPER<T>*)data;
    for(int i=0;i<helper.depths.m;i++) if(helper.depths(i).x.Intersection(box)) return helper.depths(i).y;
    return helper.default_depth;
}

template<class T,class RW> void Convert(const std::string& input_filename,int boundary_cells,PARSE_ARGS &parse_args)
{
    typedef VECTOR<T,3> TV;

    T dx=0,band=0,padding=0,thickness=(T)1e-5,offset=0;
    bool use_octree=false,mac=false,opt_heaviside(false),opt_unsigned(false),opt_only_bdy_outside(false);
    bool opt_keep_largest(false),opt_flip(false),opt_debug(false),opt_orthogonal_vote(false),use_grid_size(false);
    bool opt_domain_min(false),opt_domain_max(false),opt_path_start(false),opt_path_end(false);
    VECTOR<int,3> grid_size,path_start,path_end;
    int depth=1,positive_boundary_band=0;
    RANGE<TV> domain((TV()),TV());
    std::string output_filename,exact_grid;
    parse_args.Add("-heaviside",&opt_heaviside,"compute heaviside function");
    parse_args.Add("-unsigned",&opt_unsigned,"compute unsigned distance");
    parse_args.Add("-mac",&mac,"use MAC grid");
    parse_args.Add("-only_bdy_outside",&opt_only_bdy_outside,"only boundary region is outside");
    parse_args.Add("-keep_largest",&opt_keep_largest,"keep only largest inside region");
    parse_args.Add("-flip",&opt_flip,"flip sign if corners are inside");
    parse_args.Add("-debug",&opt_debug,"write debug data");
    parse_args.Add("-g",&grid_size,&use_grid_size,"grid size","suggested grid size - actual size may be larger");
    parse_args.Add("-dx",&dx,"grid dx","specify grid cell size");
    parse_args.Add("-domain_min",&domain.min_corner,&opt_domain_min,"vector","domain min");
    parse_args.Add("-domain_max",&domain.max_corner,&opt_domain_max,"vector","domain max");
    parse_args.Add("-band",&band,"half_band_width","half band width for fast marching method (in grid cells)");
    parse_args.Add("-padding",&padding,"padding","padding for flood fill");
    parse_args.Add("-thickness",&thickness,"thickness","surface thickness");
    parse_args.Add("-offset",&offset,"offset","subtract offset from phi");
    parse_args.Add("-positive_boundary",&positive_boundary_band,"positive boundary band","ensure this band around boundary is positive");
    parse_args.Add("-octree",&use_octree,"generate an octree levelset");
    parse_args.Add("-depth",&depth,"maximum_depth","maximum depth if using an octree");
    parse_args.Add("-o",&output_filename,"filename","output_filename");
    parse_args.Add("-exact_grid",&exact_grid,"exact_grid","exact_grid");
    parse_args.Add("-path_start",&path_start,&opt_path_start,"path_start","path_start");
    parse_args.Add("-path_end",&path_end,&opt_path_end,"path_end","path_end");
    parse_args.Add("-orthogonal_vote",&opt_orthogonal_vote,"orthogonal_vote");
    parse_args.Set_Extra_Arguments(1,"<tri file>","<tri file> tri file to convert");

    parse_args.Parse();

    TRIANGULATED_SURFACE<T>* triangulated_surface=0;FILE_UTILITIES::Create_From_File<RW>(input_filename,triangulated_surface);
    triangulated_surface->Update_Bounding_Box();
    RANGE<TV> box=*triangulated_surface->bounding_box;

    if(use_grid_size && dx){LOG::cerr<<"Only one of -g and -dx is allowed."<<std::endl;exit(1);}
    if(grid_size.Contains(0)){std::cerr<<"Invalid suggested grid size "<<grid_size<<std::endl;exit(1);}

    if(opt_domain_min && opt_domain_max) box=domain;

    // Make a cube grid using suggested grid_size and boundary cells
    GRID<VECTOR<T,3> > original_grid,grid;
    if(dx) grid=GRID<VECTOR<T,3> >::Create_Grid_Given_Cell_Size(box,dx,mac,boundary_cells);
    else{
        original_grid=GRID<VECTOR<T,3> >(grid_size,box,mac);
        if(mac && use_octree) PHYSBAM_FATAL_ERROR("Can't use MAC grids on an octree");
        grid=Make_Cube_Grid(original_grid,boundary_cells,use_octree,mac);}

    if(!exact_grid.empty()){
        FILE_UTILITIES::Read_From_File<RW>(exact_grid,grid);
        LOG::cout<<"reading grid from "<<exact_grid<<std::endl;}

    if(output_filename.empty()){
        std::string dimensions=STRING_UTILITIES::string_sprintf("%dx%dx%d",grid.counts.x,grid.counts.y,grid.counts.z);
        output_filename=FILE_UTILITIES::Get_Basename(input_filename)+dimensions+(use_octree?".oct":".phi");}

    std::cout<<"Input filename: "<<input_filename<<std::endl;
    std::cout<<"Output filename: "<<output_filename<<std::endl;
    std::cout<<"------------------------------------"<<std::endl;
    std::cout<<"Suggested grid size = "<<grid_size<<std::endl;
    std::cout<<"Use octree: "<<use_octree<<", octree depth: "<<depth<<std::endl;
    std::cout<<"Adjusted to make cube voxels and added boundary cells ("<<boundary_cells<<")..."<<std::endl;
    std::cout<<"New number of nodes: m = "<<grid.counts.x<<", n = "<<grid.counts.y<<", mn = "<<grid.counts.z<<std::endl;
    std::cout<<"Voxel size: "<<grid.Minimum_Edge_Length()<<std::endl;
    if(opt_heaviside)
        std::cout<<"Compute heaviside function"<<std::endl;
    else
        std::cout<<"Fast marching method half band width (in grid cells): "<<band<<std::endl;
    std::cout<<"Padding for flood fill: "<<padding<<std::endl;
    std::cout<<"Positive boundary band: "<<positive_boundary_band<<std::endl;
    std::cout<<"------------------------------------"<<std::endl;
    std::cout<<"Original box: -domain_min "<<box.min_corner<<" -domain_max "<<box.max_corner<<std::endl;
    std::cout<<"New box: -domain_min "<<grid.Domain().min_corner<<" -domain_max "<<grid.Domain().max_corner<<std::endl;
    std::cout<<"------------------------------------"<<std::endl;
    std::cout<<std::endl;


    {
        LEVELSET_MAKER_UNIFORM<T> levelset_maker;
        levelset_maker.Verbose_Mode(true);
        levelset_maker.Write_Debug_Data(opt_debug);
        if(opt_path_start && opt_path_end)
            levelset_maker.Write_Debug_Path(true,path_start,path_end);
        if(opt_orthogonal_vote) levelset_maker.Use_Orthogonal_Vote(true);
        levelset_maker.Set_Surface_Padding_For_Flood_Fill(padding);
        levelset_maker.Set_Surface_Thickness(thickness);
        if(opt_heaviside)
            levelset_maker.Compute_Heaviside_Function();
        else{
            if(opt_unsigned)
                levelset_maker.Compute_Unsigned_Distance_Function();
            else
                levelset_maker.Compute_Signed_Distance_Function();
            levelset_maker.Use_Fast_Marching_Method(true,band);}
            levelset_maker.Only_Boundary_Region_Is_Outside(opt_only_bdy_outside);
            levelset_maker.Keep_Only_Largest_Inside_Region(opt_keep_largest);
            levelset_maker.Flip_Sign_If_Corners_Are_Inside(opt_flip);
            levelset_maker.Set_Positive_Boundary_Band(positive_boundary_band);
            levelset_maker.Set_Phi_Offset(offset);
            ARRAY<T,VECTOR<int,3> > phi(grid.Domain_Indices());
            levelset_maker.Compute_Level_Set(*triangulated_surface,grid,phi);
            LEVELSET_IMPLICIT_OBJECT<TV> levelset_implicit_surface(grid,phi);
            //phi+=(T)1*grid.Maximum_Edge_Length();
        FILE_UTILITIES::Write_To_File<RW>(output_filename,levelset_implicit_surface);}
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false,type_float=false,compute_using_doubles=false;
    VECTOR<double,3> grid_size(50,50,50);
    int boundary_cells=3;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-float",&type_float,"data is in float format");
    parse_args.Add("-double",&type_double,"data is in double format");
    parse_args.Add("-b",&boundary_cells,"boundary cells","number of cells outside bounding box");
    parse_args.Add("-compute_using_doubles",&compute_using_doubles,"perform computations using doubles");
    parse_args.Parse(true);

    std::string input_filename=parse_args.Extra_Arg(0);

    if(!FILE_UTILITIES::Is_Tri_File(input_filename)){
        std::cerr<<"Not a tri file: "<<input_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(compute_using_doubles){
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
            std::cout<<"COMPUTING USING DOUBLES!"<<std::endl;
            Convert<double,float>(input_filename,boundary_cells,parse_args);
#else
            std::cerr<<"Double support not enabled."<<std::endl;exit(1);
#endif
        }else{Convert<float,float>(input_filename,boundary_cells,parse_args);}}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double,double>(input_filename,boundary_cells,parse_args);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
