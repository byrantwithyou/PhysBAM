#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
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
    for(int i=1;i<=helper.depths.m;i++) if(helper.depths(i).x.Intersection(box)) return helper.depths(i).y;
    return helper.default_depth;
}

template<class T,class RW> void Convert(const std::string& input_filename,int boundary_cells,PARSE_ARGS &parse_args)
{
    typedef VECTOR<T,3> TV;

    TRIANGULATED_SURFACE<T>* triangulated_surface=0;FILE_UTILITIES::Create_From_File<RW>(input_filename,triangulated_surface);
    triangulated_surface->Update_Bounding_Box();
    RANGE<TV> box=*triangulated_surface->bounding_box;

    T dx=(T)parse_args.Get_Double_Value("-dx");
    if(parse_args.Is_Value_Set("-g") && dx){LOG::cerr<<"Only one of -g and -dx is allowed."<<std::endl;exit(1);}
    VECTOR<int,3> grid_size=VECTOR<int,3>(parse_args.Get_Vector_3D_Value("-g"));
    if(grid_size.Contains(0)){std::cerr<<"Invalid suggested grid size "<<grid_size<<std::endl;exit(1);}
    bool use_octree=parse_args.Get_Option_Value("-octree");
    int depth=parse_args.Get_Integer_Value("-depth");
    bool mac=parse_args.Get_Option_Value("-mac");

    if(parse_args.Is_Value_Set("-domain_min") && parse_args.Is_Value_Set("-domain_max")){
        box=RANGE<TV>((TV)parse_args.Get_Vector_3D_Value("-domain_min"),(TV)parse_args.Get_Vector_3D_Value("-domain_max"));}

    // Make a cube grid using suggested grid_size and boundary cells
    GRID<VECTOR<T,3> > original_grid,grid;
    if(dx) grid=GRID<VECTOR<T,3> >::Create_Grid_Given_Cell_Size(box,dx,mac,boundary_cells);
    else{
        original_grid=GRID<VECTOR<T,3> >(grid_size,box,mac);
        if(mac && use_octree) PHYSBAM_FATAL_ERROR("Can't use MAC grids on an octree");
        grid=Make_Cube_Grid(original_grid,boundary_cells,use_octree,mac);}

    if(parse_args.Is_Value_Set("-exact_grid")){
        std::string filename=parse_args.Get_String_Value("-exact_grid");
        FILE_UTILITIES::Read_From_File<RW>(filename,grid);
        LOG::cout<<"reading grid from "<<filename<<std::endl;}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    OCTREE_GRID<T> octree_grid;DEPTH_HELPER<T> helper;
    if(use_octree){
        helper.default_depth=depth;
        //helper.depths.Append(PAIR<ORIENTED_BOX_3D<T>,int>(ORIENTED_BOX_3D<T>(TV(-.039,.089,-.07),TV(.076,0,0),TV(0,.04,0),TV(0,0,.026)),depth));
        octree_grid.Initialize(grid,depth,3,true,true);}
#endif

    std::string output_filename="";
    if(parse_args.Is_Value_Set("-o")) output_filename=parse_args.Get_String_Value("-o");
    else{
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
    if(parse_args.Get_Option_Value("-heaviside"))
        std::cout<<"Compute heaviside function"<<std::endl;
    else
        std::cout<<"Fast marching method half band width (in grid cells): "<<parse_args.Get_Double_Value("-band")<<std::endl;
    std::cout<<"Padding for flood fill: "<<parse_args.Get_Double_Value("-padding")<<std::endl;
    std::cout<<"Positive boundary band: "<<parse_args.Get_Integer_Value("-positive_boundary")<<std::endl;
    std::cout<<"------------------------------------"<<std::endl;
    std::cout<<"Original box: -domain_min "<<box.min_corner<<" -domain_max "<<box.max_corner<<std::endl;
    std::cout<<"New box: -domain_min "<<grid.Domain().min_corner<<" -domain_max "<<grid.Domain().max_corner<<std::endl;
    std::cout<<"------------------------------------"<<std::endl;
    std::cout<<std::endl;


#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(!use_octree)
#endif
    {
        LEVELSET_MAKER_UNIFORM<T> levelset_maker;
        levelset_maker.Verbose_Mode(true);
        levelset_maker.Write_Debug_Data(parse_args.Get_Option_Value("-debug"));
        if(parse_args.Is_Value_Set("-path_start")&&parse_args.Is_Value_Set("-path_end"))
            levelset_maker.Write_Debug_Path(true,VECTOR<int,3>(parse_args.Get_Vector_3D_Value("-path_start")),VECTOR<int,3>(parse_args.Get_Vector_3D_Value("-path_end")));
        if(parse_args.Get_Option_Value("-orthogonal_vote")) levelset_maker.Use_Orthogonal_Vote(true);
        levelset_maker.Set_Surface_Padding_For_Flood_Fill((T)parse_args.Get_Double_Value("-padding"));
        levelset_maker.Set_Surface_Thickness((T)parse_args.Get_Double_Value("-thickness"));
        if(parse_args.Get_Option_Value("-heaviside"))
            levelset_maker.Compute_Heaviside_Function();
        else{
            if(parse_args.Get_Option_Value("-unsigned"))
                levelset_maker.Compute_Unsigned_Distance_Function();
            else
                levelset_maker.Compute_Signed_Distance_Function();
            levelset_maker.Use_Fast_Marching_Method(true,(T)parse_args.Get_Double_Value("-band"));}
            levelset_maker.Only_Boundary_Region_Is_Outside(parse_args.Get_Option_Value("-only_bdy_outside"));
            levelset_maker.Keep_Only_Largest_Inside_Region(parse_args.Get_Option_Value("-keep_largest"));
            levelset_maker.Flip_Sign_If_Corners_Are_Inside(parse_args.Get_Option_Value("-flip"));
            levelset_maker.Set_Positive_Boundary_Band(parse_args.Get_Integer_Value("-positive_boundary"));
            levelset_maker.Set_Phi_Offset((T)parse_args.Get_Double_Value("-offset"));
            ARRAY<T,VECTOR<int,3> > phi(grid.Domain_Indices());
            levelset_maker.Compute_Level_Set(*triangulated_surface,grid,phi);
            LEVELSET_IMPLICIT_OBJECT<TV> levelset_implicit_surface(grid,phi);
            //phi+=(T)1*grid.Maximum_Edge_Length();
        FILE_UTILITIES::Write_To_File<RW>(output_filename,levelset_implicit_surface);}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    else{
        LEVELSET_MAKER_DYADIC<T> levelset_maker;
        levelset_maker.Verbose_Mode(true);
        levelset_maker.Write_Debug_Data(parse_args.Get_Option_Value("-debug"));
        if(parse_args.Is_Value_Set("-path_start")&&parse_args.Is_Value_Set("-path_end"))
            levelset_maker.Write_Debug_Path(true,VECTOR<int,3>(parse_args.Get_Vector_3D_Value("-path_start")),VECTOR<int,3>(parse_args.Get_Vector_3D_Value("-path_end")));
        if(parse_args.Get_Option_Value("-orthogonal_vote")) levelset_maker.Use_Orthogonal_Vote(true);
        levelset_maker.Set_Surface_Padding_For_Flood_Fill((T)parse_args.Get_Double_Value("-padding"));
        levelset_maker.Set_Surface_Thickness((T)parse_args.Get_Double_Value("-thickness"));
        if(parse_args.Get_Option_Value("-heaviside"))
            levelset_maker.Compute_Heaviside_Function();
        else{
            if(parse_args.Get_Option_Value("-unsigned"))
                levelset_maker.Compute_Unsigned_Distance_Function();
            else
                levelset_maker.Compute_Signed_Distance_Function();
            levelset_maker.Use_Fast_Marching_Method(true,(T)parse_args.Get_Double_Value("-band"));}
            levelset_maker.Only_Boundary_Region_Is_Outside(parse_args.Get_Option_Value("-only_bdy_outside"));
            levelset_maker.Keep_Only_Largest_Inside_Region(parse_args.Get_Option_Value("-keep_largest"));
            levelset_maker.Flip_Sign_If_Corners_Are_Inside(parse_args.Get_Option_Value("-flip"));
            levelset_maker.Set_Positive_Boundary_Band(parse_args.Get_Integer_Value("-positive_boundary"));
            levelset_maker.Set_Phi_Offset((T)parse_args.Get_Double_Value("-offset"));
            ARRAY<T> phi,phi_nodes;
            levelset_maker.Compute_Level_Set(*triangulated_surface,octree_grid,phi,&phi_nodes,0,coarsen_bandwidth,Depth,&helper);
            DYADIC_IMPLICIT_OBJECT<TV> implicit_surface(octree_grid,phi,&phi_nodes);
            FILE_UTILITIES::Write_To_File<T>("octree_grid.0",octree_grid);
        FILE_UTILITIES::Write_To_File<T>("octree_levelset.0",phi);
        FILE_UTILITIES::Write_To_File<T>("octree_levelset_nodes.0",phi_nodes);
        FILE_UTILITIES::Write_To_File<T>(output_filename,octree_grid,phi);}
#endif
}

int main(int argc,char *argv[])
{
    Initialize_Geometry_Particle();
    Initialize_Read_Write_General_Structures();
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false;
    VECTOR<double,3> grid_size(50,50,50);
    int boundary_cells=3;
    int positive_boundary_band=0;

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float","data is in float format");
    parse_args.Add_Option_Argument("-double","data is in double format");
    parse_args.Add_Option_Argument("-heaviside","compute heaviside function");
    parse_args.Add_Option_Argument("-unsigned","compute unsigned distance");
    parse_args.Add_Option_Argument("-mac","use MAC grid");
    parse_args.Add_Option_Argument("-only_bdy_outside","only boundary region is outside");
    parse_args.Add_Option_Argument("-keep_largest","keep only largest inside region");
    parse_args.Add_Option_Argument("-flip","flip sign if corners are inside");
    parse_args.Add_Option_Argument("-debug","write debug data");
    parse_args.Add_Vector_3D_Argument("-g",grid_size,"grid size","suggested grid size - actual size may be larger");
    parse_args.Add_Double_Argument("-dx",0,"grid dx","specify grid cell size");
    parse_args.Add_Vector_3D_Argument("-domain_min",VECTOR<double,3>());
    parse_args.Add_Vector_3D_Argument("-domain_max",VECTOR<double,3>());
    parse_args.Add_Double_Argument("-coarsen",0,"coarsen_band_width","coarsen band width for fast marching method (in grid cells)");
    parse_args.Add_Double_Argument("-band",0,"half_band_width","half band width for fast marching method (in grid cells)");
    parse_args.Add_Double_Argument("-padding",0,"padding for flood fill");
    parse_args.Add_Double_Argument("-thickness",1e-5,"surface thickness");
    parse_args.Add_Double_Argument("-offset",0,"subtract offset from phi");
    parse_args.Add_Integer_Argument("-b",boundary_cells,"boundary cells","number of cells outside bounding box");
    parse_args.Add_Integer_Argument("-positive_boundary",positive_boundary_band,"positive boundary band","ensure this band around boundary is positive");
    parse_args.Add_Option_Argument("-octree","generate an octree levelset");
    parse_args.Add_Integer_Argument("-depth",1,"maximum_depth","maximum depth if using an octree");
    parse_args.Add_String_Argument("-o","","output filename","output_filename");
    parse_args.Add_String_Argument("-exact_grid","");
    parse_args.Add_Vector_3D_Argument("-path_start",VECTOR<double,3>());
    parse_args.Add_Vector_3D_Argument("-path_end",VECTOR<double,3>());
    parse_args.Add_Option_Argument("-orthogonal_vote");
    parse_args.Add_Option_Argument("-compute_using_doubles");
    parse_args.Set_Extra_Arguments(1,"<tri file>","<tri file> tri file to convert");

    parse_args.Parse(argc, argv);

    if(parse_args.Is_Value_Set("-b")) boundary_cells=parse_args.Get_Integer_Value("-b");
    std::string input_filename=parse_args.Extra_Arg(1);

    if(parse_args.Get_Option_Value("-float")) type_double=false;
    if(parse_args.Get_Option_Value("-double")) type_double=true;

    if(!FILE_UTILITIES::Is_Tri_File(input_filename)){
        std::cerr<<"Not a tri file: "<<input_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(parse_args.Get_Option_Value("-compute_using_doubles")){
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
