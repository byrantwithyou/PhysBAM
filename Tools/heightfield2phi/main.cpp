//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Math_Tools/minmag.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_MAKER.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Heat_Flows/HEAT_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_RLE_Computations/DUALCONTOUR_RLE_3D.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_SURFACE_MAKER.h>

using namespace PhysBAM;
using namespace FILE_UTILITIES;

//#####################################################################
// Function Triangulate
//#####################################################################
template<class T> void Triangulate(const GRID<VECTOR<T,2> >& image_grid,const ARRAYS<VECTOR<T,2> >& heights,TRIANGLE_MESH& triangle_mesh,DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<T,3> >& particles,const T extension_distance)
{
    bool extend=extension_distance!=0;
    triangle_mesh.Initialize_Square_Mesh(image_grid.counts.x+2*extend,image_grid.counts.y+2*extend);
    particles.Preallocate(triangle_mesh.number_nodes);
    for(int j=1-extend;j<=image_grid.counts.y+extend;j++)for(int i=1-extend;i<=image_grid.counts.x+extend;i++)
        particles.X(particles.Add_Particle())=VECTOR<T,3>(image_grid.Axis_X(i,1),heights(i,j),image_grid.Axis_X(j,2));
    if(extend){
        T thickness_over_2=image_grid.dX.x/2;
        T multiplier=extension_distance/image_grid.dX.x;
        BOX<VECTOR<T,2> > domain=image_grid.Domain();
        for(int p=0;p<particles.number;p++){
            VECTOR<T,2> horizontal_X=particles.X(p).Horizontal_Vector();
            if(domain.Outside(horizontal_X,thickness_over_2)){
                VECTOR<T,2> dX=multiplier*(horizontal_X-domain.Surface(horizontal_X));
                particles.X(p).x+=dX.x;particles.X(p).z+=dX.y;}}}
}
//#####################################################################
// Function Rasterize_Uniform
//#####################################################################
template<class T> void Rasterize_Uniform(const GRID<VECTOR<T,2> >& image_grid,const ARRAYS<VECTOR<T,2> >& heights_ghost,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset)
{
    LOG::Time("rasterizing");
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,T> interpolation;
    levelset.phi.Resize(levelset.grid.Domain_Indices(),false,false);
    for(int i=0;i<levelset.phi.m;i++)for(int ij=0;ij<levelset.phi.mn;ij++){
        VECTOR<T,3> X=levelset.grid.Center(i,1,ij);
        T height=interpolation.Clamped_To_Array(image_grid,heights_ghost,X.Horizontal_Vector());
        for(int j=0;j<levelset.phi.n;j++) levelset.phi(i,j,ij)=levelset.grid.Axis_X(j,2)-height;}
    LOG::Time("reinitializing");
    levelset.Fast_Marching_Method();
}
//#####################################################################
// Function Rasterize_RLE
//#####################################################################
template<class T> struct HELPER_IMPLICIT_SURFACE:public NONCOPYABLE
{
    const GRID<VECTOR<T,2> >& grid;
    ARRAYS<VECTOR<BOX<VECTOR<T,1> >,2> > height_boxes;
    int bandwidth;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,T> interpolation;
    
    HELPER_IMPLICIT_SURFACE(const GRID<VECTOR<T,2> >& grid_input,const ARRAYS<VECTOR<T,2> >& heights,const int bandwidth_input)
        :grid(grid_input),height_boxes(grid.Domain_Indices(3)),bandwidth(bandwidth_input)
    {
        for(int i=heights.domain.min_corner.x;i<=heights.domain.max_corner.x;i++)for(int j=heights.domain.min_corner.y;j<=heights.domain.max_corner.y;j++){
            height_boxes(i,j)=BOX<VECTOR<T,1> >(VECTOR<T,1>(heights(i,j)));
            for(int ii=-bandwidth;ii<=bandwidth;ii++)for(int jj=-bandwidth;jj<=bandwidth;jj++)if(heights.Valid_Index(i+ii,j+jj))
                height_boxes(i,j).Enlarge_To_Include_Point(VECTOR<T,1>(heights(i+ii,j+jj)));}
    }
    
    T operator()(const VECTOR<T,3>& X) const
    {VECTOR<int,2> index=INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,BOX<VECTOR<T,1> > >::Clamped_Index(grid,height_boxes,X.Horizontal_Vector());
    BOX<VECTOR<T,1> > box=VECTOR<T,1>(X.y)-height_boxes(index);
    return box.Lazy_Inside(VECTOR<T,1>((T)0))?0:minmag(box.min_corner.x,box.max_corner.x);}
};
template<class T> void Rasterize_RLE(const GRID<VECTOR<T,2> >& image_grid,const ARRAYS<VECTOR<T,2> >& heights_ghost,RLE_GRID_3D<T>& grid,LEVELSET_RLE<RLE_GRID_3D<T> >& levelset)
{
    typedef RLE_GRID_3D<T> T_GRID;
    LOG::Time("building grid");
    int helper_bandwidth=(int)ceil(max(grid.negative_bandwidth,grid.positive_bandwidth)/grid.Minimum_Edge_Length());
    LOG::cout<<"helper bandwidth = "<<helper_bandwidth<<std::endl;
    HELPER_IMPLICIT_SURFACE<T> helper(image_grid,heights_ghost,helper_bandwidth);
    grid.Initialize(helper,0);
    LOG::cout<<"initial grid: cell = "<<grid.number_of_cells<<std::endl;
    ARRAY<T>& phi=levelset.phi;
    phi.Resize(grid.number_of_cells);

    LOG::Time("rasterizing");
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,T> interpolation;
    for(typename T_GRID::CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){int c=cell.Cell();
        VECTOR<T,3> X=cell.Center();
        phi(c)=X.y-interpolation.Clamped_To_Array(image_grid,heights_ghost,X.Horizontal_Vector());
        if(cell.Long()) phi(c+1)=phi(c);}

    LOG::Time("reinitializing");
    levelset.Fast_Marching_Method(0);

    LOG::Time("rebuilding grid");
    ARRAY<bool> cell_should_be_long(grid.number_of_cells);
    for(typename T_GRID::CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++)if(cell.Short()){int c=cell.Cell();
        cell_should_be_long(c)=phi(c)>=grid.positive_bandwidth || phi(c)<=-grid.negative_bandwidth;}
    T_GRID new_grid;
    T_GRID::Rebuild(grid,new_grid,cell_should_be_long,1,0);
    LOG::cout<<"new grid: cells = "<<new_grid.number_of_cells<<std::endl;
    levelset.Transfer_Phi(new_grid);
    T_GRID::Transfer(new_grid,grid);
    grid.Adjust_Vertical_Space();
}
//#####################################################################
// Function Convert
//#####################################################################
template<class T,class RW> void Convert(PARSE_ARGS& parse_args,const std::string& input_file,const std::string& output_base)
{
    T ymax=(T)parse_args.Get_Double_Value("-ymax");
    T slope=(T)parse_args.Get_Double_Value("-slope");
    bool flip=parse_args.Get_Option_Value("-flip");
    T extension_distance=parse_args.Get_Double_Value("-extend");
    VECTOR<T,2> extend_upwards(parse_args.Get_Vector_2D_Value("-extend_upwards"));
    VECTOR<T,3> color(parse_args.Get_Vector_3D_Value("-color"));
    bool already_heightfield=parse_args.Get_Option_Value("-already_heightfield");

    GRID<VECTOR<T,2> > image_grid;
    ARRAYS<VECTOR<T,2> > heights;

    if(!already_heightfield){
        // read image
        ARRAYS<VECTOR<VECTOR<T,3> ,2> > image;IMAGE<T>::Read(input_file,image);
        // figure out grids
        VECTOR<T,2> size(image.m,image.n);size/=size.Min();
        image_grid=GRID<VECTOR<T,2> >(image.m,image.n,BOX<VECTOR<T,2> >(VECTOR<T,2>(),size));
        LOG::cout<<"image grid = "<<image_grid<<std::endl;

        LOG::Time("making heightfield");
        heights.Resize(image_grid.Domain_Indices(),false,false);
        for(typename GRID<VECTOR<T,2> >::NODE_ITERATOR iterator(image_grid);iterator.Valid();iterator.Next()){VECTOR<int,2> node=iterator.Node_Index();
            T value=VECTOR<T,3>::Dot_Product(image(node),color);
            if(flip) value=1-value;
            if(extend_upwards.x<value){
                T t=(value-extend_upwards.x)/(1-extend_upwards.x);
                value=(1-t)*value+t*((1-t)*extend_upwards.x+t*extend_upwards.y);}
            heights(node)=value*ymax+slope*iterator.Location().x;}

        LOG::Time("smoothing");
        int smoothing_steps=parse_args.Get_Integer_Value("-smooth");
        LOG::cout<<"steps = "<<smoothing_steps<<std::endl;
        HEAT_UNIFORM<GRID<VECTOR<T,2> > >::Smooth(heights,smoothing_steps);

        LOG::Time("writing heightfield");
        FILE_UTILITIES::Write_To_File<RW>(output_base+"_grid",image_grid);
        FILE_UTILITIES::Write_To_File<RW>(output_base+"_height",heights);}
    else{
        FILE_UTILITIES::Read_From_File<RW>(input_file+"_grid",image_grid);
        FILE_UTILITIES::Read_From_File<RW>(input_file+"_height",heights);}
        
    ARRAYS<VECTOR<T,2> > heights_ghost(image_grid.Domain_Indices(3),false);
    BOUNDARY_UNIFORM<GRID<VECTOR<T,2> >,T> boundary;boundary.Fill_Ghost_Cells(image_grid,heights,heights_ghost,0,0);

    LOG::Time("triangulating");
    TRIANGLE_MESH triangle_mesh;DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<T,3> > particles;
    TRIANGULATED_SURFACE<T> triangulated_surface(triangle_mesh,particles);
    Triangulate(image_grid,heights_ghost,triangle_mesh,particles,extension_distance);
    LOG::cout<<"triangles = "<<triangle_mesh.elements.m<<", particles = "<<particles.number<<std::endl;
    FILE_UTILITIES::Write_To_File<RW>(output_base+".tri",triangulated_surface);

    LOG::SCOPE scope("RASTERIZE","rasterizing");
    T min_height=heights.Min(),max_height=heights.Max();
    BOX<VECTOR<T,3> > box(0,image_grid.domain.max_corner.x,min_height,max_height,0,image_grid.domain.max_corner.y);
    GRID<VECTOR<T,3> > grid=GRID<VECTOR<T,3> >::Create_Grid_Given_Cell_Size(box,image_grid.dX.x,false,5);
    if(!parse_args.Get_Option_Value("-rle")){
        LOG::cout<<"rasterizing onto uniform grid "<<grid<<std::endl;
        ARRAYS<VECTOR<T,3> > phi;LEVELSET_3D<GRID<VECTOR<T,3> > > levelset(grid,phi);
        Rasterize_Uniform(image_grid,heights_ghost,levelset);
        FILE_UTILITIES::Write_To_File<RW>(output_base+".phi",levelset);}
    else{
        LOG::cout<<"rasterizing onto rle grid "<<grid<<std::endl;
        int negative_bandwidth=parse_args.Get_Integer_Value("-negative_bandwidth");
        int positive_bandwidth=parse_args.Get_Integer_Value("-positive_bandwidth");
        RLE_GRID_3D<T> rle_grid;
        rle_grid.Set_Uniform_Grid(grid);
        rle_grid.Set_Positive_Bandwidth_In_Cells(positive_bandwidth);
        rle_grid.Set_Negative_Bandwidth_In_Cells(negative_bandwidth);
        LOG::cout<<"rle positive bandwidth = "<<rle_grid.positive_bandwidth<<", negative_bandwidth = "<<rle_grid.negative_bandwidth<<std::endl;
        rle_grid.Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile();
        ARRAY<T> phi;LEVELSET_RLE<RLE_GRID_3D<T> > levelset(rle_grid,phi);
        Rasterize_RLE(image_grid,heights_ghost,rle_grid,levelset);
        FILE_UTILITIES::Write_To_File<RW>(output_base+"_rle_grid",rle_grid);
        FILE_UTILITIES::Write_To_File<RW>(output_base+"_rle_levelset",phi);
        FILE_UTILITIES::Write_To_File<RW>(output_base+"_rle_grid_and_levelset",levelset);
        LOG::Time("contouring");
        TRIANGULATED_SURFACE<T>* contour_surface=DUALCONTOUR_RLE_3D<T>::Create_Triangulated_Surface_From_Levelset(levelset);
        FILE_UTILITIES::Write_To_File<RW>(output_base+"_rle_levelset.tri",*contour_surface);
        delete contour_surface;}
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char *argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-o","","output base",std::string(""));
    parse_args.Add_Integer_Argument("-smooth",0,"smoothing steps");
    parse_args.Add_Double_Argument("-ymax",1,"maximum height");
    parse_args.Add_Double_Argument("-slope",0,"slope");
    parse_args.Add_Option_Argument("-flip");
    parse_args.Add_Option_Argument("-rle");
    parse_args.Add_Integer_Argument("-negative_bandwidth",3,"rle levelset negative_bandwidth");
    parse_args.Add_Integer_Argument("-positive_bandwidth",3,"rle levelset positive_bandwidth");
    parse_args.Add_Double_Argument("-extend",0,"triangulated surface extension distance");
    parse_args.Add_Vector_2D_Argument("-extend_upwards",VECTOR<double,2>(2,1),"apply nonlinear function to grow top of heightfield");
    parse_args.Add_Vector_3D_Argument("-color",VECTOR<double,3>(1,0,0),"color to map to height");
    parse_args.Add_Option_Argument("-already_heightfield");
    parse_args.Set_Extra_Arguments(1,"<heighfield>","heighfield to convert");
    parse_args.Parse();

    std::string input_file,output_base;
    input_file=parse_args.Extra_Arg(0);
    if(parse_args.Is_Value_Set("-o")) output_base=parse_args.Get_String_Value("-o");
    else output_base=Get_Basename(input_file);

    Convert<float,float>(parse_args,input_file,output_base);
    LOG::Finish_Logging();
    return 0;
}
//#####################################################################
