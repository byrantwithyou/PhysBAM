//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Grid_Tools/Arrays/FLOOD_FILL.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Geometry/Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_MAKER.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Process_Segment
//#####################################################################
static void Process_Segment(const int m,const ARRAY_VIEW<bool,VECTOR<int,3> >& edge_is_blocked,const ARRAY<bool,VECTOR<int,3> >& is_inside,const VECTOR<int,3>& start_index,const int axis,ARRAY<char,VECTOR<int,3> >& vote)
{
    typedef VECTOR<int,3> TV_INT;
    TV_INT index=start_index,increment;
    increment[axis]=1;
    int segment_start=1;
    bool segment_starts_inside=false;
    for(index[axis]=0;index[axis]<m;index[axis]++){
        if(index[axis]>0 && edge_is_blocked(index)){segment_start=index[axis];segment_starts_inside=is_inside(index);}
        if(index[axis]==m-1 || (index[axis]<m-1 && edge_is_blocked(index+increment))){
            int segment_end=index[axis];
            bool segment_ends_inside=index[axis]<m?is_inside(index):false;
            int vote_increment=(int)segment_starts_inside+(int)segment_ends_inside;
            TV_INT t=start_index;
            for(t[axis]=segment_start;t[axis]<=segment_end;t[axis]++) vote(t)+=vote_increment;}}
}
//#####################################################################
// Function Compute_Level_Set
//#####################################################################
// Assumes triangulated surface is closed
template<class T> bool LEVELSET_MAKER<T>::
Compute_Level_Set(TRIANGULATED_SURFACE<T>& triangulated_surface,GRID<TV>& grid,ARRAY<T,TV_INT>& phi,ARRAY<TV,TV_INT>* velocity)
{
    PHYSBAM_ASSERT(grid.counts.x && grid.counts.y && grid.counts.z,"Attempting to rasterize surface onto grid with no indices");
    if(remove_degenerate_triangles_area_threshold){
        if(verbose) LOG::Time("Remove Degenerate Triangles");
        triangulated_surface.Remove_Degenerate_Triangles(remove_degenerate_triangles_area_threshold);}

    // initialize acceleration structures
    if(verbose) LOG::Time("Initialize Acceleration Structures");
    bool triangle_list_defined=triangulated_surface.triangle_list!=0;if(!triangulated_surface.triangle_list) triangulated_surface.Update_Triangle_List();
    bool hierarchy_defined=triangulated_surface.hierarchy!=0;if(!triangulated_surface.hierarchy) triangulated_surface.Initialize_Hierarchy();
    bool bounding_box_defined=triangulated_surface.bounding_box!=0;if(!triangulated_surface.bounding_box) triangulated_surface.Update_Bounding_Box();
    bool adjacent_triangles_defined=triangulated_surface.mesh.adjacent_elements!=0;if(!adjacent_triangles_defined) triangulated_surface.mesh.Initialize_Adjacent_Elements();
    bool incident_triangles_defined=triangulated_surface.mesh.incident_elements!=0;if(!incident_triangles_defined) triangulated_surface.mesh.Initialize_Incident_Elements();
 
    bool compute_velocity=velocity && triangulated_surface.particles.store_velocity;
    if(compute_velocity) velocity->Resize(grid.Domain_Indices());

    phi.Resize(grid.Domain_Indices(),false);phi.Fill(FLT_MAX);
   
    if(use_fmm && compute_velocity && extrapolate_velocity && fmm_one_sided_band_width && fmm_one_sided_band_width<velocity_extrapolation_one_sided_band_width+1){
        LOG::cerr<<"Extending FMM band width to one more than velocity extrapolation band width"<<std::endl;
        fmm_one_sided_band_width=velocity_extrapolation_one_sided_band_width+1;}
    T fmm_stopping_distance=fmm_one_sided_band_width*grid.dX.Max();

    bool need_flood_fill=compute_signed_distance_function || compute_heaviside_function;
    ARRAY<bool,FACE_INDEX<TV::m> > edge_is_blocked;
    if(need_flood_fill) edge_is_blocked.Resize(grid);

    bool store_closest_triangle_index=need_flood_fill && !only_boundary_region_is_outside;
    ARRAY<int,TV_INT> closest_triangle_index;
    if(store_closest_triangle_index) closest_triangle_index.Resize(grid.Domain_Indices());

    bool store_initialized_indices=use_fmm;
    if(store_initialized_indices){initialized_indices.Exact_Resize(0);initialized_indices.Preallocate(20);}

    const T surface_thickness_over_two=Surface_Thickness_Over_Two(grid),surface_padding_for_flood_fill=Surface_Padding_For_Flood_Fill(grid);

    const RANGE<TV>& grid_domain=grid.domain;
    if(verbose) LOG::Time("Rasterizing Triangles");
    for(int t=0;t<triangulated_surface.mesh.elements.m;t++){
        const TRIANGLE_3D<T>& triangle=(*triangulated_surface.triangle_list)(t);
        TRIANGLE_3D<T> enlarged_triangle=triangle;if(surface_padding_for_flood_fill) enlarged_triangle.Change_Size(surface_padding_for_flood_fill);
        RANGE<TV> triangle_bounding_box=enlarged_triangle.Bounding_Box();
        triangle_bounding_box.Change_Size(surface_thickness_over_two);
        if(!grid_domain.Lazy_Intersection(triangle_bounding_box)) continue;
        RANGE<TV_INT> range(grid.Clamped_Index(triangle_bounding_box.Minimum_Corner()),grid.Clamped_Index_End_Minus_One(triangle_bounding_box.Maximum_Corner())+2);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            TV grid_position=grid.X(it.index),weights,closest_point=triangle.Closest_Point(grid_position,weights);
            T distance_squared=(grid_position-closest_point).Magnitude_Squared();
            if(phi(it.index)==FLT_MAX || distance_squared<sqr(phi(it.index))){
                if(store_initialized_indices && phi(it.index)==FLT_MAX) initialized_indices.Append(it.index);
                phi(it.index)=sqrt(distance_squared);
                if(store_closest_triangle_index) closest_triangle_index(it.index)=t;
                if(compute_velocity){
                    int node1,node2,node3;triangulated_surface.mesh.elements(t).Get(node1,node2,node3);
                    (*velocity)(it.index)=weights.x*triangulated_surface.particles.V(node1)+weights.y*triangulated_surface.particles.V(node2)+weights.z*triangulated_surface.particles.V(node3);}}}
        if(need_flood_fill){
            for(int d=0;d<TV::m;d++){
                RANGE<TV_INT> range_d(range);
                range_d.min_corner(d)++;
                for(RANGE_ITERATOR<TV::m> it(range_d);it.Valid();it.Next())
                    if(!edge_is_blocked.Component(d)(it.index)){
                        TV_INT ind(it.index);
                        ind(d)--;
                        edge_is_blocked.Component(d)(it.index)=INTERSECTION::Intersects(SEGMENT_3D<T>(grid.X(it.index),grid.X(ind)),enlarged_triangle,surface_thickness_over_two);}}}}

    if((compute_signed_distance_function || compute_unsigned_distance_function) && use_fmm && fmm_stopping_distance)
        for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) phi(it.index)=min(phi(it.index),fmm_stopping_distance);
    else if(compute_heaviside_function) phi.Fill(grid.dX.Max());

    if(write_debug_data){
        GRID<TV> output_grid=grid;
        if(!grid.Is_MAC_Grid())
            output_grid=GRID<TV>(grid.counts,RANGE<TV>(grid.domain.min_corner-grid.dX/2,grid.domain.max_corner+grid.dX/2),true);
        // TODO: put this back if you need it
        /*FILE_UTILITIES::Write_To_File<T>("grid.debug",output_grid);
          FILE_UTILITIES::Write_To_File<T>("triangulated_surface.debug",triangulated_surface);*/}

    if(use_orthogonal_vote){
        // TODO: put this back if you need it
        /*if(write_debug_data){FILE_UTILITIES::Write_To_File<T>("edge_is_blocked.debug",edge_is_blocked.Component(0),edge_is_blocked.Component(1),edge_is_blocked.Component(2));}*/
        if(verbose) LOG::Time("Computing sign using orthogonal vote");
        ARRAY<bool,TV_INT> is_inside(grid.Domain_Indices());
        for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
            if(closest_triangle_index(it.index)) is_inside(it.index)=triangulated_surface.Inside_Relative_To_Triangle(grid.X(it.index),closest_triangle_index(it.index),surface_thickness_over_two);}
        ARRAY<char,TV_INT> vote(grid.Domain_Indices());
        for(int j=0;j<grid.counts.y;j++) for(int k=0;k<grid.counts.z;k++) Process_Segment(grid.counts.x,edge_is_blocked.Component(0),is_inside,TV_INT(0,j,k),0,vote);
        for(int i=0;i<grid.counts.x;i++) for(int k=0;k<grid.counts.z;k++) Process_Segment(grid.counts.y,edge_is_blocked.Component(1),is_inside,TV_INT(i,0,k),1,vote);
        for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++) Process_Segment(grid.counts.z,edge_is_blocked.Component(2),is_inside,TV_INT(i,j,0),2,vote);
        is_inside.Clean_Memory();
        for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) if(vote(it.index)>=3) phi(it.index)*=-1;

        if(keep_only_largest_inside_region){
            ARRAY<int,TV_INT> colors(grid.Domain_Indices());colors.Fill(-1);ARRAY<bool,FACE_INDEX<TV::m> > null_edge_is_blocked(grid,1);
            for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) if(phi(it.index)>0) colors(it.index)=-2; // make outside regions uncolorable
            FLOOD_FILL<TV::m> flood_fill;flood_fill.Optimize_Fill_For_Single_Cell_Regions(true);
            int number_of_colors=flood_fill.Flood_Fill(colors,null_edge_is_blocked);
            // TODO: put this back if you need it
            //if(write_debug_data){FILE_UTILITIES::Write_To_File<T>("colors.debug",colors);}
            ARRAY<int> region_size(number_of_colors);
            for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) if(colors(it.index)>0) region_size(colors(it.index))++;
            int max_region_size=region_size.Max();
            if(verbose) LOG::cout<<"Keeping only largest inside region (max region size = "<<max_region_size<<")... "<<std::flush;
            // flip smaller regions back to positive sign
            for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) if(colors(it.index)>0 && region_size(colors(it.index))<max_region_size) phi(it.index)*=-1;}}
    else if(need_flood_fill){ // Need flood fill to determine sign (inside/outside)
        if(verbose) LOG::Time("Flood Fill");
        ARRAY<int,TV_INT> colors(grid.Domain_Indices());colors.Fill(-1);
        FLOOD_FILL<TV::m> flood_fill;flood_fill.Optimize_Fill_For_Single_Cell_Regions(true);
        int number_of_colors=flood_fill.Flood_Fill(colors,edge_is_blocked);
        if(verbose) LOG::cout<<"(got "<<number_of_colors<<" colors)... "<<std::endl;
        if(number_of_colors==1 && !phi_offset){ // there is only one color. check if the whole domain is inside or outside then return
            if(triangulated_surface.Inside(grid.X(TV_INT(1,1,1)))) phi.Fill(-FLT_MAX);
            else phi.Fill(FLT_MAX);
            return false;}
        // TODO: put this back if you need it
        /*if(write_debug_data){
            FILE_UTILITIES::Write_To_File<T>("colors.debug",colors);
            FILE_UTILITIES::Write_To_File<T>("edge_is_blocked.debug",edge_is_blocked.Component(0),edge_is_blocked.Component(1),edge_is_blocked.Component(2));}
        if(write_debug_path){
            ARRAY<TV_INT> path_nodes;
            bool path_exists=flood_fill.Path_Between_Nodes(RANGE<TV_INT>(1,grid.counts.x,1,grid.counts.y,1,grid.counts.z),path_start_node,path_end_node,
                                                           edge_is_blocked.Component(0),edge_is_blocked.Component(1),edge_is_blocked.Component(2),&path_nodes);
            if(verbose){LOG::cout<<"Path between "<<path_start_node<<" and "<<path_end_node<<" "<<(path_exists?"exists":"doesn't exist")<<std::endl;}
            ARRAY<bool,TV_INT> path(grid.Domain_Indices());
            for(int i=0;i<path_nodes.m;i++){path(path_nodes(i))=true;}
            FILE_UTILITIES::Write_To_File<T>("path.debug",path);}*/
        ARRAY<bool,VECTOR<int,1> > color_is_inside(0,number_of_colors);
        if(only_boundary_region_is_outside){
            ARRAY<bool> color_touches_boundary(number_of_colors);
            if(verbose) LOG::Time("Marking boundary region as outside");
            flood_fill.Identify_Colors_Touching_Boundary(number_of_colors,colors,edge_is_blocked.data,color_touches_boundary);
            if(verbose && color_touches_boundary.Number_True()>1) LOG::cerr<<"Warning: Got "<<color_touches_boundary.Number_True()<<" colors touching boundary"<<std::endl;
            for(int i=0;i<number_of_colors;i++) color_is_inside(i)=!color_touches_boundary(i);}
        else{
            ARRAY<T> color_maximum_distance(number_of_colors,false);color_maximum_distance.Fill(-1);
            ARRAY<TV_INT> color_representatives(number_of_colors);
            for(int i=0;i<grid.counts.x;i++)
                for(int j=0;j<grid.counts.y;j++)
                    for(int k=0;k<grid.counts.z;k++)
                        if(closest_triangle_index(i,j,k)>=0 && color_maximum_distance(colors(i,j,k))<phi(i,j,k)){
                            color_maximum_distance(colors(i,j,k))=phi(i,j,k);
                            color_representatives(colors(i,j,k))=TV_INT(i,j,k);}
            for(int color=0;color<number_of_colors;color++){
                if(color_maximum_distance(color)<0){LOG::cerr<<"Error: could not determine inside/outside for color "<<color<<std::endl;PHYSBAM_FATAL_ERROR();}
                else color_is_inside(color)=triangulated_surface.Inside_Relative_To_Triangle(grid.X(color_representatives(color)),
                    closest_triangle_index(color_representatives(color)),surface_thickness_over_two);}}
        if(keep_only_largest_inside_region){
            ARRAY<int> region_size(number_of_colors);
            for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) if(colors(it.index)>0 && color_is_inside(colors(it.index))) region_size(colors(it.index))++;
            int max_region_size=region_size.Max();
            if(verbose) LOG::cout<<"Keeping only largest inside region (max region size = "<<max_region_size<<")... "<<std::flush;
            for(int i=0;i<number_of_colors;i++) if(color_is_inside(i) && region_size(i)<max_region_size) color_is_inside(i)=false;}
        if(flip_sign_if_corners_are_inside){ // If the majority of corners are labelled as inside then we flip signs
            int num_corners_inside=(int)color_is_inside(colors(1,1,1))+(int)color_is_inside(colors(1,1,grid.counts.z))+
                                   (int)color_is_inside(colors(1,grid.counts.y,1))+(int)color_is_inside(colors(1,grid.counts.y,grid.counts.z))+
                                   (int)color_is_inside(colors(grid.counts.x,1,1))+(int)color_is_inside(colors(grid.counts.x,1,grid.counts.z))+
                                   (int)color_is_inside(colors(grid.counts.x,grid.counts.y,1))+(int)color_is_inside(colors(grid.counts.x,grid.counts.y,grid.counts.z));
            if(num_corners_inside>4){
                if(verbose) LOG::cout<<"Majority of corners are inside -- flipping sign!"<<std::endl;
                for(int i=0;i<number_of_colors;i++) color_is_inside(i)=!color_is_inside(i);}}
        for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) if(color_is_inside(colors(it.index))) phi(it.index)*=-1;}

    if(positive_boundary_band){
        RANGE<TV> clip(grid.X(TV_INT()+positive_boundary_band)+(T).5*grid.dX,grid.X(grid.counts-positive_boundary_band)+(T).5*grid.dX);
        for(int j=0;j<grid.counts.y;j++)
            for(int k=0;k<grid.counts.z;k++){
                for(int i=0;i<positive_boundary_band+1;i++) phi(i,j,k)=max(phi(i,j,k),clip.min_corner.x-grid.X(TV_INT(i,j,k)).x);
                for(int i=grid.counts.x-positive_boundary_band;i<grid.counts.x;i++) phi(i,j,k)=max(phi(i,j,k),grid.X(TV_INT(i,j,k)).x-clip.max_corner.x);}

        for(int i=0;i<grid.counts.x;i++)
            for(int k=0;k<grid.counts.z;k++){
                for(int j=0;j<positive_boundary_band+1;j++) phi(i,j,k)=max(phi(i,j,k),clip.min_corner.y-grid.X(TV_INT(i,j,k)).y);
                for(int j=grid.counts.y-positive_boundary_band;j<grid.counts.y;j++) phi(i,j,k)=max(phi(i,j,k),grid.X(TV_INT(i,j,k)).y-clip.max_corner.y);}

        for(int i=0;i<grid.counts.x;i++)
            for(int j=0;j<grid.counts.y;j++){
                for(int k=0;k<positive_boundary_band+1;k++) phi(i,j,k)=max(phi(i,j,k),clip.min_corner.z-grid.X(TV_INT(i,j,k)).z);
                for(int k=grid.counts.z-positive_boundary_band;k<grid.counts.z;k++) phi(i,j,k)=max(phi(i,j,k),grid.X(TV_INT(i,j,k)).z-clip.max_corner.z);}}

    if(use_fmm && (compute_unsigned_distance_function || compute_signed_distance_function)){
        for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) 
            phi(it.index)=clamp(phi(it.index),-10*grid.dX.Min(),10*grid.dX.Min()); // clamp away from FLT_MAX to avoid floating point exceptions
        if(verbose) LOG::Time(LOG::sprintf("Fast Marching (one sided band width=%f)",fmm_one_sided_band_width));
        GRID<TV> grid_copy=grid;LEVELSET<TV> levelset(grid_copy,phi);
        if(compute_unsigned_distance_function) levelset.Fast_Marching_Method(0,fmm_stopping_distance,&initialized_indices);
        else if(compute_signed_distance_function) levelset.Fast_Marching_Method(0,fmm_stopping_distance,phi_offset?&initialized_indices:0);}

    if(compute_velocity && extrapolate_velocity){
        if(!compute_signed_distance_function || !use_fmm){LOG::cerr<<"Can only extrapolate velocity if computing signed distance function"<<std::endl;}
        else{
            if(verbose) LOG::Time(LOG::sprintf("Extrapolating velocity (one sided band width=%f)",velocity_extrapolation_one_sided_band_width));
            GRID<TV> grid_copy=grid;
            EXTRAPOLATION_UNIFORM<TV,TV> extrapolation(grid_copy,phi,*velocity,3);
            extrapolation.Set_Custom_Seed_Indices(&initialized_indices);
            extrapolation.Set_Band_Width(velocity_extrapolation_one_sided_band_width);
            extrapolation.Extrapolate(0);
            phi*=-1;
            extrapolation.Extrapolate(0);
            phi*=-1;}}

    if(phi_offset){
        phi-=phi_offset;
        if(use_fmm && compute_signed_distance_function)
            LEVELSET<TV>(grid,phi).Fast_Marching_Method(0,fmm_stopping_distance);}

    // TODO: put this back if you need it
    /*if(write_debug_data){
        GRID<TV> grid_copy=grid;LEVELSET<TV> levelset(grid_copy,phi);
        FILE_UTILITIES::Write_To_File<T>("levelset.debug",levelset);
        if(compute_velocity)FILE_UTILITIES::Write_To_File<T>("velocity.debug",velocity);}*/

    if(verbose) LOG::cout<<"Done"<<std::endl;

    // delete acceleration structures if defined in this function
    if(!incident_triangles_defined){delete triangulated_surface.mesh.incident_elements;triangulated_surface.mesh.incident_elements=0;}
    if(!adjacent_triangles_defined){delete triangulated_surface.mesh.adjacent_elements;triangulated_surface.mesh.adjacent_elements=0;}
    if(!bounding_box_defined){delete triangulated_surface.bounding_box;triangulated_surface.bounding_box=0;}
    if(!hierarchy_defined){delete triangulated_surface.hierarchy;triangulated_surface.hierarchy=0;}
    if(!triangle_list_defined){delete triangulated_surface.triangle_list;triangulated_surface.triangle_list=0;}
    return true;
}
//#####################################################################
namespace PhysBAM{
template class LEVELSET_MAKER<float>;
template class LEVELSET_MAKER<double>;
}
