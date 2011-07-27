//#####################################################################
// Copyright 2005-2006, William Fong, Michael Lentine, Andrew Selle, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include "VH_LEVELSET_BUILDER.h"
#include <Fl/Fl_Hold_Browser.h>
using namespace PhysBAM;
//#####################################################################
// Function VH_LEVELSET_BUILDER
//#####################################################################
template<class T> VH_LEVELSET_BUILDER<T>::
VH_LEVELSET_BUILDER(std::string data_directory_in)
    :data_directory(data_directory_in),box_padding_voxels(0),box_padding((T)0),include_all_tissues(true),levelset(grid,phi)
{
    FILE_UTILITIES::Read_From_File<T>(data_directory+"/dimensions",image_resolution,image_size);
    std::cout<<"Got data set dimensions resolution=("<<image_resolution<<") cell size=("<<image_size<<")"<<std::endl;
    FILE_UTILITIES::Read_From_File<T>(data_directory+"/labels",labels);
    Compute_Bounding_Boxes();
    //for(int i=0;i<0xffff;i++) if(labels(i)!="") std::cout<<"i="<<i<<" "<<labels(i)<<std::endl;
}
//#####################################################################
// Function Compute_Bounding_Boxes
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Compute_Bounding_Boxes()
{
    std::cout<<"Trying to read tissue bounding box"<<std::endl;
    try{FILE_UTILITIES::Read_From_File<T>(data_directory+"/tissue_bounding_boxes",bounding_boxes);}
    catch(FILESYSTEM_ERROR&){
        std::cout<<"Compute Bounding Boxes..."<<std::endl;
        ARRAYS<VECTOR<bool,1> > initialized(0,0xffff,false);initialized.Fill(false);
        bounding_boxes.Resize(0,0xffff,false,false);bounding_boxes.Fill((RANGE<VECTOR<int,3> >(-1,-1,-1,-1,-1,-1)));
        for(int j=1;j<=image_resolution.y;j++){ 
            std::cout<<"Slice "<<j<<" of "<<image_resolution.y<<std::endl;
            Read_Slice(j);
            for(int i=1;i<=image_resolution.x;i++) for(int ij=1;ij<=image_resolution.z;ij++){
                if(initialized(image(i,ij))) bounding_boxes(image(i,ij)).Enlarge_To_Include_Point(VECTOR<int,3>(i,j,ij));
                else{bounding_boxes(image(i,ij)).Reset_Bounds(VECTOR<int,3>(i,j,ij));initialized(image(i,ij))=true;}}}
        FILE_UTILITIES::Write_To_File<T>(data_directory+"tissue_bounding_boxes",bounding_boxes);}
}
//#####################################################################
// Function Create_Levelset_Interp
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Create_Levelset_Interp(const ARRAY<int>& tissues,RANGE<VECTOR<int,3> >& box,VECTOR<int,3>& granularity,const std::string filename)
{
    assert(tissues.m>0);
    for(int t=1;t<=tissues.m;t++){
        if(bounding_boxes(tissues(t)).min_corner.x==-1){std::cout<<"Fatal tissue "<<tissues(t)<<" has degenerate bounding box"<<std::endl;exit(1);}
        else std::cout<<"Tissue "<<t<<" is "<<labels(tissues(t))<<std::endl;}
    
    //BOX_3D<int> box=bounding_boxes(tissues(1));
    //for(int t=2;t<=tissues.m;t++) box.Enlarge_To_Include_Box(bounding_boxes(tissues(t)));
    std::cout<<"Building level set from image bounded by "<<box<<std::endl;
    Initialize_Grid(box,granularity);
    std::cout<<"Granularity: "<<grid.m<<" "<<grid.n<<" "<<grid.mn<<std::endl;
    phi.Resize(grid.Domain_Indices(),false,false);
    phi.Fill(grid.min_dx_dy_dz);
    T image_start,image_x_next,image_y_next,image_z_next,image_xy_next,image_yz_next,image_xz_next,image_xyz_next;
    VECTOR<T,3> alpha,interp,interp_next;
    VECTOR<int,3> image_index,image_index_next,image_index_offset(box.min_corner),phi_point,levelset_bounds(box.max_corner-box.min_corner);
    VECTOR<double,3> count_index(1,1,1);
    VECTOR<double,3> scale((double)levelset_bounds.x/grid.m,(double)levelset_bounds.y/grid.n,(double)levelset_bounds.z/grid.mn);
    std::cout<<grid.m<<" "<<grid.n<<" "<<grid.mn<<" "<<scale<<std::endl;
    for(count_index.y=box.min_corner.y;count_index.y<=box.max_corner.y;count_index.y+=scale.y){
        image_index.y=(int)floor(count_index.y);
        image_index_next.y=image_index.y+1;
        alpha.y=1-(count_index.y-image_index.y);
        std::cout<<"Slice "<<count_index.y-box.min_corner.y+1<<" of "<<box.max_corner.y-box.min_corner.y+1<<" "<<std::endl;
        Read_Slice(image_index.y,image);
        Read_Slice(image_index_next.y,image_next);
        for(count_index.x=box.min_corner.x;count_index.x<=box.max_corner.x;count_index.x+=scale.x) for(count_index.z=box.min_corner.z;count_index.z<=box.max_corner.z;count_index.z+=scale.z){
            image_index.z=(int)floor(count_index.z);image_index.x=(int)floor(count_index.x);
            image_index_next.z=image_index.z+1;image_index_next.x=image_index.x+1;
            alpha.z=1-(count_index.z-image_index.z);alpha.x=1-(count_index.x-image_index.x);
            phi_point.x=max((int)((count_index.x-image_index_offset.x)/scale.x),1);phi_point.y=max((int)((count_index.y-image_index_offset.y)/scale.y),1);phi_point.z=max((int)((count_index.z-image_index_offset.z)/scale.z),1);
            //Find Interp Value (Trilinear interpolation)
            image_start=image_x_next=image_y_next=image_z_next=image_xy_next=image_yz_next=image_xz_next=image_xyz_next=grid.min_dx_dy_dz;
            for(int t=1;t<=tissues.m;t++){
                if(image(image_index.x,image_index.z)==tissues(t)) image_start=-grid.min_dx_dy_dz;
                if(image(image_index_next.x,image_index.z)==tissues(t)) image_x_next=-grid.min_dx_dy_dz;
                if(image(image_index.x,image_index_next.z)==tissues(t)) image_z_next=-grid.min_dx_dy_dz;
                if(image(image_index_next.x,image_index_next.z)==tissues(t)) image_xz_next=-grid.min_dx_dy_dz;
                if(image_next(image_index.x,image_index.z)==tissues(t)) image_y_next=-grid.min_dx_dy_dz;
                if(image_next(image_index_next.x,image_index.z)==tissues(t)) image_xy_next=-grid.min_dx_dy_dz;
                if(image_next(image_index.x,image_index_next.z)==tissues(t)) image_yz_next=-grid.min_dx_dy_dz;
                if(image_next(image_index_next.x,image_index_next.z)==tissues(t)) image_xyz_next=-grid.min_dx_dy_dz;}
            interp.x=alpha.x*image_start+(1-alpha.x)*image_x_next;
            interp_next.x=alpha.x*image_z_next+(1-alpha.x)*image_xz_next;
            interp.z=alpha.z*interp.x+(1-alpha.z)*interp_next.x;
            interp.x=alpha.x*image_y_next+(1-alpha.x)*image_xy_next;
            interp_next.x=alpha.x*image_yz_next+(1-alpha.x)*image_xyz_next;
            interp_next.z=alpha.z*interp.x+(1-alpha.z)*interp_next.x;
            interp.y=alpha.y*interp.z+(1-alpha.y)*interp_next.z;
            //std::cout<<"Phi at "<<phi_point<<" is "<<interp.y<<std::endl;
            if(phi_point.z<=phi.mn_end&&phi_point.y<=phi.n_end&&phi_point.x<=phi.m_end)phi(phi_point)=interp.y;}}
    Close_Boundaries();
    std::cout<<"Making signed distance.."<<std::endl;
    levelset.Fast_Marching_Method(0,grid.dx*8);
    std::cout<<"Writing to file..."<<std::endl;
    FILE_UTILITIES::Write_To_File<T>(filename,levelset);
}
//#####################################################################
// Function Create_Levelset
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Create_Levelset(const ARRAY<int>& tissues,RANGE<VECTOR<int,3> >& box,VECTOR<int,3>& granularity,const std::string filename)
{
    assert(tissues.m>0);
    for(int t=1;t<=tissues.m;t++){
        if(bounding_boxes(tissues(t)).min_corner.x==-1){std::cout<<"Fatal tissue "<<tissues(t)<<" has degenerate bounding box"<<std::endl;exit(1);}
        else std::cout<<"Tissue "<<t<<" is "<<labels(tissues(t))<<std::endl;}
    
    //BOX_3D<int> box=bounding_boxes(tissues(1));
    //for(int t=2;t<=tissues.m;t++) box.Enlarge_To_Include_Box(bounding_boxes(tissues(t)));
    std::cout<<"Building level set from image bounded by "<<box<<std::endl;
    Initialize_Grid(box,granularity);
    std::cout<<"Granularity: "<<grid.m<<" "<<grid.n<<" "<<grid.mn<<std::endl;
    phi.Resize(grid.Domain_Indices(),false,false);
    phi.Fill(grid.min_dx_dy_dz);
    VECTOR<int,3> image_index,image_index_offset(box.min_corner),phi_point,levelset_bounds(box.max_corner-box.min_corner);
    VECTOR<double,3> count_index(1,1,1);
    VECTOR<double,3> scale((double)levelset_bounds.x/grid.m,(double)levelset_bounds.y/grid.n,(double)levelset_bounds.z/grid.mn);
    std::cout<<grid.m<<" "<<grid.n<<" "<<grid.mn<<" "<<scale<<std::endl;
    for(count_index.y=box.min_corner.y;count_index.y<=box.max_corner.y;count_index.y+=scale.y){
        image_index.y=(int)floor(count_index.y);
        std::cout<<"Slice "<<image_index.y-box.min_corner.y+1<<" of "<<box.max_corner.y-box.min_corner.y+1<<" "<<std::endl;
        Read_Slice(image_index.y);
        for(count_index.x=box.min_corner.x;count_index.x<=box.max_corner.x;count_index.x+=scale.x) for(count_index.z=box.min_corner.z;count_index.z<=box.max_corner.z;count_index.z+=scale.z){
            image_index.z=(int)floor(count_index.z);image_index.x=(int)floor(count_index.x);
            phi_point.x=max((int)((image_index.x-image_index_offset.x)/scale.x),1);phi_point.y=max((int)((image_index.y-image_index_offset.y)/scale.y),1);phi_point.z=max((int)((image_index.z-image_index_offset.z)/scale.z),1);
            if(phi_point.z<=phi.mn_end&&phi_point.y<=phi.n_end&&phi_point.x<=phi.m_end)phi(phi_point)=grid.min_dx_dy_dz;
            if(image(image_index.x,image_index.z)==tissues(9179)&&image(image_index.x,image_index.z)==tissues(680)) continue;
            for(int t=1;t<=tissues.m;t++) if(image(image_index.x,image_index.z)==tissues(t)){
                if(phi_point.z<=phi.mn_end&&phi_point.y<=phi.n_end&&phi_point.x<=phi.m_end) phi(phi_point)=-grid.min_dx_dy_dz;break;}}}
    Close_Boundaries();
    std::cout<<"Making signed distance.."<<std::endl;
    levelset.Fast_Marching_Method(0,grid.dx*8);
    std::cout<<"Writing to file..."<<std::endl;
    FILE_UTILITIES::Write_To_File<T>(filename,levelset);
}
//#####################################################################
// Function Close_Boundaries
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Close_Boundaries()
{
    VECTOR<int,3> phi_point;
    for(int i=phi.m_start;i<phi.m_end;i++) for(int j=phi.mn_start;j<phi.mn_end;j++)
    {phi_point.x=i;phi_point.y=phi.n_start;phi_point.z=j;phi(phi_point)=grid.min_dx_dy_dz;
    phi_point.y=phi.n_end;phi(phi_point)=grid.min_dx_dy_dz;}
    for(int i=phi.n_start;i<phi.n_end;i++) for(int j=phi.mn_start;j<phi.mn_end;j++)
    {phi_point.x=phi.m_start;phi_point.y=i;phi_point.z=j;phi(phi_point)=grid.min_dx_dy_dz;
    phi_point.x=phi.m_end;phi(phi_point)=grid.min_dx_dy_dz;}
    for(int i=phi.n_start;i<phi.n_end;i++) for(int j=phi.m_start;j<phi.m_end;j++)
    {phi_point.x=j;phi_point.y=i;phi_point.z=phi.mn_start;phi(phi_point)=grid.min_dx_dy_dz;
    phi_point.z=phi.mn_end;phi(phi_point)=grid.min_dx_dy_dz;}
}
//#####################################################################
// Function All_Tissues
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
All_Tissues(Fl_Hold_Browser* tissues, HASHTABLE<int,int>& tissues_hash)
{
#ifdef USE_FLTK
   for (int i = 1; i < bounding_boxes.m; i++){
        if (bounding_boxes(i).min_corner.x != -1){
            int dummy;
            if(!tissues_hash.Get(i,dummy)){
                std::string s=STRING_UTILITIES::string_sprintf("(%d) %s",i,labels(i).c_str());
                tissues->add(s.c_str(),(void*)i);
                tissues_hash.Insert(i,0);}}}
#endif
}
//#####################################################################
// Function All_Tissues_Array
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
All_Tissues_Array(ARRAY<int>& tissues)
{
    for (int i = 1; i < bounding_boxes.m; i++)
        if (bounding_boxes(i).min_corner.x != -1) tissues.Append(i);
}
//#####################################################################
// Function Read_Slice
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Read_Slice(int slice_j)
{
    assert(slice_j >= 1 && slice_j <= image_resolution.y);
    FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/slice%05d.img",data_directory.c_str(),slice_j),image);
}
//#####################################################################
// Function Read_Slice
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Read_Slice(int slice_j,ARRAYS<VECTOR<unsigned short,2> >& image)
{
    assert(slice_j >= 1 && slice_j <= image_resolution.y);
    FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/slice%05d.img",data_directory.c_str(),slice_j),image);
}
//#####################################################################
// Function Dump_Images
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Dump_Images()
{
    ARRAYS<VECTOR<VECTOR<T,3> ,2> > image;
    VECTOR<int,3> image_index;
    for(image_index.y=1;image_index.y<=image_resolution.y;image_index.y++){
        std::cout<<"Slice "<<image_index.y<<" of "<<image_resolution.y<<std::endl;
        Read_Slice(image_index.y);
        Convert_To_Image(image);
        IMAGE<T>::Write(STRING_UTILITIES::string_sprintf("slice%05d.jpg",image_index.y),image);}
}
//#####################################################################
// Function Convert_To_Image
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Convert_To_Image(ARRAYS<VECTOR<VECTOR<T,3> ,2> >& float_image)
{
    float_image.Resize(1,image_resolution.x,1,image_resolution.z);

    static VECTOR<unsigned short,3> masks(0x001f,0x07e0,0xf800);
    static VECTOR<T,3> maxes=VECTOR<T,3>(1<<5,1<<6,1<<5)-VECTOR<T,3>(1,1,1);
    static VECTOR<T,3> color_scale=VECTOR<T,3>(255,255,255)/maxes;
    static T color_map_factor=(T)1/(T)256;

    for(int i=1;i<=float_image.m;i++) for(int j=1;j<=float_image.n;j++){
        float_image(i,j)=color_map_factor*(color_scale*VECTOR<T,3>((image(i,j)&masks.x),(image(i,j)&masks.y)>>5,(image(i,j)&masks.z)>>11));
        float_image(i,j).x-=floor(float_image(i,j).x);float_image(i,j).y-=floor(float_image(i,j).y);float_image(i,j).z-=floor(float_image(i,j).z);}
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Initialize_Grid(RANGE<VECTOR<int,3> >& image_box,VECTOR<int,3>& grid_size)
{
    grid.Initialize(grid_size.x,grid_size.y,grid_size.z,RANGE<VECTOR<T,3> >(image_size*VECTOR<T,3>(image_box.Minimum_Corner()-box_padding_voxels*VECTOR<int,3>(1,1,1)),
                                                                  image_size*VECTOR<T,3>(image_box.Maximum_Corner()+box_padding_voxels*VECTOR<int,3>(1,1,1))));
}
//#####################################################################
// Function Get_Grid_Size
//#####################################################################
template<class T> VECTOR<int,3> VH_LEVELSET_BUILDER<T>::
Get_Grid_Size(RANGE<VECTOR<int,3> >& image_box)
{
    VECTOR<int,3> grid_size=image_box.Edge_Lengths()+VECTOR<int,3>(1,1,1)+2*box_padding_voxels*VECTOR<int,3>(1,1,1);
    return grid_size;
}
//#####################################################################
// Function Get_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<int,3> > VH_LEVELSET_BUILDER<T>::
Get_Bounding_Box(const ARRAY<int>& tissues)
{
    if (tissues.m>0){
        for(int t=1;t<=tissues.m;t++){
            if(bounding_boxes(tissues(t)).min_corner.x==-1){std::cout<<"Fatal tissue "<<tissues(t)<<" has degenerate bounding box"<<std::endl;exit(1);}}
        RANGE<VECTOR<int,3> > box=bounding_boxes(tissues(1));
        for(int t=2;t<=tissues.m;t++) box.Enlarge_To_Include_Box(bounding_boxes(tissues(t)));return box;}
    else{std::cout<<"No slices selected"<<std::endl;return RANGE<VECTOR<int,3> >(0,0,0,0,0,0);}
}
//#####################################################################
// Function Tissue_Index
//#####################################################################
template<class T> int VH_LEVELSET_BUILDER<T>::
Tissue_Index(std::string tissue_label)
{
    for(int i=1;i<=labels.m;i++) if(labels(i)==tissue_label) return i;
    return 0;
}
//#####################################################################
// Function Make_Sub_Data_Set
//#####################################################################
template<class T> void VH_LEVELSET_BUILDER<T>::
Make_Sub_Data_Set(std::string output_directory,RANGE<VECTOR<int,3> > box)
{
    int sub_slice=1;
    VECTOR<int,3> size=box.Maximum_Corner()-box.Minimum_Corner()+VECTOR<int,3>(1,1,1);
    std::cout<<"Bounding box "<<box<<std::endl;
    std::cout<<"Dimension "<<size<<std::endl;

    ARRAYS<VECTOR<unsigned short,2> > sub_image(1,size.x,1,size.z);
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Write_To_File<T>(output_directory+"/dimensions",size,image_size);
    for(int slice=box.min_corner.y;slice<=box.max_corner.y;slice++){
        std::cout<<"Doing slice "<<slice<<" which will be new slice "<<sub_slice<<std::endl;
        Read_Slice(slice);
        for(int i=1;i<=size.x;i++) for(int ij=1;ij<=size.z;ij++) sub_image(i,ij)=image(i+box.min_corner.x-1,ij+box.min_corner.z-1);
        FILE_UTILITIES::Write_To_File<T>(STRING_UTILITIES::string_sprintf("%s/slice%05d.img",output_directory.c_str(),sub_slice),sub_image);
        sub_slice++;}
}
//#####################################################################
template class VH_LEVELSET_BUILDER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VH_LEVELSET_BUILDER<double>;
#endif
//#####################################################################
