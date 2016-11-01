//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <OpenGL/OpenGL/OPENGL_TEXTURED_RECT.h>
#include <OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> OPENGL_SCALAR_FIELD_3D<T,T2>::
OPENGL_SCALAR_FIELD_3D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,ARRAY<T2,VECTOR<int,3> > &values_input,OPENGL_COLOR_MAP<T2> *color_map_input,DRAW_MODE draw_mode_input)
    :OPENGL_OBJECT<T>(stream_type),grid(grid_input),values(values_input),current_color_map(0),opengl_textured_rect(0),opengl_points(0),smooth_slice_texture(false),scale_range(false),selected_cell(-1,-1,-1),selected_node(-1,-1,-1),selected_point(-1)
{
    PHYSBAM_ASSERT(color_map_input);
    Initialize_Color_Maps(color_map_input);
    Set_Draw_Mode(draw_mode_input);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2> OPENGL_SCALAR_FIELD_3D<T,T2>::
~OPENGL_SCALAR_FIELD_3D()
{
    Delete_Textured_Rect();
    Delete_Points();
    color_maps.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<float,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<double,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<float,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<double,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input)
{
    color_maps.Append(color_map_input);
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Jet(8e4,1e6));
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Hot(8e4,1e6));
}
//#####################################################################
// Function Delete_Textured_Rect
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Delete_Textured_Rect()
{
    if(opengl_textured_rect){delete opengl_textured_rect->texture;}
    delete opengl_textured_rect;opengl_textured_rect=0;
}
//#####################################################################
// Function Delete_Points
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Delete_Points()
{
    if(opengl_points) delete &opengl_points->points;
    delete opengl_points;opengl_points=0;
}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<float,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<double,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<float,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<double,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Set_Scale_Range(const T2 range_min,const T2 range_max)
{
    scale_range=true;
    scale_range_min=range_min;
    T2 range_length=(range_max-range_min);
    scale_range_dx=range_length>0?(T2)1/range_length:(T2)0;
}
//#####################################################################
// Function Reset_Scale_Range
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Reset_Scale_Range()
{
    scale_range=false;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<> bool OPENGL_SCALAR_FIELD_3D<float,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<> bool OPENGL_SCALAR_FIELD_3D<double,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}

//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<class T,class T2> T2 OPENGL_SCALAR_FIELD_3D<T,T2>::
Pre_Map_Value(const T2 value) const
{
    if(!scale_range) return value;
    else return (value-scale_range_min)*scale_range_dx; 
}
//#####################################################################
// Function Do_Color
//#####################################################################
template<class T,class T2> OPENGL_COLOR OPENGL_SCALAR_FIELD_3D<T,T2>::
Do_Color(const int i, const int j, const int k) const
{
    return color_maps(current_color_map)->Lookup(Pre_Map_Value(values(i,j,k)));
}
//#####################################################################
// Function Do_Color
//#####################################################################
template<class T,class T2> OPENGL_COLOR OPENGL_SCALAR_FIELD_3D<T,T2>::
Do_Color(const TV_INT& index) const
{
    return color_maps(current_color_map)->Lookup(Pre_Map_Value(values(index)));
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Display() const
{
    if(values.domain.Empty()) return;
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)this->slice;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    if(draw_mode==DRAW_TEXTURE){
        if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE) Display_3D();
        else if(opengl_textured_rect)
        opengl_textured_rect->Display();
    }
    else{
        PHYSBAM_ASSERT(opengl_points);
        if(slice && slice->Is_Slice_Mode()){
            glPushAttrib(GL_ENABLE_BIT);
            slice->Enable_Clip_Planes();}
        opengl_points->Display();
        if(slice && slice->Is_Slice_Mode()){
            glPopAttrib();}}

    glPopMatrix();
}
//#####################################################################
// Function Display_3D
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Display_3D() const
{
    glPushAttrib(GL_ENABLE_BIT|GL_DEPTH_BUFFER_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);

    TV view_forward,view_up,view_right;
    OPENGL_WORLD<T>::Singleton()->Get_View_Frame(view_forward,view_up,view_right);
    int dominant_axis=1;

    if(dominant_axis==0){
        if(view_forward[0]>0){
            for(int i=grid.counts.x-1;i>=0;i--) for(int j=0;j<grid.counts.y;j++) for(int k=0;k<grid.counts.z;k++){
                OpenGL_Begin(GL_TRIANGLE_STRIP);
                OpenGL_Color(Do_Color(i,j,k).rgba);
                TV pos=grid.X(TV_INT(i,j,k));
                OpenGL_Vertex(TV(pos.x,pos.y-0.5*grid.dX.y,pos.z-0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x,pos.y-0.5*grid.dX.y,pos.z+0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x,pos.y+0.5*grid.dX.y,pos.z-0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x,pos.y+0.5*grid.dX.y,pos.z+0.5*grid.dX.z));
                OpenGL_End();}}
        else{
            for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
                OpenGL_Begin(GL_TRIANGLE_STRIP);
                OpenGL_Color(Do_Color(it.index).rgba);
                TV pos=grid.X(it.index);
                OpenGL_Vertex(TV(pos.x,pos.y-0.5*grid.dX.y,pos.z+0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x,pos.y-0.5*grid.dX.y,pos.z-0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x,pos.y+0.5*grid.dX.y,pos.z+0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x,pos.y+0.5*grid.dX.y,pos.z-0.5*grid.dX.z));
                OpenGL_End();}}}
    else if(dominant_axis==1){
        if(view_forward[1]>0){
            for(int j=grid.counts.y-1;j>=0;j--) for(int i=0;i<grid.counts.x;i++) for(int k=0;k<grid.counts.z;k++){
                OpenGL_Begin(GL_TRIANGLE_STRIP);
                OpenGL_Color(Do_Color(i,j,k).rgba);
                TV pos=grid.X(TV_INT(i,j,k));
                OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z));
                OpenGL_End();}}
        else{
            for(int j=0;j<grid.counts.y;j++) for(int i=0;i<grid.counts.x;i++) for(int k=0;k<grid.counts.z;k++){
                OpenGL_Begin(GL_TRIANGLE_STRIP);
                OpenGL_Color(Do_Color(i,j,k).rgba);
                TV pos=grid.X(TV_INT(i,j,k));
                OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z));
                OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z));
                OpenGL_End();}}}
    else if(dominant_axis==2){
        if(view_forward[2]>0){
            for(int k=grid.counts.z-1;k>=0;k--) for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
                OpenGL_Begin(GL_TRIANGLE_STRIP);
                OpenGL_Color(Do_Color(i,j,k).rgba);
                TV pos=grid.X(TV_INT(i,j,k));
                OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z));
                OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z));
                OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z));
                OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z));
                OpenGL_End();}}
        else{
            for(int k=0;k<grid.counts.z;k++) for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
                OpenGL_Begin(GL_TRIANGLE_STRIP);
                OpenGL_Color(Do_Color(i,j,k).rgba);
                TV pos=grid.X(TV_INT(i,j,k));
                OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z));
                OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z));
                OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z));
                OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z));
                OpenGL_End();}}}

    glPopAttrib();
}
//#####################################################################
// Function Display_3D_Slice
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Display_3D_Slice() const
{
    glPushAttrib(GL_ENABLE_BIT|GL_DEPTH_BUFFER_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);

    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)this->slice;

    if(slice->axis==0){
        int i=slice->index;
        for(int j=0;j<grid.counts.y;j++) for(int k=0;k<grid.counts.z;k++){
            OpenGL_Begin(GL_TRIANGLE_STRIP);
            OpenGL_Color(Do_Color(i,j,k).rgba);
            TV pos=grid.X(TV_INT(i,j,k));
            OpenGL_Vertex(TV(pos.x,pos.y-0.5*grid.dX.y,pos.z-0.5*grid.dX.z));
            OpenGL_Vertex(TV(pos.x,pos.y-0.5*grid.dX.y,pos.z+0.5*grid.dX.z));
            OpenGL_Vertex(TV(pos.x,pos.y+0.5*grid.dX.y,pos.z-0.5*grid.dX.z));
            OpenGL_Vertex(TV(pos.x,pos.y+0.5*grid.dX.y,pos.z+0.5*grid.dX.z));
            OpenGL_End();}}
    else if(slice->axis==1){
        int j=slice->index;
        for(int i=0;i<grid.counts.x;i++) for(int k=0;k<grid.counts.z;k++){
            OpenGL_Begin(GL_TRIANGLE_STRIP);
            OpenGL_Color(Do_Color(i,j,k).rgba);
            TV pos=grid.X(TV_INT(i,j,k));
            OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z));
            OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z));
            OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z));
            OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z));
            OpenGL_End();}}
    else if(slice->axis==2){
        int k=slice->index;
        for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
            OpenGL_Begin(GL_TRIANGLE_STRIP);
            OpenGL_Color(Do_Color(i,j,k).rgba);
            TV pos=grid.X(TV_INT(i,j,k));
            OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z));
            OpenGL_Vertex(TV(pos.x-0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z));
            OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z));
            OpenGL_Vertex(TV(pos.x+0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z));
            OpenGL_End();}}

    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    if(slice && slice->Is_Slice_Mode() && opengl_textured_rect) return opengl_textured_rect->Bounding_Box();
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Function Set_Draw_Mode
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Set_Draw_Mode(DRAW_MODE draw_mode_input)
{
    draw_mode=draw_mode_input;

    if(draw_mode==DRAW_TEXTURE){
        Delete_Points();
        if(!opengl_textured_rect) opengl_textured_rect=new OPENGL_TEXTURED_RECT<T>(stream_type);}
    else{
        Delete_Textured_Rect();
        if(!opengl_points) opengl_points=new OPENGL_POINTS_3D<T>(stream_type,*new ARRAY<TV>);}

    Update();
}
//#####################################################################
// Function Update
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Update()
{
    if(draw_mode==DRAW_TEXTURE){if(slice && slice->Is_Slice_Mode()) Update_Slice();}
    else Update_Points();
}
//#####################################################################
// Function Print_Selection_Info_Helper
//#####################################################################
template<class T2,class TV> static void
Print_Selection_Info_Helper(std::ostream& output_stream,const TV& location,const GRID<TV>& grid,ARRAY<T2,VECTOR<int,3> >& values)
{
    output_stream<<" @ particle = "<<LINEAR_INTERPOLATION_UNIFORM<TV,T2>().Clamped_To_Array(grid,values,location);
}
//#####################################################################
// Function Print_Selection_Info_Helper
//#####################################################################
// no interpolation for bool's and int's
template<class TV> static void
Print_Selection_Info_Helper(std::ostream& output_stream,const TV& location,const GRID<TV>&,ARRAY<bool,VECTOR<int,3> >& values){}
//#####################################################################
// Function Print_Selection_Info_Helper
//#####################################################################
template<class TV> static void
Print_Selection_Info_Helper(std::ostream& output_stream,const TV& location,const GRID<TV>&,ARRAY<int,VECTOR<int,3> >& values){}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Print_Selection_Info(std::ostream& output_stream) const
{
    if(selected_cell.x>=0 && grid.Is_MAC_Grid() && values.Valid_Index(selected_cell)) output_stream<<values(selected_cell);
    if(selected_node.x>=0 && !grid.Is_MAC_Grid() && values.Valid_Index(selected_node)) output_stream<<values(selected_node);
    // if(current_selection && current_selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_3D){
    //     OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T> *selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T>*)current_selection;
    //     Print_Selection_Info_Helper(output_stream,selection,grid,values);}
    output_stream<<std::endl;
}
namespace{
//#####################################################################
// Function Update_Slice_Helper
//#####################################################################
template<class T> void
Update_Slice_Helper(OPENGL_SCALAR_FIELD_3D<T,bool>* self,int tex_width,int tex_height){
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)self->slice;
    VECTOR<int,3> domain_start(self->values.domain.min_corner.x,self->values.domain.min_corner.y,self->values.domain.min_corner.z),domain_end(self->values.domain.max_corner.x,self->values.domain.max_corner.y,self->values.domain.max_corner.z);
    if(!self->opengl_textured_rect->texture || self->opengl_textured_rect->texture->width!=tex_width || self->opengl_textured_rect->texture->height!=tex_height){
        delete self->opengl_textured_rect->texture;
        self->opengl_textured_rect->texture=new OPENGL_TEXTURE();
        self->opengl_textured_rect->texture->Initialize(tex_width,tex_height);
        self->opengl_textured_rect->texture->Set_Smooth_Shading(self->smooth_slice_texture);}

    OPENGL_COLOR* bitmap=new OPENGL_COLOR[tex_width*tex_height];
    OPENGL_COLOR_MAP<bool>* color_map=self->color_maps(self->current_color_map);
    for(int i=0;i<tex_width;i++)
        for(int j=0;j<tex_height;j++){
            bool value=bool();
            switch (slice->axis){
                case 0: value=self->values(slice->index,domain_start.y+j,domain_start.z+tex_width-i);break;
                case 1: value=self->values(domain_start.x+i,slice->index,domain_start.z+tex_height-j);break;
                case 2: value=self->values(domain_start.x+i,domain_start.y+j,slice->index);break;}
            int idx=j*tex_width+i;
            bitmap[idx]=color_map->Lookup(self->Pre_Map_Value(value));}

    self->opengl_textured_rect->texture->Update_Texture(bitmap);

    delete[]bitmap;
}
//#####################################################################
// Function Update_Slice_Helper
//#####################################################################
template<class T> void
Update_Slice_Helper(OPENGL_SCALAR_FIELD_3D<T,int>* self,int tex_width,int tex_height){
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)self->slice;
    VECTOR<int,3> domain_start(self->values.domain.min_corner.x,self->values.domain.min_corner.y,self->values.domain.min_corner.z),domain_end(self->values.domain.max_corner.x,self->values.domain.max_corner.y,self->values.domain.max_corner.z);
    if(!self->opengl_textured_rect->texture || self->opengl_textured_rect->texture->width!=tex_width || self->opengl_textured_rect->texture->height!=tex_height){
        delete self->opengl_textured_rect->texture;
        self->opengl_textured_rect->texture=new OPENGL_TEXTURE();
        self->opengl_textured_rect->texture->Initialize(tex_width,tex_height);
        self->opengl_textured_rect->texture->Set_Smooth_Shading(self->smooth_slice_texture);}

    OPENGL_COLOR* bitmap=new OPENGL_COLOR[tex_width*tex_height];
    OPENGL_COLOR_MAP<int>* color_map=self->color_maps(self->current_color_map);
    for(int i=0;i<tex_width;i++)
        for(int j=0;j<tex_height;j++){
            int value=int();
            switch (slice->axis){
                case 0: value=self->values(slice->index,domain_start.y+j,domain_start.z+tex_width-i);break;
                case 1: value=self->values(domain_start.x+i,slice->index,domain_start.z+tex_height-j);break;
                case 2: value=self->values(domain_start.x+i,domain_start.y+j,slice->index);break;}
            int idx=j*tex_width+i;
            bitmap[idx]=color_map->Lookup(self->Pre_Map_Value(value));}

    self->opengl_textured_rect->texture->Update_Texture(bitmap);

    delete[]bitmap;
}
//#####################################################################
// Function Update_Slice_Helper
//#####################################################################
template<class T,class T2> void
Update_Slice_Helper(OPENGL_SCALAR_FIELD_3D<T,T2>* self,int tex_width,int tex_height){
    typedef VECTOR<T,3> TV;
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)self->slice;
    int k=1;
    tex_width=k*tex_width;
    tex_height=k*tex_height;

    if(!self->opengl_textured_rect->texture || self->opengl_textured_rect->texture->width!=tex_width || self->opengl_textured_rect->texture->height!=tex_height){
        delete self->opengl_textured_rect->texture;
        self->opengl_textured_rect->texture=new OPENGL_TEXTURE();
        self->opengl_textured_rect->texture->Initialize(tex_width,tex_height);
        self->opengl_textured_rect->texture->Set_Smooth_Shading(self->smooth_slice_texture);}

    OPENGL_COLOR* bitmap=new OPENGL_COLOR[tex_width*tex_height];
    OPENGL_COLOR_MAP<T2>* color_map=self->color_maps(self->current_color_map);

    LINEAR_INTERPOLATION_UNIFORM<TV,T2> interpolation;
    for(int i=0;i<tex_width;i++)
        for(int j=0;j<tex_height;j++){
            T2 value=T2();
            TV location;
            switch (slice->axis){
                case 0: location=TV(self->grid.X(VECTOR<int,3>(slice->index,0,0)).x,
                    self->grid.domain.min_corner.y+j*(self->grid.domain.max_corner.y-self->grid.domain.min_corner.y)/tex_height,
                    self->grid.domain.max_corner.z-i*(self->grid.domain.max_corner.z-self->grid.domain.min_corner.z)/tex_width);
                    break;
                case 1: location=TV(self->grid.domain.min_corner.x+i*(self->grid.domain.max_corner.x-self->grid.domain.min_corner.x)/tex_width,
                    self->grid.X(VECTOR<int,3>(0,slice->index,0)).y,
                    self->grid.domain.max_corner.z-j*(self->grid.domain.max_corner.z-self->grid.domain.min_corner.z)/tex_height);
                    break;
                case 2: location=TV(self->grid.domain.min_corner.x+i*(self->grid.domain.max_corner.x-self->grid.domain.min_corner.x)/tex_width,
                    self->grid.domain.min_corner.y+j*(self->grid.domain.max_corner.y-self->grid.domain.min_corner.y)/tex_height,
                    self->grid.X(VECTOR<int,3>(0,0,slice->index)).z);
                    break;}
            value=interpolation.Clamped_To_Array(self->grid,self->values,location);
            int idx=j*tex_width+i;
            bitmap[idx]=color_map->Lookup(self->Pre_Map_Value(value));}

    self->opengl_textured_rect->texture->Update_Texture(bitmap);
    delete []bitmap;
}
}
//#####################################################################
// Function Update_Slice
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Update_Slice()
{
    PHYSBAM_ASSERT(this->slice);
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)this->slice;
    VECTOR<int,3> domain_start(values.domain.min_corner.x,values.domain.min_corner.y,values.domain.min_corner.z),domain_end(values.domain.max_corner.x,values.domain.max_corner.y,values.domain.max_corner.z);
    if((slice->mode==OPENGL_SLICE::CELL_SLICE && (!grid.Is_MAC_Grid() || slice->index<domain_start[slice->axis] || slice->index>=domain_end[slice->axis])) ||
        (slice->mode==OPENGL_SLICE::NODE_SLICE && (grid.Is_MAC_Grid() || slice->index<domain_start[slice->axis] || slice->index>=domain_end[slice->axis]))){
        // Currently we don't draw anything ifthe slice doesn't match where the scalar field lives
        Delete_Textured_Rect();
        return;}

    if(!opengl_textured_rect) opengl_textured_rect=new OPENGL_TEXTURED_RECT<T>(stream_type);

    // Handle values arrays which are not (1,m)(1,n)
    TV half_dX=(T)0.5*grid.dX;
    RANGE<VECTOR<int,3> > domain_indices(values.Domain_Indices());
    RANGE<TV> domain(grid.X(domain_indices.min_corner)-half_dX,grid.X(domain_indices.max_corner)+half_dX);
    opengl_textured_rect->frame->t=domain.Center();

    // rectangle will face you ifyou're looking down the positive axis
    int tex_width=0,tex_height=0; // texture width and height
    switch (slice->axis){
        case 0:
            opengl_textured_rect->frame->t.x=grid.X(TV_INT()+slice->index)(slice->axis);
            opengl_textured_rect->frame->r=ROTATION<TV>(0.5*pi,TV(0,1,0));
            opengl_textured_rect->width=domain.Edge_Lengths().z;
            opengl_textured_rect->height=domain.Edge_Lengths().y;
            tex_width=values.Size().z;
            tex_height=values.Size().y;
            break;

        case 1:
            opengl_textured_rect->frame->t.y=grid.X(TV_INT()+slice->index)(slice->axis);
            opengl_textured_rect->frame->r=ROTATION<TV>(-0.5*pi,TV(1,0,0));
            opengl_textured_rect->width=domain.Edge_Lengths().x;
            opengl_textured_rect->height=domain.Edge_Lengths().z;
            tex_width=values.Size().x;
            tex_height=values.Size().z;
            break;

        case 2:
            opengl_textured_rect->frame->t.z=grid.X(TV_INT()+slice->index)(slice->axis);
            opengl_textured_rect->frame->r=ROTATION<TV>();
            opengl_textured_rect->width=domain.Edge_Lengths().x;
            opengl_textured_rect->height=domain.Edge_Lengths().y;
            tex_width=values.Size().x;
            tex_height=values.Size().y;
            break;}

    Update_Slice_Helper(this,tex_width,tex_height);
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Slice_Has_Changed()
{
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)this->slice;
    if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE) // In 3D mode now, can erase textured rect
        Delete_Textured_Rect();
    Update();
}
//#####################################################################
// Function Update_Points
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Update_Points()
{
    PHYSBAM_ASSERT(opengl_points);
    opengl_points->points.Resize(values.Size().Product());
    int index=0;
    for(int i=values.domain.min_corner.x;i<values.domain.max_corner.x;i++) for(int j=values.domain.min_corner.y;j<values.domain.max_corner.y;j++) for(int k=values.domain.min_corner.z;k<values.domain.max_corner.z;k++){
        opengl_points->points(index)=grid.X(TV_INT(i,j,k));
        opengl_points->Set_Point_Color(index,color_maps(current_color_map)->Lookup(values(i,j,k)));
        index++;}
}
//#####################################################################
// Specialization for bool scalars: only draw point if value is "true"
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<float,bool>::
Update_Points()
{
    PHYSBAM_ASSERT(opengl_points);
    opengl_points->color=color_maps(current_color_map)->Lookup(true);
    opengl_points->points.Resize(values.Size().Product());
    int index=0;
    for(int i=values.domain.min_corner.x;i<values.domain.max_corner.x;i++) for(int j=values.domain.min_corner.y;j<values.domain.max_corner.y;j++) for(int k=values.domain.min_corner.z;k<values.domain.max_corner.z;k++)
        if(values(i,j,k)) opengl_points->points(index++)=grid.X(TV_INT(i,j,k));
    opengl_points->points.Resize(index);
}
//#####################################################################
// Function Update_Points
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<double,bool>::
Update_Points()
{
    PHYSBAM_ASSERT(opengl_points);
    opengl_points->color=color_maps(current_color_map)->Lookup(true);
    opengl_points->points.Resize(values.Size().Product());
    int index=0;
    for(int i=values.domain.min_corner.x;i<values.domain.max_corner.x;i++) for(int j=values.domain.min_corner.y;j<values.domain.max_corner.y;j++) for(int k=values.domain.min_corner.z;k<values.domain.max_corner.z;k++)
        if(values(i,j,k)) opengl_points->points(index++)=grid.X(TV_INT(i,j,k));
    opengl_points->points.Resize(index);
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Toggle_Draw_Mode()
{
    DRAW_MODE new_draw_mode=(DRAW_MODE)(((int)draw_mode+1)%2);
    Set_Draw_Mode(new_draw_mode);
}
//#####################################################################
// Function Set_Smooth_Slice_Texture
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Set_Smooth_Slice_Texture(bool smooth_slice_texture_input)
{
    smooth_slice_texture=smooth_slice_texture_input;
    if(draw_mode==DRAW_TEXTURE && opengl_textured_rect && opengl_textured_rect->texture) opengl_textured_rect->texture->Set_Smooth_Shading(smooth_slice_texture);
}
//#####################################################################
// Function Toggle_Smooth_Slice_Texture
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Toggle_Smooth_Slice_Texture()
{
    Set_Smooth_Slice_Texture(!smooth_slice_texture);
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Toggle_Color_Map()
{
    current_color_map=(current_color_map+1)%color_maps.m;
    Update();
}
template class OPENGL_SCALAR_FIELD_3D<float,int>;
template class OPENGL_SCALAR_FIELD_3D<float,bool>;
template class OPENGL_SCALAR_FIELD_3D<float,float>;
template class OPENGL_SCALAR_FIELD_3D<double,int>;
template class OPENGL_SCALAR_FIELD_3D<double,bool>;
template class OPENGL_SCALAR_FIELD_3D<double,double>;
}
