//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <OpenGL/OpenGL/OPENGL_POINTS_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T_ARRAY> OPENGL_POINTS_3D<T,T_ARRAY>::
OPENGL_POINTS_3D(STREAM_TYPE stream_type,T_ARRAY& points_input,const OPENGL_COLOR& color_input,const T point_size_input)
    :OPENGL_OBJECT<T>(stream_type),points(points_input),color(color_input),point_size(point_size_input),draw_point_numbers(false),point_colors(0),point_ids(0),draw_mask(0),selected_index(-1),selected_id(-1),selected_old_color(color)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T_ARRAY> OPENGL_POINTS_3D<T,T_ARRAY>::
~OPENGL_POINTS_3D()
{
    delete point_colors;
    delete point_ids;
    delete draw_mask;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T_ARRAY> RANGE<VECTOR<T,3> > OPENGL_POINTS_3D<T,T_ARRAY>::
Bounding_Box() const
{
    if(!points.Size()) return RANGE<VECTOR<T,3> >::Empty_Box();
    return World_Space_Box(RANGE<TV>::Bounding_Box(points));
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Display() const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
    glPointSize(point_size);
    glDisable(GL_LIGHTING);

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE,&mode);
    OPENGL_COLOR current_color;
    glGetFloatv(GL_CURRENT_COLOR,current_color.rgba);
    if(mode==GL_SELECT){
        glPushName(0);
        color.Send_To_GL_Pipeline();
        for(int i=0;i<points.Size();i++)if(!draw_mask || (*draw_mask)(i)){
            glLoadName(i);
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(points(i));
            OpenGL_End();}
        glPopName();}
    else
    {
        color.Send_To_GL_Pipeline(); // set color outside OpenGL_Begin/OpenGL_End to work around apparent driver bug
        OpenGL_Begin(GL_POINTS);
        if(draw_mask){
            if(point_colors){for(int i=0;i<points.Size();i++)if((*draw_mask)(i)){(*point_colors)(i).Send_To_GL_Pipeline();OpenGL_Vertex(points(i));}}
            else{color.Send_To_GL_Pipeline();for(int i=0;i<points.Size();i++)if((*draw_mask)(i)) OpenGL_Vertex(points(i));}}
        else if(point_colors) for(int i=0;i<points.Size();i++){(*point_colors)(i).Send_To_GL_Pipeline();OpenGL_Vertex(points(i));}
        else{color.Send_To_GL_Pipeline();for(int i=0;i<points.Size();i++) OpenGL_Vertex(points(i));}
        OpenGL_End();
        if(draw_point_numbers)
            for(int i=0;i<points.Size();i++)if(!draw_mask || (*draw_mask)(i)){
                OPENGL_COLOR label_color=(point_colors)?((*point_colors)(i)*0.8):(color*0.8);
                label_color.Send_To_GL_Pipeline();
                OpenGL_String(points(i),point_ids?LOG::sprintf("%d [id=%d] [%g %g]",i,(*point_ids)(i),points(i).x,points(i).y):LOG::sprintf("%d",i));}
        }
    current_color.Send_To_GL_Pipeline();

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Set_Point_Color
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Set_Point_Color(int index,const OPENGL_COLOR &point_color)
{
    Store_Point_Colors(true);
    (*point_colors)(index)=point_color;
}
//#####################################################################
// Function Set_Point_Colors
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Set_Point_Colors(const ARRAY<int> &indices, const OPENGL_COLOR &point_color)
{
    Store_Point_Colors(true);
    for(int i=0;i<indices.m;i++) (*point_colors)(indices(i))=point_color;
}
//#####################################################################
// Function Reset_Point_Colors
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Reset_Point_Colors()
{
    if(point_colors) point_colors->Fill(color);
}    
//#####################################################################
// Function Store_Point_Colors
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Store_Point_Colors(const bool store_point_colors)
{
    if(store_point_colors){
        if(!point_colors){
            point_colors=new ARRAY<OPENGL_COLOR>(points.Size(),false);
            point_colors->Fill(color);}
        else if(point_colors->m!=points.Size()){
            point_colors->Resize(points.Size());
            point_colors->Fill(color);}}
    else{delete point_colors;point_colors=0;}
}
//#####################################################################
// Function Store_Point_Ids
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Store_Point_Ids(bool store_ids)
{
    if(store_ids){
        if(!point_ids) point_ids=new ARRAY<int>();
        point_ids->Resize(points.Size());}
    else{delete point_ids;point_ids=0;}
}
//#####################################################################
// Function Store_Draw_Mask
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Store_Draw_Mask(bool store_draw_mask)
{
    if(store_draw_mask){if(!draw_mask) draw_mask=new ARRAY<bool>();draw_mask->Resize(points.Size());}
    else{delete draw_mask;draw_mask=0;}
}
//#####################################################################
// Function Filter_By_Slice_Planes
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Filter_By_Slice_Planes(const PLANE<T> &leftplane, const PLANE<T> &rightplane)
{
    Store_Draw_Mask(true);
    for(int i=0;i<points.Size();i++) (*draw_mask)(i)=leftplane.Lazy_Inside(points(i)) && rightplane.Lazy_Inside(points(i));
}
//#####################################################################
// Function Clear_Filter
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Clear_Filter()
{
    Store_Draw_Mask(false);
}
//#####################################################################
// Function Select_Point
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Select_Point(int index)
{
    Set_Point_Color(index,OPENGL_COLOR::Yellow());
}
//#####################################################################
// Function Select_Points
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Select_Points(const ARRAY<int> &indices)
{
    Set_Point_Colors(indices,OPENGL_COLOR::Yellow());
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T,class T_ARRAY> int OPENGL_POINTS_3D<T,T_ARRAY>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    PHYSBAM_ASSERT(indices.m==1);
    return 100;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class T_ARRAY> bool OPENGL_POINTS_3D<T,T_ARRAY>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    selected_index=indices(0);
    if(point_ids) selected_id=(*point_ids)(indices(0));
    Store_Point_Colors(true);
    selected_old_color=(*point_colors)(selected_index);
    (*point_colors)(selected_index)=OPENGL_COLOR::Yellow();
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Clear_Selection()
{
    (*point_colors)(selected_index)=selected_old_color;
    selected_index=-1;
    selected_id=-1;
    selected_old_color=color;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T_ARRAY> RANGE<VECTOR<T,3> > OPENGL_POINTS_3D<T,T_ARRAY>::
Selection_Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<T,3> >(points(selected_index)));
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Print_Selection_Info(std::ostream &output_stream) const
{
    output_stream<<"Free particle "<<selected_index<<std::endl;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_POINTS_3D<float,ARRAY<VECTOR<float,3> > >;
template class OPENGL_POINTS_3D<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,3> > > >;
template class OPENGL_POINTS_3D<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,3> > > >;
template class OPENGL_POINTS_3D<double,ARRAY<VECTOR<double,3> > >;
template class OPENGL_POINTS_3D<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,3> > > >;
template class OPENGL_POINTS_3D<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,3> > > >;
}
