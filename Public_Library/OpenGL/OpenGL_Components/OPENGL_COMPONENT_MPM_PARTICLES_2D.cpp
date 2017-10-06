//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_MPM_PARTICLES_2D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
OPENGL_COMPONENT_MPM_PARTICLES_2D(STREAM_TYPE stream_type,const std::string &filename_input)
    :OPENGL_COMPONENT<T>(stream_type,"Particles 2D"), color_map(OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map()),
    particles(*new MPM_PARTICLES<TV>),
    default_color(OPENGL_COLOR::Yellow()),velocity_color(OPENGL_COLOR(1,(T).078,(T).576)),
    draw_velocities(false),draw_phases(false),draw_arrows(true),draw_B(false),draw_F(false),
    B_color(OPENGL_COLOR::Red(),OPENGL_COLOR::Green()),
    F_color(OPENGL_COLOR::Red(),OPENGL_COLOR::Green()),scale_velocities((T).025),
    filename(filename_input),frame_loaded(-1),valid(false),selected_index(-1)
{
    // We use white color for the default phase(0).
    color_map->Set_Color(0,OPENGL_COLOR(1, 1, 1));

    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("toggle_draw_phases",{[this](){Toggle_Draw_Phases();},"Toggle draw phases"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrow heads"});
    viewer_callbacks.Set("toggle_F",{[this](){draw_F=!draw_F;},"Toggle F display"});
    viewer_callbacks.Set("toggle_B",{[this](){draw_B=!draw_B;},"Toggle B display"});

    is_animation=Is_Animated(filename);
    // Don't need to call Reinitialize here because it will be called in first call to Set_Frame
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
~OPENGL_COMPONENT_MPM_PARTICLES_2D()
{
    delete &particles;
    delete color_map;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Valid_Frame(int frame_input) const
{
    return Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Display() const
{
    if(!valid || !draw) return;
    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}

    GLint mode=0;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
    glEnable(GL_CULL_FACE);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
    glPointSize(5);
    glDisable(GL_LIGHTING);

    glGetIntegerv(GL_RENDER_MODE,&mode);

    ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >("color");
    ARRAY_VIEW<MATRIX<T,TV::m> >* B=particles.template Get_Array<MATRIX<T,TV::m> >("B");
    if(!B) B=particles.template Get_Array<MATRIX<T,TV::m> >("C");

    ARRAY_VIEW<PHASE_ID>* phase=particles.template Get_Array<PHASE_ID>("phase");

    if(draw_velocities && mode!=GL_SELECT){
        glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int i=0;i<particles.X.m;i++){
            TV X=particles.X(i);
            TV Y=X+particles.V(i)*scale_velocities;
            if(draw_arrows) OPENGL_SHAPES::Draw_Arrow(X,Y);
            else OpenGL_Line(X,Y);}
        OpenGL_End();
        glPopAttrib();}

    if(draw_F && mode!=GL_SELECT){
        glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        OpenGL_Begin(GL_LINES);
        for(int i=0;i<particles.F.m;i++){
            TV X=particles.X(i);
            MATRIX<T,TV::m> F=particles.F(i);
            for(int a=0;a<TV::m;a++){
                F_color(a).Send_To_GL_Pipeline();
                OpenGL_Line(X,X+scale_velocities*F.Column(a));}}
        OpenGL_End();
        glPopAttrib();}

    if(draw_B && B && mode!=GL_SELECT){
        glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        OpenGL_Begin(GL_LINES);
        for(int i=0;i<B->m;i++){
            TV X=particles.X(i);
            for(int a=0;a<TV::m;a++){
                B_color(a).Send_To_GL_Pipeline();
                OpenGL_Line(X,X+scale_velocities*(*B)(i).Column(a));}}
        OpenGL_End();
        glPopAttrib();}

    if(mode==GL_SELECT) glPushName(0);
    for(int i=0;i<particles.X.m;i++){
        if(mode==GL_SELECT) glLoadName(i);

        if(draw_phases && phase) color_map->Lookup(Value((*phase)(i))).Send_To_GL_Pipeline();
        else if(colors) OPENGL_COLOR((*colors)(i)).Send_To_GL_Pipeline();
        else default_color.Send_To_GL_Pipeline();

        OpenGL_Begin(GL_POINTS);
        OpenGL_Vertex(particles.X(i));
        OpenGL_End();}
    if(mode==GL_SELECT) glPopName();

    if(mode!=GL_SELECT && selected_index>=0)
        OPENGL_SELECTION::Draw_Highlighted_Vertex(particles.X(selected_index),selected_index);
    
    glPopAttrib();
    glPopMatrix();
    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Use_Bounding_Box() const
{
    return draw && valid && particles.X.Size()>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(RANGE<TV>::Bounding_Box(particles.X));
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    return 110;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    selected_index=indices(0);
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Clear_Selection()
{
    selected_index=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    if(selected_index>=0){
        output_stream<<"MPM particle"<<std::endl;
        output_stream<<"(total number = "<<particles.Size()<<")"<<std::endl;
        output_stream<<"current index = "<<selected_index<<std::endl;
        particles.Print(output_stream,selected_index);}
}
//#####################################################################
// Function Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Destroy_Selection_After_Frame_Change()
{
    return false;
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Selection_Bounding_Box() const
{
    return World_Space_Box(RANGE<TV>(particles.X(selected_index)));
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))) return;
    valid=true;

    std::string frame_filename;
    frame_filename=Get_Frame_Filename(filename,frame);
        
    try{
        std::istream* input_file=Safe_Open_Input(frame_filename);
        TYPED_ISTREAM typed_input(*input_file,stream_type);
        Read_Binary(typed_input,particles);
        delete input_file;}
    catch(FILESYSTEM_ERROR&){valid=false;}
    frame_loaded=frame;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
}
//#####################################################################
// Function Toggle_Draw_Phases
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Toggle_Draw_Phases()
{
    draw_phases=!draw_phases;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Increase_Vector_Size()
{
    scale_velocities*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Decrease_Vector_Size()
{
    scale_velocities/=(T)1.1;
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Toggle_Arrowhead()
{
    draw_arrows=!draw_arrows;
}
//#####################################################################
// Function Set_Slice
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_2D<T>::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;
    Slice_Has_Changed();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_MPM_PARTICLES_2D<float>;
template class OPENGL_COMPONENT_MPM_PARTICLES_2D<double>;
}
