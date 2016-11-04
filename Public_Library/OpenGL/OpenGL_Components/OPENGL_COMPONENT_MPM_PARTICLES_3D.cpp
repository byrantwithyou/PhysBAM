//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_MPM_PARTICLES_3D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
OPENGL_COMPONENT_MPM_PARTICLES_3D(STREAM_TYPE stream_type,const std::string &filename_input)
    :OPENGL_COMPONENT<T>(stream_type,"Particles 3D"),particles(*new MPM_PARTICLES<TV>),
    default_color(OPENGL_COLOR::Yellow()),velocity_color(OPENGL_COLOR(1,(T).078,(T).576)),
    draw_velocities(false),draw_arrows(true),draw_B(false),draw_F(false),
    B_color(OPENGL_COLOR::Red()*.5,OPENGL_COLOR::Green()*.5,OPENGL_COLOR::Blue()*.5),
    F_color(OPENGL_COLOR::Red(),OPENGL_COLOR::Green(),OPENGL_COLOR::Blue()),scale_velocities((T).025),
    filename(filename_input),frame_loaded(-1),valid(false),
    selected_index(-1)
{
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrow heads"});

    is_animation=FILE_UTILITIES::Is_Animated(filename);
    // Don't need to call Reinitialize here because it will be called in first call to Set_Frame
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
~OPENGL_COMPONENT_MPM_PARTICLES_3D()
{
    delete &particles;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
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

    ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    ARRAY_VIEW<T>* sizes=particles.template Get_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);

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

    if(draw_B && particles.store_B && mode!=GL_SELECT){
        glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        OpenGL_Begin(GL_LINES);
        for(int i=0;i<particles.B.m;i++){
            TV X=particles.X(i);
            MATRIX<T,TV::m> B=particles.B(i);
            for(int a=0;a<TV::m;a++){
                B_color(a).Send_To_GL_Pipeline();
                OpenGL_Line(X,X+scale_velocities*B.Column(a));}}
        OpenGL_End();
        glPopAttrib();}

    if(mode==GL_SELECT) glPushName(0);
    for(int i=0;i<particles.X.m;i++){
        if(mode==GL_SELECT) glLoadName(i);

        if(colors) OPENGL_COLOR((*colors)(i)).Send_To_GL_Pipeline();
        else default_color.Send_To_GL_Pipeline();

        if(sizes && (*sizes)(i)) OPENGL_SHAPES::Draw_Circle(particles.X(i),(*sizes)(i)*scale_velocities,20,false);
        else{
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(particles.X(i));
            OpenGL_End();}}
    if(mode==GL_SELECT) glPopName();

    glPopAttrib();
    glPopMatrix();
    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Use_Bounding_Box() const
{
    return draw && valid && particles.X.Size()>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(RANGE<TV>::Bounding_Box(particles.X));
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    return 110;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    selected_index=indices(0);
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Clear_Selection()
{
    selected_index=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    if(selected_index>=0){
        output_stream<<"(total number = "<<particles.Size()<<")"<<std::endl;
        output_stream<<"current index = "<<selected_index<<std::endl;
        particles.Print(output_stream,selected_index);}
}
//#####################################################################
// Function Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Destroy_Selection_After_Frame_Change()
{
    return true;
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Selection_Bounding_Box() const
{
    return World_Space_Box(RANGE<TV>(particles.X(selected_index)));
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))) return;
    valid=true;

    std::string frame_filename;
    frame_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
        
    try{
        std::istream* input_file=FILE_UTILITIES::Safe_Open_Input(frame_filename);
        TYPED_ISTREAM typed_input(*input_file,stream_type);
        Read_Binary(typed_input,particles);
        delete input_file;}
    catch(FILESYSTEM_ERROR&){valid=false;}
    frame_loaded=frame;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Increase_Vector_Size()
{
    scale_velocities*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Decrease_Vector_Size()
{
    scale_velocities/=(T)1.1;
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Toggle_Arrowhead()
{
    draw_arrows=!draw_arrows;
}
//#####################################################################
// Function Set_Slice
//#####################################################################
template<class T> void OPENGL_COMPONENT_MPM_PARTICLES_3D<T>::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;
    Slice_Has_Changed();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_MPM_PARTICLES_3D<float>;
template class OPENGL_COMPONENT_MPM_PARTICLES_3D<double>;
}
