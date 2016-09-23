//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_HEIGHTFIELD_1D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
OPENGL_COMPONENT_HEIGHTFIELD_1D(STREAM_TYPE stream_type,const GRID<TV> &grid_input, 
                                const std::string& height_filename_input,
                                const std::string& x_filename_input,
                                const std::string& ground_filename_input,
                                const std::string& u_filename_input)
    :OPENGL_COMPONENT<T>(stream_type,"Heightfield 1D"), grid(grid_input), opengl_vector_field(stream_type,vector_field,vector_locations), scale(1), displacement_scale(1), valid(false), 
      draw_velocities(true), draw_points(true), selected_index(0)
{
    viewer_callbacks.Set("increase_scale",{[this](){Increase_Scale();},"Increase scale"});
    viewer_callbacks.Set("decrease_scale",{[this](){Decrease_Scale();},"Decrease scale"});
    viewer_callbacks.Set("increase_displacement_scale",{[this](){Increase_Displacement_Scale();},"Increase displacement scale"});
    viewer_callbacks.Set("decrease_displacement_scale",{[this](){Decrease_Displacement_Scale();},"Decrease displacement space"});
    viewer_callbacks.Set("increase_velocity_scale",{[this](){Increase_Velocity_Scale();},"Increase velocity scale"});
    viewer_callbacks.Set("decrease_velocity_scale",{[this](){Decrease_Velocity_Scale();},"Decrease velocity scale"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("toggle_draw_points",{[this](){Toggle_Draw_Points();},"Toggle draw points"});

    height_filename=height_filename_input;
    if(x_filename_input.length()){x=new ARRAY<T,TV_INT>;x_filename=x_filename_input;}
    else{x=0;x_filename="";}
    if(ground_filename_input.length()){ground=new ARRAY<T,TV_INT>;ground_filename=ground_filename_input;}
    else{ground=0;ground_filename="";}
    if(u_filename_input.length()){
        u=new ARRAY<T,TV_INT>;
        vector_field.Resize(grid.counts.x);
        vector_locations.Resize(grid.counts.x);
        u_filename=u_filename_input;}
    else{u=0;u_filename="";}
    is_animation=height_filename.find("%d")!=std::string::npos;
    frame_loaded=-1;

    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
~OPENGL_COMPONENT_HEIGHTFIELD_1D()
{
    delete x;
    delete ground;
    delete u;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(is_animation?LOG::sprintf(height_filename.c_str(),frame_input):height_filename);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Display() const
{
    OPENGL_COLOR water_color(0,0,1);
    OPENGL_COLOR ground_color=OPENGL_COLOR::Gray(0.8);
    OPENGL_COLOR displaced_water(1,0,0);
    OPENGL_COLOR points_color(0,1,1);
    OPENGL_COLOR selected_point_color=OPENGL_COLOR::Yellow();
    GLfloat point_size=3.0;

    GLint mode=0;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    if(valid && draw)
    {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);

        if(mode == GL_SELECT)
        {
            glPushName(0);
            glPushAttrib(GL_POINT_BIT);
            glPointSize(point_size);
            points_color.Send_To_GL_Pipeline();
            for(int i=0;i<grid.counts.x;i++)
            {
                glLoadName(i);
                OpenGL_Begin(GL_POINTS);
                OpenGL_Vertex(VECTOR<T,2>(grid.X(TV_INT(i)).x, scale*height(i)));
                OpenGL_End();
            }
            glPopAttrib();
            glPopName();
        }
        else
        {
            if(x)
            {
                displaced_water.Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINE_STRIP);
                if(displacement_scale == 1)
                {
                    for(int i=0;i<x->Size().x;i++)
                        OpenGL_Vertex(VECTOR<T,2>((*x)(i), scale*height(i)));       
                }
                else
                {
                    for(int i=0;i<x->Size().x;i++)
                        OpenGL_Vertex(VECTOR<T,2>(grid.X(TV_INT(i)).x + displacement_scale*((*x)(i)-grid.X(TV_INT(i)).x),
                                scale*height(i)));       
                }
                OpenGL_End();
            }

            // Water
            water_color.Send_To_GL_Pipeline();
            OpenGL_Begin(GL_LINE_STRIP);
            for(int i=0;i<grid.counts.x;i++)
                OpenGL_Vertex(VECTOR<T,2>(grid.X(TV_INT(i)).x, scale*height(i)));
            OpenGL_End();

            if(draw_points)
            {
                glPushAttrib(GL_POINT_BIT);
                glPointSize(point_size);
                points_color.Send_To_GL_Pipeline();
                OpenGL_Begin(GL_POINTS);
                for(int i=0;i<grid.counts.x;i++)
                    OpenGL_Vertex(VECTOR<T,2>(grid.X(TV_INT(i)).x, scale*height(i)));
                OpenGL_End();

                for(int i=0;i<grid.counts.x;i++)
                    OpenGL_String(VECTOR<T,2>(grid.X(TV_INT(i)).x,scale*height(i)),LOG::sprintf("%d",i));

                if(selected_index>=0){
                    selected_point_color.Send_To_GL_Pipeline();
                    int i=selected_index;
                    OpenGL_Begin(GL_POINTS);
                    OpenGL_Vertex(VECTOR<T,2>(grid.X(TV_INT(i)).x, scale*height(i)));
                    OpenGL_End();
                    OpenGL_String(VECTOR<T,2>(grid.X(TV_INT(i)).x,scale*height(i)),LOG::sprintf("%d",i));
                }

                glPopAttrib();
            }

            // Ground
            ground_color.Send_To_GL_Pipeline();
            if(ground)
            {
                OpenGL_Begin(GL_LINE_STRIP);
                for(int i=0;i<grid.counts.x;i++)
                    OpenGL_Vertex(VECTOR<T,2>(grid.X(TV_INT(i)).x, scale*(*ground)(i)));      
                OpenGL_End();
            }
            else
            {
                OpenGL_Begin(GL_LINE_STRIP);
                OpenGL_Vertex(VECTOR<T,2>(grid.X(TV_INT()).x, 0));
                OpenGL_Vertex(VECTOR<T,2>(grid.X(grid.counts-1).x, 0));
                OpenGL_End();
            }
        }

        glEnable(GL_LIGHTING);
        glEnable(GL_DEPTH_TEST);

        if(mode != GL_SELECT && draw_velocities) opengl_vector_field.Display();
    }
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Bounding_Box() const
{
    if(valid && draw)
    {
        T min_height=min((T)0,height.Min());
        T max_height=max((T)0,height.Max());

        T min_x=grid.domain.min_corner.x, max_x=grid.domain.max_corner.x;
        if(x) { 
            min_x=min(min_x,x->Min());
            max_x=max(max_x,x->Max());
        }

        return RANGE<VECTOR<T,3> >(VECTOR<T,3>(min_x,min_height,0),VECTOR<T,3>(max_x,max_height,0));
    }
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Reinitialize(bool force)
{
    if(draw){
        if(force || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            bool success=true;
            valid=false;

            if(success){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(height_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File(stream_type,filename,height);
                    if(height.Size().x != grid.counts.x) success=false;}
                else success=false;}

            if(success && x){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(x_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File(stream_type,filename,*x);
                    if(height.Size().x != x->Size().x) success=false;}
                else success=false;}
    
            if(success && ground){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(ground_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File(stream_type,filename,*ground);
                    if(height.Size().x != ground->Size().x) success=false;
                    else height += (*ground);}
                else success=false;}

            if(success && draw_velocities && u_filename.length()){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(u_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File(stream_type,filename,*u);
                    if(height.Size().x != u->Size().x) success=false;
                    else for(int i=0;i<grid.counts.x;i++){
                            vector_field(i)=VECTOR<T,2>((*u)(i),0);
                            vector_locations(i)=VECTOR<T,2>(grid.X(TV_INT(i)).x, scale*height(i));}}
                else success=false;}

            if(success){
                frame_loaded=frame;
                valid=true;}}}
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    if(buffer_size == 1)
    {
        OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T> *selection=new OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T>(this);
        selection->index=buffer[0];
        return selection;
    }
    else return 0;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    if(selection->type != OPENGL_SELECTION<T>::COMPONENT_HEIGHTFIELD_1D) return;
    OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T>*)selection;
    selected_index=real_selection->index;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Clear_Highlight()
{
    selected_index=-1;
}
//#####################################################################
// Function Set_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Set_Scale(T scale_input)
{
    scale=scale_input;
    Reinitialize(true);// To recompute velocity vector positions correctly
}
//#####################################################################
// Function Increase_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Increase_Scale()
{
    scale *= 1.1;
    Reinitialize(true);// To recompute velocity vector positions correctly
}
//#####################################################################
// Function Decrease_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Decrease_Scale()
{
    scale *= 1/1.1;
    Reinitialize(true);// To recompute velcity vector positions correctly
}
//#####################################################################
// Function Increase_Displacement_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Increase_Displacement_Scale()
{
    displacement_scale *= 1.1;
}
//#####################################################################
// Function Decrease_Displacement_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Decrease_Displacement_Scale()
{
    displacement_scale *= 1/1.1;
}
//#####################################################################
// Function Increase_Velocity_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Increase_Velocity_Scale()
{
    opengl_vector_field.size *= 1.1;
}
//#####################################################################
// Function Decrease_Velocity_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Decrease_Velocity_Scale()
{
    opengl_vector_field.size *= 1/1.1;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
    Reinitialize(true);
}
//#####################################################################
// Function Toggle_Draw_Points
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Toggle_Draw_Points()
{
    draw_points=!draw_points;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<T,3> >::Empty_Box();
}

namespace PhysBAM{
template class OPENGL_COMPONENT_HEIGHTFIELD_1D<float>;
template class OPENGL_COMPONENT_HEIGHTFIELD_1D<double>;
}
