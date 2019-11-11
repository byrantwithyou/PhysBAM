//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_HEIGHTFIELD_1D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
OPENGL_COMPONENT_HEIGHTFIELD_1D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid_input,
                                const std::string& height_filename_input,
                                const std::string& x_filename_input,
                                const std::string& ground_filename_input,
                                const std::string& u_filename_input)
    :OPENGL_COMPONENT<T>(viewer_dir,"Heightfield 1D"), grid(grid_input), opengl_vector_field(vector_field,vector_locations), scale(1), displacement_scale(1), valid(false), 
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
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Set_Frame()
{
    
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
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Selection_Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<T,3> >(VECTOR<T,3>(grid.X(TV_INT(selected_index)).x, scale*height(selected_index),0)));
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Reinitialize()
{
    if(!draw) return;
    bool success=true;
    valid=false;

    if(success){
        std::string filename=viewer_dir.current_directory+"/"+height_filename;
        if(File_Exists(filename)){
            Read_From_File(filename,height);
            if(height.Size().x != grid.counts.x) success=false;}
        else success=false;}

    if(success && x){
        std::string filename=viewer_dir.current_directory+"/"+x_filename;
        if(File_Exists(filename)){
            Read_From_File(filename,*x);
            if(height.Size().x != x->Size().x) success=false;}
        else success=false;}
    
    if(success && ground){
        std::string filename=viewer_dir.current_directory+"/"+ground_filename;
        if(File_Exists(filename)){
            Read_From_File(filename,*ground);
            if(height.Size().x != ground->Size().x) success=false;
            else height += (*ground);}
        else success=false;}

    if(success && draw_velocities && u_filename.length()){
        std::string filename=viewer_dir.current_directory+"/"+u_filename;
        if(File_Exists(filename)){
            Read_From_File(filename,*u);
            if(height.Size().x != u->Size().x) success=false;
            else for(int i=0;i<grid.counts.x;i++){
                    vector_field(i)=VECTOR<T,2>((*u)(i),0);
                    vector_locations(i)=VECTOR<T,2>(grid.X(TV_INT(i)).x, scale*height(i));}}
        else success=false;}

    if(success) valid=true;
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    return 10;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    selected_index=indices(0);
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Clear_Selection()
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
    Reinitialize();// To recompute velocity vector positions correctly
}
//#####################################################################
// Function Increase_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Increase_Scale()
{
    scale *= 1.1;
    Reinitialize();// To recompute velocity vector positions correctly
}
//#####################################################################
// Function Decrease_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Decrease_Scale()
{
    scale *= 1/1.1;
    Reinitialize();// To recompute velcity vector positions correctly
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
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Points
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T>::
Toggle_Draw_Points()
{
    draw_points=!draw_points;
}

namespace PhysBAM{
template class OPENGL_COMPONENT_HEIGHTFIELD_1D<float>;
template class OPENGL_COMPONENT_HEIGHTFIELD_1D<double>;
}
