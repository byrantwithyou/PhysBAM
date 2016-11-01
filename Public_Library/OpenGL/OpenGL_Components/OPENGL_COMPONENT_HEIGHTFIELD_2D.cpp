//#####################################################################
// Copyright 2004-2016, Eran Guendelman, craig.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_HEIGHTFIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
OPENGL_COMPONENT_HEIGHTFIELD_2D(STREAM_TYPE stream_type,const GRID<TV2> &grid_input, 
                                const std::string& height_filename_input,
                                const std::string& xz_filename_input,
                                const std::string& uv_filename_input,
                                int m_start_input,int m_end_input,int n_start_input,int n_end_input)
    :OPENGL_COMPONENT<T>(stream_type,"Heightfield 2D"), 
    triangulated_surface(*TRIANGULATED_SURFACE<T>::Create()),
    opengl_triangulated_surface(stream_type,triangulated_surface,false,OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Cyan())),
    vertical_offset(0), allow_smooth_shading(true), subdivide_surface(false),
    initial_grid(grid_input), grid(grid_input), height(grid.Domain_Indices()), 
    opengl_vector_field(stream_type,vector_field, vector_locations, OPENGL_COLOR::Green(), 0.025, true, false, true),
    grid_filename(""), scale(1), displacement_scale(1), valid(false), draw_velocities(true), use_triangle_strip(false)
{
    viewer_callbacks.Set("increase_scale",{[this](){Increase_Scale();},"Increase scale"});
    viewer_callbacks.Set("decrease_scale",{[this](){Decrease_Scale();},"Decrease scale"});
    viewer_callbacks.Set("increase_displacement_scale",{[this](){Increase_Displacement_Scale();},"Increase displacement scale"});
    viewer_callbacks.Set("decrease_displacement_scale",{[this](){Decrease_Displacement_Scale();},"Decrease displacement scale"});
    viewer_callbacks.Set("increase_velocity_scale",{[this](){Increase_Velocity_Scale();},"Increase velocity scale"});
    viewer_callbacks.Set("decrease_velocity_scale",{[this](){Decrease_Velocity_Scale();},"Decrease velocity scale"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("toggle_subdivision",{[this](){Toggle_Subdivision();},"Toggle subdivision"});

    if(m_start_input == 0 && m_end_input == 0) {
        domain=grid.Domain_Indices();
    } else {
        domain.min_corner.x = m_start_input; domain.max_corner.x = m_end_input; domain.min_corner.y = n_start_input; domain.max_corner.y = n_end_input;
    }
    counts=domain.Edge_Lengths();
    PHYSBAM_ASSERT(counts.Min()>0);

    opengl_triangulated_surface.Set_Two_Sided();
    opengl_triangulated_surface.Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)0,0.6,(float)0.9)));
    opengl_triangulated_surface.Set_Back_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)0.8,(float)0,0)));

    triangulated_surface.mesh.Initialize_Square_Mesh(counts.x,counts.y);
    for(int i=domain.min_corner.x;i<domain.max_corner.x;i++) for(int j=domain.min_corner.y;j<domain.max_corner.y;j++)
    {
        triangulated_surface.particles.Add_Element();
    }

    height_filename=height_filename_input;
    if(xz_filename_input.length()){xz=new ARRAY<TV2,TV_INT2>;xz_filename=xz_filename_input;}
    else{xz=0;xz_filename="";}
    if(uv_filename_input.length()){
        uv = new ARRAY<TV2,TV_INT2>;
        vector_field.Resize(counts.Product());
        vector_locations.Resize(counts.Product());
        uv_filename=uv_filename_input;}
    else{uv=0;uv_filename="";}

    is_animation = FILE_UTILITIES::Is_Animated(height_filename);
    frame_loaded = -1;

    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
~OPENGL_COMPONENT_HEIGHTFIELD_2D()
{
    delete &triangulated_surface.mesh;
    delete &triangulated_surface.particles;
    delete &triangulated_surface;
    delete xz;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(is_animation?LOG::sprintf(height_filename.c_str(),frame_input):height_filename);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Display() const
{
    if(valid && draw)
    {
        if(use_triangle_strip)
        {
            opengl_triangulated_surface.front_material.Send_To_GL_Pipeline();
            TV v1,v2,v3,v4;
            for(int i = domain.min_corner.x; i <= domain.max_corner.x-1; i++)
            {
                OpenGL_Begin(GL_TRIANGLE_STRIP);
                TV2 lo=grid.X(TV_INT2(i,1)),hi=grid.X(TV_INT2(i+1,1));
                v1=TV(hi.x, scale*(height(i+1,0)+vertical_offset), hi.y);
                v2=TV(lo.x, scale*(height(i,0)+vertical_offset), lo.y);
                OpenGL_Vertex(v1);
                OpenGL_Vertex(v2);
                for(int j = domain.min_corner.y+1; j <= domain.max_corner.y; j++)
                {
                    TV2 nx=grid.X(TV_INT2(i,j));
                    v3=TV(hi.x, scale*(height(i+1,j)+vertical_offset), nx.y);
                    v4=TV(lo.x, scale*(height(i,j)+vertical_offset), nx.y);

                    OpenGL_Normal(PLANE<T>(v1,v2,v3).Normal());
                    OpenGL_Vertex(v3);
                    OpenGL_Normal(PLANE<T>(v3,v2,v4).Normal());
                    OpenGL_Vertex(v4);
                    v1=v3;
                    v2=v4;
                }
                OpenGL_End();
            }
        }
        else
            opengl_triangulated_surface.Display();

        if(draw_velocities) opengl_vector_field.Display();
    }
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> auto OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Bounding_Box() const -> RANGE<TV>
{
    if(valid && draw) return opengl_triangulated_surface.Bounding_Box();
    else return RANGE<TV>::Centered_Box();
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Turn_Smooth_Shading_On()
{
    if(allow_smooth_shading) opengl_triangulated_surface.Turn_Smooth_Shading_On();
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Turn_Smooth_Shading_Off()
{
    opengl_triangulated_surface.Turn_Smooth_Shading_Off();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Reinitialize(bool force)
{
    if(draw){
        if(force || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            bool success = true;
            valid = false;

            if(success && !grid_filename.empty()){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(grid_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File(stream_type,filename,grid);
                    PHYSBAM_ASSERT(grid.counts.x==initial_grid.counts.x && grid.counts.y==initial_grid.counts.y);}
                else success=false;}

            if(success){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(height_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File(stream_type,filename,height);
                    if(!height.Size().All_Greater_Equal(counts)) success = false;}
                else success=false;}

            if(success && xz){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(xz_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File(stream_type,filename,*xz);
                    if(!xz->Size().All_Greater_Equal(counts)) success = false;}
                else success=false;}

            if(success && draw_velocities && uv_filename.length()){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(uv_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File(stream_type,filename,*uv);
                    if(!uv->Size().All_Greater_Equal(counts))
                        success = false;
                    else{
                        int idx = 1;
                        for(int i=domain.min_corner.x;i<domain.max_corner.x;i++)for(int j=domain.min_corner.y;j<domain.max_corner.y;j++){
                            vector_field(idx) = TV((*uv)(i,j).x,0,(*uv)(i,j).y);
                            TV2 pt=grid.X(TV_INT2(i,j));
                            vector_locations(idx) = TV(pt.x, scale*(height(i,j)+vertical_offset), pt.y);
                            idx++;}}}
                else success=false;}

            if(success){
                Update_Surface();
                frame_loaded = frame;
                valid = true;}}}
}
//#####################################################################
// Function Update_Surface
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Update_Surface()
{
    if(use_triangle_strip) return;

    if(xz)
    {
        for(int i = domain.min_corner.x; i <= domain.max_corner.x; i++) for(int j = domain.min_corner.y; j <= domain.max_corner.y; j++)
        {
            TV2 pt=grid.X(TV_INT2(i,j));
            triangulated_surface.particles.X(To_Linear_Index(i,j)) = 
                TV(pt.x + displacement_scale*((*xz)(i,j).x-pt.x), 
                    scale*(height(i,j)+vertical_offset), 
                    pt.y + displacement_scale*((*xz)(i,j).y-pt.y));
        }
    }
    else
    {
        if(subdivide_surface)
        {
            triangulated_surface.mesh.Initialize_Square_Mesh(grid.counts.x,grid.counts.y);
            triangulated_surface.particles.Clean_Memory();
            triangulated_surface.particles.Preallocate(grid.counts.x*grid.counts.y);
            for(int i = domain.min_corner.x; i <= domain.max_corner.x; i++) for(int j = domain.min_corner.y; j <= domain.max_corner.y; j++)
                triangulated_surface.particles.Add_Element();
        }

        for(int i = domain.min_corner.x; i <= domain.max_corner.x; i++) for(int j = domain.min_corner.y; j <= domain.max_corner.y; j++)
        {
            TV2 pt=grid.X(TV_INT2(i,j));
            triangulated_surface.particles.X(To_Linear_Index(i,j)) = 
                TV(pt.x, scale*(height(i,j)+vertical_offset), pt.y);
        }

        if(subdivide_surface) triangulated_surface.Loop_Subdivide();
    }
    if(opengl_triangulated_surface.Is_Smooth_Normals())
        opengl_triangulated_surface.Initialize_Vertex_Normals();
    else
        opengl_triangulated_surface.Delete_Vertex_Normals();
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    return opengl_triangulated_surface.Get_Selection_Priority(indices);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(use_triangle_strip) return false;
    return opengl_triangulated_surface.Set_Selection(indices,modifiers);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Clear_Selection()
{
    opengl_triangulated_surface.Clear_Selection();
}
//#####################################################################
// Function Set_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Set_Scale(T scale_input)
{
    scale = scale_input;
    if(valid) Update_Surface();
}
//#####################################################################
// Function Use_Triangle_Strip
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Use_Triangle_Strip(bool use_triangle_strip_input)
{
    use_triangle_strip=use_triangle_strip_input;
}
//#####################################################################
// Function Increase_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Increase_Scale()
{
    scale *= 1.1;
    if(valid) Update_Surface();
}
//#####################################################################
// Function Decrease_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Decrease_Scale()
{
    scale *= 1/1.1;
    if(valid) Update_Surface();
}
//#####################################################################
// Function Increase_Displacement_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Increase_Displacement_Scale()
{
    displacement_scale *= 1.1;
    if(valid) Update_Surface();
}
//#####################################################################
// Function Decrease_Displacement_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Decrease_Displacement_Scale()
{
    displacement_scale *= 1/1.1;
    if(valid) Update_Surface();
}
//#####################################################################
// Function Increase_Velocity_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Increase_Velocity_Scale()
{
    opengl_vector_field.size *= 1.1;
}
//#####################################################################
// Function Decrease_Velocity_Scale
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Decrease_Velocity_Scale()
{
    opengl_vector_field.size *= 1/1.1;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities = !draw_velocities;
    Reinitialize(true);
}
//#####################################################################
// Function Toggle_Subdivision
//#####################################################################
template<class T> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Toggle_Subdivision()
{
    subdivide_surface = !subdivide_surface;
    Reinitialize(true);
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_HEIGHTFIELD_2D<T>::
Selection_Bounding_Box() const
{
    return opengl_triangulated_surface.Selection_Bounding_Box();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_HEIGHTFIELD_2D<float>;
template class OPENGL_COMPONENT_HEIGHTFIELD_2D<double>;
}
