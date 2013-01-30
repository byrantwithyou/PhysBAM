//#####################################################################
// Copyright 2004, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_LEVELSET_2D<T>::
OPENGL_LEVELSET_2D(LEVELSET<TV>& levelset_input,const OPENGL_COLOR& inside_color_input,const OPENGL_COLOR& outside_color_input,ARRAY<bool,VECTOR<int,2> > *active_cells_input)
    :OPENGL_SCALAR_FIELD_2D<T,T>(levelset_input.grid,levelset_input.phi,OPENGL_COLOR_RAMP<T>::Levelset_Color_Constant_Ramp(inside_color_input,outside_color_input),active_cells_input),
    levelset(levelset_input),opengl_triangulated_area(0),opengl_segmented_curve_2d(0),color_mode(COLOR_SOLID),inside_color(inside_color_input),outside_color(outside_color_input), 
    gradient_color_map(0),draw_cells(false),draw_area(false),draw_curve(true),dominant_sign(-1),draw_normals(false)
{
    Update();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_LEVELSET_2D<T>::
~OPENGL_LEVELSET_2D()
{
    delete opengl_segmented_curve_2d;
    delete opengl_triangulated_area;
    delete solid_color_map;
    delete gradient_color_map;
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_LEVELSET_2D<T>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    if(draw_cells) OPENGL_SCALAR_FIELD_2D<T,T>::Display(in_color);
    if(draw_normals){
        levelset.Compute_Normals();
        glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        OPENGL_COLOR::White().Send_To_GL_Pipeline();
        ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
        for(int i=0;i<levelset.grid.counts.x;i++)
            for(int j=0;j<levelset.grid.counts.y;j++)
                if(!active_cells || (*active_cells)(i,j))
                    OpenGL_Line(grid.X(TV_INT(i,j)),grid.X(TV_INT(i,j))+T(0.01)*(*levelset.normals)(i,j),vertices);
        OpenGL_Draw_Arrays(GL_LINES,2,vertices);
        glPopAttrib();}
    if(draw_area && opengl_triangulated_area){
        glDepthMask(GL_FALSE);
        opengl_triangulated_area->Display(in_color);
        glDepthMask(GL_TRUE);}
    if(draw_curve && opengl_segmented_curve_2d) opengl_segmented_curve_2d->Display(in_color);
    glPopMatrix();
}
//#####################################################################
// Function Set_Inside_And_Outside_Colors
//#####################################################################
template<class T> void OPENGL_LEVELSET_2D<T>::
Set_Inside_And_Outside_Colors(const OPENGL_COLOR& inside_color_input,const OPENGL_COLOR& outside_color_input)
{
    inside_color=inside_color_input;
    outside_color=outside_color_input;
    this->color_maps(0);
    this->color_maps(0)=OPENGL_COLOR_RAMP<T>::Levelset_Color_Constant_Ramp(inside_color,outside_color);
}
//#####################################################################
// Function Toggle_Normals
//#####################################################################
template<class T> void OPENGL_LEVELSET_2D<T>::
Toggle_Normals()
{
    if(draw_normals)draw_normals=false;else draw_normals=true;
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_LEVELSET_2D<T>::
Update()
{
    OPENGL_SCALAR_FIELD_2D<T,T>::Update();
    if(levelset.phi.Size().Min()>1){
        if(opengl_triangulated_area){
            delete &opengl_triangulated_area->triangulated_area;
            delete opengl_triangulated_area;
            opengl_triangulated_area=0;}
        if(opengl_segmented_curve_2d){
            delete &opengl_segmented_curve_2d->curve;
            delete opengl_segmented_curve_2d;
            opengl_segmented_curve_2d=0;}
        if(draw_area){
            TRIANGULATED_AREA<T>& ta=*TRIANGULATED_AREA<T>::Create();
            MARCHING_CUBES<TV>::Create_Interior(ta,levelset.grid,levelset.phi,true);
            opengl_triangulated_area=new OPENGL_TRIANGULATED_AREA<T>(ta,false,OPENGL_COLOR::Red(),inside_color,inside_color);}
        if(draw_curve){
            SEGMENTED_CURVE_2D<T>& sc=*SEGMENTED_CURVE_2D<T>::Create();
            MARCHING_CUBES<TV>::Create_Surface(sc,levelset.grid,levelset.phi);
            opengl_segmented_curve_2d=new OPENGL_SEGMENTED_CURVE_2D<T>(sc);}}
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_LEVELSET_2D<float>;
template class OPENGL_LEVELSET_2D<double>;
}
