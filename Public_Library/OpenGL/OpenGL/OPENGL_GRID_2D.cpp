//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Geoffrey Irving, Sndrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
using namespace PhysBAM;
template<class T> OPENGL_GRID_2D<T>::
OPENGL_GRID_2D(STREAM_TYPE stream_type,GRID<TV> &grid_input,const OPENGL_COLOR &color_input,
    const std::string basedir_input,const int frame_input)
    :OPENGL_OBJECT<T>(stream_type),grid(grid_input),active_cell_mask(0),ghost_cell_mask(0),
    active_face_mask(0),ghost_face_mask(0),active_node_mask(0),ghost_node_mask(0),color(color_input),
    draw(true),draw_ghost_values(true),draw_mask_type(0),select_type(SELECT_TYPE::NONE),selected_cell(-1,-1),selected_node(-1,-1),
    basedir(basedir_input),frame(frame_input)
{
    viewer_callbacks.Set("toggle_draw_ghost_values",{[this](){Toggle_Draw_Ghost_Values();},"toggle_draw_ghost_values"});
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Display() const
{
    if(!draw)return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    color.Send_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    int ghost_cells=draw_ghost_values?3:0;

    if(mode == GL_SELECT){
        glPushAttrib(GL_ENABLE_BIT|GL_POINT_BIT);glDisable(GL_CULL_FACE);
        TV min_corner=grid.Node(TV_INT(-ghost_cells,-ghost_cells)),X;int i,j;

        // Draw grid cells for selection
        glPushName(0);
        for(i=-ghost_cells,X.x=min_corner.x;i<grid.numbers_of_cells.x+ghost_cells;i++,X.x+=grid.dX.x){
            glPushName(i);
            for(j=-ghost_cells,X.y=min_corner.y;j<grid.numbers_of_cells.y+ghost_cells;j++,X.y+=grid.dX.y){
                glPushName(j);
                OpenGL_Begin(GL_QUADS);
                OpenGL_Quad_2D(X,X+grid.dX);
                OpenGL_End();
                glPopName();}
            glPopName();}

        // Draw grid nodes for selection
        glLoadName(1);
        glPointSize(OPENGL_PREFERENCES::selection_point_size);
        for(i=-ghost_cells,X.x=min_corner.x;i<grid.numbers_of_cells.x+ghost_cells+1;i++,X.x+=grid.dX.x){
            glPushName(i);
            for(j=-ghost_cells,X.y=min_corner.y;j<grid.numbers_of_cells.y+ghost_cells+1;j++,X.y+=grid.dX.y){
                glPushName(j);
                OpenGL_Begin(GL_POINTS);
                OpenGL_Vertex(X);
                OpenGL_End();
                glPopName();}
            glPopName();}
        glPopName();glPopAttrib();}
    else
    {
        // Draw masks
        T x,y;int i,j,i_mask=1,j_mask=1;
        TV min_corner=grid.Node(TV_INT()-ghost_cells);
        TV max_corner=grid.Node(grid.numbers_of_cells+ghost_cells);

        ARRAY<bool,TV_INT> *cell_mask=0;
        ARRAY<bool,FACE_INDEX<TV::m> > *face_mask=0;
        ARRAY<bool,TV_INT> *node_mask=0;
        if(draw_mask_type==0) cell_mask=active_cell_mask;
        else if(draw_mask_type==1){cell_mask=ghost_cell_mask;ghost_cells=ghost_cells?4:0;}
        else if(draw_mask_type==2) face_mask=active_face_mask;
        else if(draw_mask_type==3){face_mask=ghost_face_mask;ghost_cells=ghost_cells?4:0;}
        else if(draw_mask_type==4) node_mask=active_node_mask;
        else if(draw_mask_type==5){node_mask=ghost_node_mask;ghost_cells=ghost_cells?4:0;}

        if(ghost_cells==4){min_corner=min_corner-grid.dX;max_corner=max_corner+grid.dX;}
        if(cell_mask){
            int mask_m_start=cell_mask->domain.min_corner.x,mask_n_start=cell_mask->domain.min_corner.y;
            glPushAttrib(GL_ALL_ATTRIB_BITS);glLineWidth(4*OPENGL_PREFERENCES::line_width);OPENGL_COLOR::Cyan().Send_To_GL_Pipeline();
            OpenGL_Begin(GL_LINES);
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells;j++,y+=grid.dX.y,j_mask++)
                    if((*cell_mask)(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(TV(x,y),TV(x+grid.dX.x,y+grid.dX.y));}
            OpenGL_End();
            glPopAttrib();}

        if(face_mask){
            int mask_m_start,mask_n_start;
            glPushAttrib(GL_ALL_ATTRIB_BITS);glLineWidth(4*OPENGL_PREFERENCES::line_width);OPENGL_COLOR::Cyan().Send_To_GL_Pipeline();
            OpenGL_Begin(GL_LINES);
            mask_m_start=face_mask->Component(0).domain.min_corner.x,mask_n_start=face_mask->Component(0).domain.min_corner.y;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells+1;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells;j++,y+=grid.dX.y,j_mask++)
                    if((face_mask->Component(0))(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(TV(x,y),TV(x,y+grid.dX.y));}
            mask_m_start=face_mask->Component(1).domain.min_corner.x,mask_n_start=face_mask->Component(1).domain.min_corner.y;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells+1;j++,y+=grid.dX.y,j_mask++)
                    if((face_mask->Component(1))(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(TV(x,y),TV(x+grid.dX.x,y));}
            OpenGL_End();
            glPopAttrib();}

        if(node_mask){
            int mask_m_start=node_mask->domain.min_corner.x,mask_n_start=node_mask->domain.min_corner.y;
            glPushAttrib(GL_ALL_ATTRIB_BITS);glPointSize(5);OPENGL_COLOR::Cyan().Send_To_GL_Pipeline();
            OpenGL_Begin(GL_POINTS);
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells+1;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells+1;j++,y+=grid.dX.y,j_mask++)
                    if((*node_mask)(mask_m_start+i_mask-1,mask_n_start+j_mask-1)) OpenGL_Vertex(TV(x,y));
            OpenGL_End();
            glPopAttrib();}

        // Draw grid
        if(active_face_mask){
            if(draw_mask_type&&ghost_cells==4){ghost_cells=3;min_corner=min_corner+TV(grid.dX.x,grid.dX.y);max_corner=max_corner-TV(grid.dX.x,grid.dX.y);}
            face_mask=active_face_mask;
            int mask_m_start,mask_n_start;
            OpenGL_Begin(GL_LINES);
            mask_m_start=face_mask->Component(0).domain.min_corner.x,mask_n_start=face_mask->Component(0).domain.min_corner.y;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells+1;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells;j++,y+=grid.dX.y,j_mask++)
                    if((face_mask->Component(0))(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(TV(x,y),TV(x,y+grid.dX.y));}
            mask_m_start=face_mask->Component(1).domain.min_corner.x,mask_n_start=face_mask->Component(1).domain.min_corner.y;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells+1;j++,y+=grid.dX.y,j_mask++)
                    if((face_mask->Component(1))(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(TV(x,y),TV(x+grid.dX.x,y));}
            OpenGL_End();}
        else{
            if(draw_mask_type&&ghost_cells==3){ghost_cells=4;min_corner=min_corner-TV(grid.dX.x,grid.dX.y);max_corner=max_corner+TV(grid.dX.x,grid.dX.y);}
            OpenGL_Begin(GL_LINES);
            for(i=-ghost_cells,x=min_corner.x;i<grid.numbers_of_cells.x+ghost_cells+1;i++,x+=grid.dX.x){
                OpenGL_Line(TV(x,min_corner.y),TV(x,max_corner.y));}
            for(j=-ghost_cells,y=min_corner.y;j<grid.numbers_of_cells.y+ghost_cells+1;j++,y+=grid.dX.y){
                OpenGL_Line(TV(min_corner.x,y),TV(max_corner.x,y));}
            OpenGL_End();}

        // Outline boundary of real domain in wider line
        if(ghost_cells>0){
            glPushAttrib(GL_LINE_BIT);glLineWidth(2*OPENGL_PREFERENCES::line_width);
            OpenGL_Begin(GL_LINE_LOOP);
            OpenGL_Quad_2D(grid.domain.min_corner,grid.domain.max_corner);
            OpenGL_End();
            glPopAttrib();}

        // Highlight current selection
        if(select_type==SELECT_TYPE::CELL){
            TV min_corner=grid.Node(selected_cell),max_corner=min_corner+grid.dX;
            OPENGL_SELECTION::Draw_Highlighted_Quad(min_corner,max_corner);}
        if(select_type==SELECT_TYPE::NODE){
            OPENGL_SELECTION::Draw_Highlighted_Vertex(grid.Node(selected_node));}}

    glPopAttrib();
    glPopMatrix();
}
template<class T> void OPENGL_GRID_2D<T>::
Set_Frame(int frame_input)
{
    frame=frame_input;
    Reinitialize();
    return;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_GRID_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_GRID_2D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    PHYSBAM_ASSERT(indices.m==3);
    const static int priority[]={60,70};
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_GRID_2D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    ARRAY<int> signed_indices(indices.m);
    for(int i=1;i<indices.m;i++)
        if(indices(i)&(GLuint)(-1))
            signed_indices(i)=-~(indices(i)-1);
        else
            signed_indices(i)=indices(i);
    if(indices(0)==0){
        select_type=SELECT_TYPE::CELL;
        selected_cell=TV_INT(signed_indices(1),signed_indices(2));}
    else if(indices(0)==1){
        select_type=SELECT_TYPE::NODE;
        selected_node=TV_INT(signed_indices(1),signed_indices(2));}
    else{
        select_type=SELECT_TYPE::NONE;
        return false;}
    return true;
}
//#####################################################################
// Function Clear_Selection
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Clear_Selection()
{
    select_type=SELECT_TYPE::NONE;
}
//#####################################################################
// Function Toggle_Draw_Ghost_Values
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Toggle_Draw_Ghost_Values()
{
    draw_ghost_values=!draw_ghost_values;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Reinitialize()
{
    std::string filename=LOG::sprintf("%s/%d/active_cell_mask",basedir.c_str(),frame);
    if(File_Exists(filename)){
        if(!active_cell_mask) active_cell_mask=new ARRAY<bool,TV_INT>();
        active_cell_mask->Clean_Memory();
        Read_From_File<bool>(filename,*active_cell_mask);}

    filename=LOG::sprintf("%s/%d/ghost_cell_mask",basedir.c_str(),frame);
    if(File_Exists(filename)){
        if(!ghost_cell_mask) ghost_cell_mask=new ARRAY<bool,TV_INT>();
        ghost_cell_mask->Clean_Memory();
        Read_From_File<bool>(filename,*ghost_cell_mask);}

    filename=LOG::sprintf("%s/%d/active_face_mask",basedir.c_str(),frame);
    if(File_Exists(filename)){
        if(!active_face_mask) active_face_mask=new ARRAY<bool,FACE_INDEX<TV::m> >();
        active_face_mask->Clean_Memory();
        Read_From_File<bool>(filename,*active_face_mask);}

    filename=LOG::sprintf("%s/%d/ghost_face_mask",basedir.c_str(),frame);
    if(File_Exists(filename)){
        if(!ghost_face_mask) ghost_face_mask=new ARRAY<bool,FACE_INDEX<TV::m> >();
        ghost_face_mask->Clean_Memory();
        Read_From_File<bool>(filename,*ghost_face_mask);}

    filename=LOG::sprintf("%s/%d/active_node_mask",basedir.c_str(),frame);
    if(File_Exists(filename)){
        if(!active_node_mask) active_node_mask=new ARRAY<bool,TV_INT>();
        active_node_mask->Clean_Memory();
        Read_From_File<bool>(filename,*active_node_mask);}

    filename=LOG::sprintf("%s/%d/ghost_node_mask",basedir.c_str(),frame);
    if(File_Exists(filename)){
        if(!ghost_node_mask) ghost_node_mask=new ARRAY<bool,TV_INT>();
        ghost_node_mask->Clean_Memory();
        Read_From_File<bool>(filename,*ghost_node_mask);}
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Print_Selection_Info(std::ostream& stream) const
{
    if(select_type==SELECT_TYPE::NODE){
        stream<<"Selected node "<<selected_node<<" ("<<grid.Get_Regular_Grid().X(selected_node)<<")"<<std::endl;
        for(int i=0;i<grid_objects.m;i++)
            grid_objects(i)->Print_Node_Selection_Info(stream,selected_node);}
    else if(select_type==SELECT_TYPE::CELL){
        stream<<"Selected cell "<<selected_cell<<" ("<<grid.Get_MAC_Grid().X(selected_cell)<<")"<<std::endl;
        for(int i=0;i<grid_objects.m;i++)
            grid_objects(i)->Print_Cell_Selection_Info(stream,selected_cell);}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_GRID_2D<T>::
Selection_Bounding_Box() const
{
    if(select_type==SELECT_TYPE::NODE) return World_Space_Box(RANGE<TV>(grid.Node(selected_node)));
    if(select_type==SELECT_TYPE::CELL) return World_Space_Box(RANGE<TV>(grid.Node(selected_cell),grid.Node(selected_cell+1)));
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_GRID_2D<float>;
template class OPENGL_GRID_2D<double>;
}
