//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Michael Lentine, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_GRID_3D<T>::
OPENGL_GRID_3D(STREAM_TYPE stream_type,GRID<TV> &grid_input,const OPENGL_COLOR &color_input) 
    :OPENGL_OBJECT<T>(stream_type),select_type(SELECT_TYPE::NONE),grid(grid_input),color(color_input),
    draw_ghost_values(true),hide_non_selected_grid(false),owns_grid(false),scale(1)
{
    viewer_callbacks.Set("toggle_draw_ghost_values",{[this](){Toggle_Draw_Ghost_Values();},"toggle_draw_ghost_values"});
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_GRID_3D<T>::
~OPENGL_GRID_3D()
{
    if(owns_grid) delete &grid;
}
//#####################################################################
// Display
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Display() const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);

    color.Send_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)this->slice;

    int ghost_cells=draw_ghost_values?3:0;

    if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE)
    {
        VECTOR<int,3> node_start(0,0,0),node_end(grid.numbers_of_cells+1);
        if(!hide_non_selected_grid) Draw_Subgrid(node_start,node_end);
    }
    else if(slice->mode==OPENGL_SLICE::NODE_SLICE)
    {
        VECTOR<int,3> node_start(0,0,0),node_end(grid.numbers_of_cells+1);
        if(slice->axis==0) { node_start.x=node_end.x=slice->index/scale; }
        else if(slice->axis==1) { node_start.y=node_end.y=slice->index/scale+1; }
        else if(slice->axis==2) { node_start.z=node_end.z=slice->index/scale+1; }

        if(mode==GL_SELECT)
        {
            // Currently only support node selection in this mode
            glPushName(1);
            Draw_Nodes_For_Selection(node_start,node_end);
            glPopName();
        }
        else Draw_Subgrid(node_start,node_end);
    }
    else if(slice->mode==OPENGL_SLICE::CELL_SLICE)
    {
        if(mode==GL_SELECT)
        {
            // Currently only support cell selection in this mode
            glPushAttrib(GL_ENABLE_BIT);
            glDisable(GL_CULL_FACE);

            TV x_vector(grid.dX.x,0,0),y_vector(0,grid.dX.y,0),z_vector(0,0,grid.dX.z);

            T x, y, z;
            int i, j, k;
            int i_start, i_end, j_start, j_end, k_start, k_end;
            TV axis_1, axis_2, axis_3;
            if(slice->axis==0) { i_start=slice->index/scale;i_end=i_start+1; axis_1=y_vector; axis_2=z_vector; axis_3=x_vector; } 
            else { i_start=1-ghost_cells; i_end=grid.numbers_of_cells.x+ghost_cells; }
            if(slice->axis==1) { j_start=slice->index/scale;j_end=j_start+1; axis_1=z_vector; axis_2=x_vector; axis_3=y_vector; } 
            else { j_start=1-ghost_cells; j_end=grid.numbers_of_cells.y+ghost_cells; }
            if(slice->axis==2) { k_start=slice->index/scale;k_end=k_start+1; axis_1=x_vector; axis_2=y_vector; axis_3=z_vector; } 
            else { k_start=1-ghost_cells; k_end=grid.numbers_of_cells.z+ghost_cells; }

            TV pos_start=grid.Node(TV_INT(i_start,j_start,k_start));

            glPushName(0);
            for(i=i_start, x=pos_start.x; i<i_end; i++, x+=grid.dX.x)
            {
                glPushName(i);
                for(j=j_start, y=pos_start.y; j<j_end; j++, y+=grid.dX.y)
                {
                    glPushName(j);
                    for(k=k_start, z=pos_start.z; k<k_end; k++, z+=grid.dX.z)
                    {
                        TV min_corner(x,y,z);
                        glPushName(k);
                        OpenGL_Begin(GL_QUADS);
                        OpenGL_Quad(min_corner,axis_1,axis_2);
                        OpenGL_Quad(min_corner+axis_3,axis_1,axis_2);
                        OpenGL_End();
                        glPopName();
                    }
                    glPopName();
                }
                glPopName();
            }

            glPopName();
            glPopAttrib();
        }
        else
        {
            VECTOR<int,3> node_start(-ghost_cells,-ghost_cells,-ghost_cells),node_end(grid.numbers_of_cells+ghost_cells+1);
            if(slice->axis==0) { node_start.x=slice->index/scale; node_end.x=node_start.x+1; }
            else if(slice->axis==1) { node_start.y=slice->index/scale; node_end.y=node_start.y+1; }
            else if(slice->axis==2) { node_start.z=slice->index/scale; node_end.z=node_start.z+1; }
            Draw_Subgrid(node_start,node_end);

            // Outline boundary of real domain in wider line
            if(ghost_cells>0){
                TV x000,x111;
                if(slice->axis==0){x000=TV(grid.domain.min_corner.x+grid.dX.x*slice->index/scale,grid.domain.min_corner.y,grid.domain.min_corner.z);x111=TV(grid.domain.min_corner.x+grid.dX.x*slice->index/scale,grid.domain.max_corner.y,grid.domain.max_corner.z);}
                else if(slice->axis==1){x000=TV(grid.domain.min_corner.x,grid.domain.min_corner.y+grid.dX.y*slice->index/scale,grid.domain.min_corner.z);x111=TV(grid.domain.max_corner.x,grid.domain.min_corner.y+grid.dX.y*slice->index/scale,grid.domain.max_corner.z);}
                else if(slice->axis==2){x000=TV(grid.domain.min_corner.x,grid.domain.min_corner.y,grid.domain.min_corner.z+grid.dX.z*slice->index/scale);x111=TV(grid.domain.max_corner.x,grid.domain.max_corner.y,grid.domain.min_corner.z+grid.dX.z*slice->index/scale);}
                glPushAttrib(GL_LINE_BIT);glLineWidth(3*OPENGL_PREFERENCES::line_width);
                TV x001(x000.x,x000.y,x111.z),x010(x000.x,x111.y,x000.z),x011(x000.x,x111.y,x111.z),
                    x100(x111.x,x000.y,x000.z),x101(x111.x,x000.y,x111.z),x110(x111.x,x111.y,x000.z);
                
                OpenGL_Begin(GL_LINES);
                OpenGL_Line(x000,x010);
                OpenGL_Line(x010,x011);
                OpenGL_Line(x011,x001);
                OpenGL_Line(x001,x000);

                OpenGL_Line(x100,x110);
                OpenGL_Line(x110,x111);
                OpenGL_Line(x111,x101);
                OpenGL_Line(x101,x100);
                
                OpenGL_Line(x000,x100);
                OpenGL_Line(x010,x110);
                OpenGL_Line(x011,x111);
                OpenGL_Line(x001,x101);
                OpenGL_End();
                glPopAttrib();}
        }
    }

    if(select_type==SELECT_TYPE::CELL)
    {
        TV min_corner=grid.Node(selected_cell),max_corner=min_corner+grid.dX;
        OPENGL_SELECTION::Draw_Highlighted_Box(min_corner,max_corner);
    }
    else if(select_type==SELECT_TYPE::NODE)
    {
        OPENGL_SELECTION::Draw_Highlighted_Vertex(grid.Node(selected_node));
    }
    else if(selected_cell_list.m)
    {
        for(int i=0;i<selected_cell_list.m;i++){
            TV min_corner=grid.Node(selected_cell_list(i)),max_corner=min_corner+grid.dX;
            OPENGL_SELECTION::Draw_Highlighted_Box(min_corner,max_corner);}
    }
    else if(selected_node_list.m)
    {
        for(int i=0;i<selected_node_list.m;i++) OPENGL_SELECTION::Draw_Highlighted_Vertex(grid.Node(selected_node_list(i)));
    }

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Draw_Subgrid
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Draw_Subgrid(const VECTOR<int,3> &node_start,const VECTOR<int,3> &node_end) const
{
    int i,j,k;
    T x,y,z;

    TV start_position=grid.Node(node_start),end_position=grid.Node(node_end-1);

    OpenGL_Begin(GL_LINES);

    if(node_start.z!=node_end.z)
        for(i=node_start.x, x=start_position.x; i<node_end.x; i++, x+=grid.dX.x)
            for(j=node_start.y, y=start_position.y; j<node_end.y; j++, y+=grid.dX.y)
                OpenGL_Line(TV(x,y,start_position.z),TV(x,y,end_position.z));

    if(node_start.y!=node_end.y)
        for(i=node_start.x, x=start_position.x; i<node_end.x; i++, x+=grid.dX.x)
            for(k=node_start.z, z=start_position.z; k<node_end.z; k++, z+=grid.dX.z)
                OpenGL_Line(TV(x,start_position.y,z),TV(x,end_position.y,z));

    if(node_start.x!=node_end.x)
        for(j=node_start.y, y=start_position.y; j<node_end.y; j++, y+=grid.dX.y)
            for(k=node_start.z, z=start_position.z; k<node_end.z; k++, z+=grid.dX.z)
                OpenGL_Line(TV(start_position.x,y,z),TV(end_position.x,y,z));

    OpenGL_End();
}
//#####################################################################
// Draw_Nodes_For_Selection
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Draw_Nodes_For_Selection(const VECTOR<int,3> &node_start,const VECTOR<int,3> &node_end) const
{
    int i,j,k;
    T x,y,z;

    TV start_position=grid.Node(node_start);

    glPushAttrib(GL_POINT_BIT);
    glPointSize(OPENGL_PREFERENCES::selection_point_size);

    for(i=node_start.x, x=start_position.x; i<node_end.x; i++, x+=grid.dX.x)
    {
        glPushName(i);
        for(j=node_start.y, y=start_position.y; j<node_end.y; j++, y+=grid.dX.y)
        {
            glPushName(j);
            for(k=node_start.z, z=start_position.z; k<node_end.z; k++, z+=grid.dX.z)
            {
                glPushName(k);
                OpenGL_Begin(GL_POINTS);
                OpenGL_Vertex(TV(x,y,z));
                OpenGL_End();
                glPopName();
            }
            glPopName();
        }
        glPopName();
    }

    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_GRID_3D<T>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_GRID_3D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    const static int priority[]={60,70,60,70};
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_GRID_3D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    ARRAY<int> signed_indices(indices.m);
    for(int i=1;i<indices.m;i++)
        if(indices(i)&(GLuint)(-1))
            signed_indices(i)=-~(indices(i)-1);
        else
            signed_indices(i)=indices(i);
    if(indices.m==4)
    {
        if(indices(0)==0){
            select_type=SELECT_TYPE::CELL;
            selected_cell=TV_INT(signed_indices(1),signed_indices(2),signed_indices(3));}
        else if(indices(0)==1){
            select_type=SELECT_TYPE::NODE;
            selected_node=TV_INT(signed_indices(1),signed_indices(2),signed_indices(3));}
        else{
            select_type=SELECT_TYPE::NONE;
            return false;}
        return true;
    } 
    else if(indices.m%3==1)
    {
        if(indices(0)==2)
        {
            selected_cell_list.Resize(indices.m/3);
            for(int i=0;i<selected_cell_list.m;i++)
                selected_cell_list(i)=TV_INT(signed_indices(3*i+1),signed_indices(3*i+2),signed_indices(3*i+3));
        }
        else if(indices(0)==3)
        {
            selected_node_list.Resize(indices.m/3);
            for(int i=0;i<selected_node_list.m;i++)
                selected_node_list(i)=TV_INT(signed_indices(3*i+1),signed_indices(3*i+2),signed_indices(3*i+3));
        }
        else return false;
        return true;
    }
    return false;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Clear_Selection()
{
    select_type=SELECT_TYPE::NONE;
    selected_cell_list.Remove_All();
    selected_node_list.Remove_All();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_GRID_3D<T>::
Selection_Bounding_Box() const
{
    if(select_type==SELECT_TYPE::CELL)
        return World_Space_Box(grid.Cell_Domain(selected_cell));
    if(select_type==SELECT_TYPE::NODE)
        return World_Space_Box(RANGE<TV>(grid.Node(selected_node)));
    RANGE<TV> range;
    for(int i=0;i<selected_cell_list.m;i++)
        range=range.Unite(grid.Cell_Domain(selected_cell_list(i)));
    for(int i=0;i<selected_node_list.m;i++)
        range.Enlarge_To_Include_Point(grid.Node(selected_node_list(i)));
    return World_Space_Box(range);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
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
template<class T> void OPENGL_GRID_3D<T>::
Toggle_Draw_Ghost_Values()
{
    draw_ghost_values=!draw_ghost_values;
}
namespace PhysBAM{
template class OPENGL_GRID_3D<float>;
template class OPENGL_GRID_3D<double>;
}
