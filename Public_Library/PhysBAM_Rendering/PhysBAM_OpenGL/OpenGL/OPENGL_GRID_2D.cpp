//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Geoffrey Irving, Sndrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Display(const int in_color) const
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
#ifndef USE_OPENGLES
    glGetIntegerv(GL_RENDER_MODE, &mode);
#endif

    int ghost_cells=draw_ghost_values?3:0;

#ifndef USE_OPENGLES
    if(mode == GL_SELECT){
        glPushAttrib(GL_ENABLE_BIT|GL_POINT_BIT);glDisable(GL_CULL_FACE);
        VECTOR<T,2> min_corner=grid.Node(-ghost_cells,-ghost_cells),X;int i,j;

        // Draw grid cells for selection
        glPushName(1);
        for(i=-ghost_cells,X.x=min_corner.x;i<grid.numbers_of_cells.x+ghost_cells;i++,X.x+=grid.dX.x){
            glPushName(i);
            for(j=-ghost_cells,X.y=min_corner.y;j<grid.numbers_of_cells.y+ghost_cells;j++,X.y+=grid.dX.y){
                glPushName(j);
                ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
                OpenGL_Quad_2D(X,X+grid.dX,vertices);
                OpenGL_Draw_Arrays(GL_QUADS,2,vertices);
                glPopName();}
            glPopName();}

        // Draw grid nodes for selection
        glLoadName(2);
        glPointSize(OPENGL_PREFERENCES::selection_point_size);
        for(i=-ghost_cells,X.x=min_corner.x;i<grid.numbers_of_cells.x+ghost_cells+1;i++,X.x+=grid.dX.x){
            glPushName(i);
            for(j=-ghost_cells,X.y=min_corner.y;j<grid.numbers_of_cells.y+ghost_cells+1;j++,X.y+=grid.dX.y){
                glPushName(j);
                ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
                OpenGL_Vertex(X,vertices);
                OpenGL_Draw_Arrays(GL_POINTS,2,vertices);
                glPopName();}
            glPopName();}
        glPopName();glPopAttrib();}
    else
#endif
    {
        // Draw masks
        T x,y;int i,j,i_mask=1,j_mask=1;
        VECTOR<T,2> min_corner=grid.Node(-ghost_cells,-ghost_cells);
        VECTOR<T,2> max_corner=grid.Node(grid.numbers_of_cells.x+ghost_cells,grid.numbers_of_cells.y+ghost_cells);

        T_ARRAYS_BOOL *cell_mask=0;
        T_FACE_ARRAYS_BOOL *face_mask=0;
        T_ARRAYS_BOOL *node_mask=0;
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
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells;j++,y+=grid.dX.y,j_mask++)
                    if((*cell_mask)(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(VECTOR<T,2>(x,y),VECTOR<T,2>(x+grid.dX.x,y+grid.dX.y),vertices);}
            OpenGL_Draw_Arrays(GL_LINES,2,vertices);
            glPopAttrib();}

        if(face_mask){
            int mask_m_start,mask_n_start;
            glPushAttrib(GL_ALL_ATTRIB_BITS);glLineWidth(4*OPENGL_PREFERENCES::line_width);OPENGL_COLOR::Cyan().Send_To_GL_Pipeline();
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            mask_m_start=face_mask->Component(0).domain.min_corner.x,mask_n_start=face_mask->Component(0).domain.min_corner.y;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells+1;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells;j++,y+=grid.dX.y,j_mask++)
                    if((face_mask->Component(0))(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(VECTOR<T,2>(x,y),VECTOR<T,2>(x,y+grid.dX.y),vertices);}
            mask_m_start=face_mask->Component(1).domain.min_corner.x,mask_n_start=face_mask->Component(1).domain.min_corner.y;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells+1;j++,y+=grid.dX.y,j_mask++)
                    if((face_mask->Component(1))(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(VECTOR<T,2>(x,y),VECTOR<T,2>(x+grid.dX.x,y),vertices);}
            OpenGL_Draw_Arrays(GL_LINES,2,vertices);
            glPopAttrib();}

        if(node_mask){
            int mask_m_start=node_mask->domain.min_corner.x,mask_n_start=node_mask->domain.min_corner.y;
            glPushAttrib(GL_ALL_ATTRIB_BITS);glPointSize(5);OPENGL_COLOR::Cyan().Send_To_GL_Pipeline();
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells+1;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells+1;j++,y+=grid.dX.y,j_mask++)
                    if((*node_mask)(mask_m_start+i_mask-1,mask_n_start+j_mask-1)) OpenGL_Vertex(VECTOR<T,2>(x,y),vertices);
            OpenGL_Draw_Arrays(GL_POINTS,2,vertices);
            glPopAttrib();}

        // Draw grid
        if(active_face_mask){
            if(draw_mask_type&&ghost_cells==4){ghost_cells=3;min_corner=min_corner+VECTOR<T,2>(grid.dX.x,grid.dX.y);max_corner=max_corner-VECTOR<T,2>(grid.dX.x,grid.dX.y);}
            face_mask=active_face_mask;
            int mask_m_start,mask_n_start;
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            mask_m_start=face_mask->Component(0).domain.min_corner.x,mask_n_start=face_mask->Component(0).domain.min_corner.y;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells+1;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells;j++,y+=grid.dX.y,j_mask++)
                    if((face_mask->Component(0))(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(VECTOR<T,2>(x,y),VECTOR<T,2>(x,y+grid.dX.y),vertices);}
            mask_m_start=face_mask->Component(1).domain.min_corner.x,mask_n_start=face_mask->Component(1).domain.min_corner.y;
            for(i=-ghost_cells,x=min_corner.x,i_mask=1;i<grid.numbers_of_cells.x+ghost_cells;i++,x+=grid.dX.x,i_mask++)
                for(j=-ghost_cells,y=min_corner.y,j_mask=1;j<grid.numbers_of_cells.y+ghost_cells+1;j++,y+=grid.dX.y,j_mask++)
                    if((face_mask->Component(1))(mask_m_start+i_mask-1,mask_n_start+j_mask-1)){OpenGL_Line(VECTOR<T,2>(x,y),VECTOR<T,2>(x+grid.dX.x,y),vertices);}
            OpenGL_Draw_Arrays(GL_LINES,2,vertices);}
        else{
            if(draw_mask_type&&ghost_cells==3){ghost_cells=4;min_corner=min_corner-VECTOR<T,2>(grid.dX.x,grid.dX.y);max_corner=max_corner+VECTOR<T,2>(grid.dX.x,grid.dX.y);}
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            for(i=-ghost_cells,x=min_corner.x;i<grid.numbers_of_cells.x+ghost_cells+1;i++,x+=grid.dX.x){
                OpenGL_Line(VECTOR<T,2>(x,min_corner.y),VECTOR<T,2>(x,max_corner.y),vertices);}
            for(j=-ghost_cells,y=min_corner.y;j<grid.numbers_of_cells.y+ghost_cells+1;j++,y+=grid.dX.y){
                OpenGL_Line(VECTOR<T,2>(min_corner.x,y),VECTOR<T,2>(max_corner.x,y),vertices);}
            OpenGL_Draw_Arrays(GL_LINES,2,vertices);}

        // Outline boundary of real domain in wider line
        if(ghost_cells>0){
            glPushAttrib(GL_LINE_BIT);glLineWidth(2*OPENGL_PREFERENCES::line_width);
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            OpenGL_Quad_2D(grid.domain.min_corner,grid.domain.max_corner,vertices);
            OpenGL_Draw_Arrays(GL_LINE_LOOP,2,vertices);
            glPopAttrib();}

        // Highlight current selection
        if(current_selection){
            if(current_selection->type == OPENGL_SELECTION::GRID_CELL_2D){
                OPENGL_SELECTION_GRID_CELL_2D<T>* real_selection=(OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection;
                VECTOR<T,2> min_corner=grid.Node(real_selection->index),max_corner=min_corner+grid.dX;
                OPENGL_SELECTION::Draw_Highlighted_Quad(min_corner,max_corner);}
            else if(current_selection->type == OPENGL_SELECTION::GRID_NODE_2D){
                OPENGL_SELECTION_GRID_NODE_2D<T>* real_selection=(OPENGL_SELECTION_GRID_NODE_2D<T>*)current_selection;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(grid.Node(real_selection->index));}}}

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
template<class T> RANGE<VECTOR<float,3> > OPENGL_GRID_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,2> >(grid.domain));
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_GRID_2D<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size == 3){
        if(buffer[0] == 1) selection=new OPENGL_SELECTION_GRID_CELL_2D<T>(this,VECTOR<int,2>(buffer[1],buffer[2]));
        else if(buffer[0] == 2) selection=new OPENGL_SELECTION_GRID_NODE_2D<T>(this,VECTOR<int,2>(buffer[1],buffer[2]));}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if (selection->type == OPENGL_SELECTION::GRID_CELL_2D){
        OPENGL_SELECTION_GRID_CELL_2D<T>* real_selection=(OPENGL_SELECTION_GRID_CELL_2D<T>*)selection;
        current_selection=new OPENGL_SELECTION_GRID_CELL_2D<T>(this,real_selection->index);}
    else if(selection->type == OPENGL_SELECTION::GRID_NODE_2D){
        OPENGL_SELECTION_GRID_NODE_2D<T>* real_selection=(OPENGL_SELECTION_GRID_NODE_2D<T>*)selection;
        current_selection=new OPENGL_SELECTION_GRID_NODE_2D<T>(this,real_selection->index);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
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
    std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/active_cell_mask",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        if(!active_cell_mask) active_cell_mask=new T_ARRAYS_BOOL();
        active_cell_mask->Clean_Memory();
        FILE_UTILITIES::Read_From_File<bool>(filename,*active_cell_mask);}

    filename=STRING_UTILITIES::string_sprintf("%s/%d/ghost_cell_mask",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        if(!ghost_cell_mask) ghost_cell_mask=new T_ARRAYS_BOOL();
        ghost_cell_mask->Clean_Memory();
        FILE_UTILITIES::Read_From_File<bool>(filename,*ghost_cell_mask);}

    filename=STRING_UTILITIES::string_sprintf("%s/%d/active_face_mask",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        if(!active_face_mask) active_face_mask=new T_FACE_ARRAYS_BOOL();
        active_face_mask->Clean_Memory();
        FILE_UTILITIES::Read_From_File<bool>(filename,*active_face_mask);}

    filename=STRING_UTILITIES::string_sprintf("%s/%d/ghost_face_mask",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        if(!ghost_face_mask) ghost_face_mask=new T_FACE_ARRAYS_BOOL();
        ghost_face_mask->Clean_Memory();
        FILE_UTILITIES::Read_From_File<bool>(filename,*ghost_face_mask);}

    filename=STRING_UTILITIES::string_sprintf("%s/%d/active_node_mask",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        if(!active_node_mask) active_node_mask=new T_ARRAYS_BOOL();
        active_node_mask->Clean_Memory();
        FILE_UTILITIES::Read_From_File<bool>(filename,*active_node_mask);}

    filename=STRING_UTILITIES::string_sprintf("%s/%d/ghost_node_mask",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        if(!ghost_node_mask) ghost_node_mask=new T_ARRAYS_BOOL();
        ghost_node_mask->Clean_Memory();
        FILE_UTILITIES::Read_From_File<bool>(filename,*ghost_node_mask);}
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_GRID_2D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_2D){
        VECTOR<int,2> index=((OPENGL_SELECTION_GRID_NODE_2D<T>*)current_selection)->index;
        stream<<"Selected node "<<index<<" ("<<grid.Get_Regular_Grid().X(index)<<")"<<std::endl;}
    else if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D){
        VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)->index;
        stream<<"Selected cell "<<index<<" ("<<grid.Get_MAC_Grid().X(index)<<")"<<std::endl;}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_GRID_CELL_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV>& grid=((OPENGL_GRID_2D<T>*)object)->grid;
    RANGE<VECTOR<T,2> > box(grid.Node(index),grid.Node(index.x+1,index.y+1));
    return object->World_Space_Box(RANGE<VECTOR<float,2> >(box));
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_GRID_NODE_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV>& grid=((OPENGL_GRID_2D<T>*)object)->grid;
    RANGE<VECTOR<T,2> > box(grid.Node(index));
    return object->World_Space_Box(RANGE<VECTOR<float,2> >(box));
}
//#####################################################################
template class OPENGL_GRID_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_GRID_2D<double>;
#endif
