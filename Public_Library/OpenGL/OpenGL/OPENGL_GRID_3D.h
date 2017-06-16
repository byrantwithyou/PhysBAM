//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Michael Lentine, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_3D
//##################################################################### 
#ifndef __OPENGL_GRID_3D__
#define __OPENGL_GRID_3D__
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_GRID_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T>
class OPENGL_GRID_3D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;

    enum class SELECT_TYPE
    {
        NONE,
        CELL,
        NODE
    };
    SELECT_TYPE select_type;
    TV_INT selected_cell;
    TV_INT selected_node;
    ARRAY<TV_INT> selected_cell_list;
    ARRAY<TV_INT> selected_node_list;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    using OPENGL_OBJECT<T>::viewer_callbacks;

    GRID<TV> &grid;
    OPENGL_COLOR color;
    bool draw_ghost_values;
    bool hide_non_selected_grid;
    bool owns_grid;
    int scale;
    ARRAY<OPENGL_GRID_OBJECT<TV>*> grid_objects;

//##################################################################### 
    OPENGL_GRID_3D(STREAM_TYPE stream_type,GRID<TV> &grid_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::White());
    virtual ~OPENGL_GRID_3D();
    void Display() const override;
    virtual RANGE<TV> Bounding_Box() const override;
    void Print_Selection_Info(std::ostream& stream) const override;
    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    virtual RANGE<TV> Selection_Bounding_Box() const override;
    void Clear_Selection() override;
    void Toggle_Draw_Ghost_Values();
private:
    void Draw_Subgrid(const TV_INT &node_start,const TV_INT &node_end) const;
    void Draw_Nodes_For_Selection(const TV_INT &node_start,const TV_INT &node_end) const;
//##################################################################### 
};

}

#endif
