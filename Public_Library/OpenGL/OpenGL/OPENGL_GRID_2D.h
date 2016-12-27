//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_2D
//##################################################################### 
#ifndef __OPENGL_GRID_2D__
#define __OPENGL_GRID_2D__

#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_GRID_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{

template<class T>
class OPENGL_GRID_2D:public OPENGL_OBJECT<T>
{
public:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    using OPENGL_OBJECT<T>::viewer_callbacks;
    GRID<TV>& grid;
    ARRAY<bool,TV_INT> *active_cell_mask,*ghost_cell_mask;
    ARRAY<bool,FACE_INDEX<TV::m> > *active_face_mask,*ghost_face_mask;
    ARRAY<bool,TV_INT> *active_node_mask,*ghost_node_mask;
    OPENGL_COLOR    color;
    bool draw;
    bool draw_ghost_values;
    int draw_mask_type;
    ARRAY<OPENGL_GRID_OBJECT<TV>*> grid_objects;
private:
    TV_INT selected_cell;
    TV_INT selected_node;
    std::string basedir;
    int frame;

public:
    OPENGL_GRID_2D(STREAM_TYPE stream_type,GRID<TV> &grid_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::White(),
        const std::string basedir_input="",const int frame_input=0);

    void Display() const override;
    virtual void Set_Frame(int frame_input);
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;

    void Toggle_Draw_Ghost_Values();
    void Reinitialize();
    void Print_Selection_Info(std::ostream& stream) const override;
};

}

#endif
