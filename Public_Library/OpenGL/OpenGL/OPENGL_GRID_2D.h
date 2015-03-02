//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_2D
//##################################################################### 
#ifndef __OPENGL_GRID_2D__
#define __OPENGL_GRID_2D__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>

namespace PhysBAM
{

template<class T>
class OPENGL_GRID_2D:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    GRID<TV>      &grid;
    ARRAY<bool,TV_INT> *active_cell_mask,*ghost_cell_mask;
    ARRAY<bool,FACE_INDEX<TV::m> > *active_face_mask,*ghost_face_mask;
    ARRAY<bool,TV_INT> *active_node_mask,*ghost_node_mask;
    OPENGL_COLOR    color;
    bool draw;
    bool draw_ghost_values;
    int draw_mask_type;
private:
    OPENGL_SELECTION<T>* current_selection;
    std::string basedir;
    int frame;

public:
    OPENGL_GRID_2D(STREAM_TYPE stream_type,GRID<TV> &grid_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::White(),const std::string basedir_input="",const int frame_input=0)
        :OPENGL_OBJECT<T>(stream_type),grid(grid_input),active_cell_mask(0),ghost_cell_mask(0),active_face_mask(0),ghost_face_mask(0),active_node_mask(0),ghost_node_mask(0),color(color_input),draw(true),draw_ghost_values(true),draw_mask_type(0),current_selection(0),basedir(basedir_input),frame(frame_input)
    {Reinitialize();}

    void Display() const PHYSBAM_OVERRIDE;
    virtual void Set_Frame(int frame_input);
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;

    void Toggle_Draw_Ghost_Values();
    void Reinitialize();
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;

    DEFINE_CALLBACK_CREATOR(OPENGL_GRID_2D, Toggle_Draw_Ghost_Values);
};

template<class T>
class OPENGL_SELECTION_GRID_CELL_2D:public OPENGL_SELECTION<T>
{
private:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    using OPENGL_SELECTION<T>::object;
    VECTOR<int,2> index;
    OPENGL_SELECTION_GRID_CELL_2D(OPENGL_OBJECT<T>* object, const VECTOR<int,2> &index=TV_INT()) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::GRID_CELL_2D, object), index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_GRID_NODE_2D:public OPENGL_SELECTION<T>
{
private:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    using OPENGL_SELECTION<T>::object;
    VECTOR<int,2> index;
    OPENGL_SELECTION_GRID_NODE_2D(OPENGL_OBJECT<T>* object, const VECTOR<int,2> &index=TV_INT()) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::GRID_NODE_2D, object), index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
