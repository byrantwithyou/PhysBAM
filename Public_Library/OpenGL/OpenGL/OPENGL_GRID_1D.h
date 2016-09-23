//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_1D
//##################################################################### 
#ifndef __OPENGL_GRID_1D__
#define __OPENGL_GRID_1D__

#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>

namespace PhysBAM
{

template<class T>
class OPENGL_GRID_1D:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;using OPENGL_OBJECT<T>::viewer_callbacks;
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
    GRID<TV> &grid;
    OPENGL_COLOR color;
    bool draw;
    bool draw_ghost_values;
private:
    std::string basedir;
    int frame;

public:
    OPENGL_GRID_1D(STREAM_TYPE stream_type,GRID<TV> &grid_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::White(),const std::string basedir_input="",const int frame_input=0)
        :OPENGL_OBJECT<T>(stream_type),grid(grid_input),color(color_input),draw(true),draw_ghost_values(true),basedir(basedir_input),frame(frame_input)
    {
        viewer_callbacks.Set("toggle_draw_ghost_values",{[this](){Toggle_Draw_Ghost_Values();},"toggle_draw_ghost_values"});
    }

    void Display() const override;
    virtual void Set_Frame(int frame_input);
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    void Toggle_Draw_Ghost_Values();
};

template<class T>
class OPENGL_SELECTION_GRID_CELL_1D:public OPENGL_SELECTION<T>
{
private:
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
public:
    using OPENGL_SELECTION<T>::object;
    TV_INT index;
    OPENGL_SELECTION_GRID_CELL_1D(OPENGL_OBJECT<T>* object,const TV_INT& index=TV_INT()) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::GRID_CELL_1D,object),index(index)
    {}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

}

#endif
