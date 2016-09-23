//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SYMMETRIC_MATRIX_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_SYMMETRIC_MATRIX_FIELD_3D__
#define __OPENGL_SYMMETRIC_MATRIX_FIELD_3D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T>
class OPENGL_SYMMETRIC_MATRIX_FIELD_3D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    GRID<TV> grid;
    const ARRAY<SYMMETRIC_MATRIX<T,3>,VECTOR<int,3> >& field;
    ARRAY<TRIPLE<VECTOR<T,3>,MATRIX<T,3>,VECTOR<bool,3> > > entries;
    OPENGL_COLOR positive_color,negative_color;
    T size;

    OPENGL_SYMMETRIC_MATRIX_FIELD_3D(STREAM_TYPE stream_type,const GRID<TV>& grid_input,const ARRAY<SYMMETRIC_MATRIX<T,3>,VECTOR<int,3> >& field_input,const T size_input=0.025,
        const OPENGL_COLOR& positive_color_input=OPENGL_COLOR::Red(),const OPENGL_COLOR& negative_color_input=OPENGL_COLOR::Yellow())
        :OPENGL_OBJECT<T>(stream_type),grid(grid_input),field(field_input),positive_color(positive_color_input),negative_color(negative_color_input),size(size_input)
    {}

    void Slice_Has_Changed() override
    {Update();}

    void Display() const override;
    virtual void Update();
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
//##################################################################### 
};
}
#endif
