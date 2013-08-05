//#####################################################################
// Copyright 2004, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_LEVELSET_2D
//#####################################################################
#ifndef __OPENGL_LEVELSET_2D__
#define __OPENGL_LEVELSET_2D__

#include <Geometry/Level_Sets/LEVELSET.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
namespace PhysBAM{

template<class T>
class OPENGL_LEVELSET_2D:public OPENGL_SCALAR_FIELD_2D<T,T>
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    using OPENGL_SCALAR_FIELD_2D<T,T>::grid;
    using OPENGL_SCALAR_FIELD_2D<T,T>::Send_Transform_To_GL_Pipeline;
    using OPENGL_SCALAR_FIELD_2D<T,T>::active_cells;

    LEVELSET<TV>& levelset;
    OPENGL_TRIANGULATED_AREA<T>* opengl_triangulated_area;
    OPENGL_SEGMENTED_CURVE_2D<T>* opengl_segmented_curve_2d;
    enum COLOR_MODE {COLOR_SOLID,COLOR_GRADIENT};
    COLOR_MODE color_mode;
    OPENGL_COLOR inside_color, outside_color;
    OPENGL_COLOR_MAP<T> *solid_color_map,*gradient_color_map;
    bool draw_cells,draw_area,draw_curve;
    int dominant_sign; // interesting sign -1 for water 1 for fire
private:
    bool draw_normals;

public:
    OPENGL_LEVELSET_2D(LEVELSET<TV>& levelset_input,const OPENGL_COLOR& inside_color_input=OPENGL_COLOR::Blue(),
        const OPENGL_COLOR& outside_color_input=OPENGL_COLOR::Red((T).5),ARRAY<bool,VECTOR<int,2> > *active_cells_input=0);

    virtual ~OPENGL_LEVELSET_2D();

//#####################################################################
    void Set_Inside_And_Outside_Colors(const OPENGL_COLOR& inside_color_input,const OPENGL_COLOR& outside_color_input);
    void Toggle_Normals();
    void Update() PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
private:
//#####################################################################
};
}
#endif
