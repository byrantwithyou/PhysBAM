//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SCALAR_FIELD_3D
//#####################################################################
#ifndef __OPENGL_SCALAR_FIELD_3D__
#define __OPENGL_SCALAR_FIELD_3D__

#include <Core/Arrays/ARRAY.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_GRID_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{
template<class T> class OPENGL_TEXTURED_RECT;
template<class T,class T_ARRAY> class OPENGL_POINTS_3D;

template<class T,class T2=T>
class OPENGL_SCALAR_FIELD_3D:public OPENGL_OBJECT<T>,public OPENGL_GRID_OBJECT<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    using OPENGL_OBJECT<T>::slice;using OPENGL_OBJECT<T>::stream_type;
    GRID<TV> grid;
    ARRAY<T2,VECTOR<int,3> > &values;

    enum DRAW_MODE { DRAW_TEXTURE, DRAW_POINTS };
    DRAW_MODE draw_mode;
    ARRAY<OPENGL_COLOR_MAP<T2>*> color_maps; // all owned by us
    int current_color_map;

    OPENGL_SCALAR_FIELD_3D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,
        ARRAY<T2,VECTOR<int,3> > &values_input,OPENGL_COLOR_MAP<T2> *color_map_input,
        DRAW_MODE draw_mode_input=DRAW_TEXTURE);
    virtual ~OPENGL_SCALAR_FIELD_3D();
    void Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map);

    void Set_Scale_Range(const T2 range_min,const T2 range_max);
    void Reset_Scale_Range();
    T2 Pre_Map_Value(const T2 value) const;

    void Display() const override;
    virtual RANGE<TV> Bounding_Box() const override;
    void Slice_Has_Changed() override;    

    void Set_Draw_Mode(DRAW_MODE draw_mode);
    virtual void Update();  // Call when values or other attributes have changed
    void Print_Cell_Selection_Info(std::ostream& stream,const TV_INT& cell) const override;
    void Print_Node_Selection_Info(std::ostream& stream,const TV_INT& node) const override;

    // convenience functions
    void Toggle_Draw_Mode();
    void Set_Smooth_Slice_Texture(bool smooth_slice_texture_input=true);
    void Toggle_Smooth_Slice_Texture();
    void Toggle_Color_Map();
    OPENGL_COLOR Do_Color(const int i,const int j,const int k) const;
    OPENGL_COLOR Do_Color(const TV_INT& index) const;

private:
    void Display_3D() const;
    void Display_3D_Slice() const;
    void Update_Slice();
    void Update_Points();
    void Delete_Textured_Rect();
    void Delete_Points();

public:
    OPENGL_TEXTURED_RECT<T>* opengl_textured_rect;
    OPENGL_POINTS_3D<T,ARRAY<VECTOR<T,3> > > *opengl_points;
    bool smooth_slice_texture;
    bool scale_range;
    T2 scale_range_min,scale_range_dx;
    int selected_point;
};

}

#endif
