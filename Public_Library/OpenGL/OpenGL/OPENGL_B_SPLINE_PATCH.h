//#####################################################################
// Copyright 2002-2004, Eran Guendelman, Eilene Hao, Neil Molino, Robert Bridson, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_B_SPLINE_PATCH
//##################################################################### 
#ifndef __OPENGL_B_SPLINE_PATCH__
#define __OPENGL_B_SPLINE_PATCH__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class TV,int d> class B_SPLINE_PATCH;
template<class T> class TRIANGULATED_PATCH;
template<class T> class OPENGL_SELECTION_B_SPLINE_PATCH_VERTEX;
template<class T> class OPENGL_SELECTION_B_SPLINE_PATCH_ELEMENT;

template<class T>
class OPENGL_B_SPLINE_PATCH:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Point;
    using OPENGL_OBJECT<T>::World_Space_Box;
    OPENGL_B_SPLINE_PATCH(STREAM_TYPE stream_type,B_SPLINE_PATCH<VECTOR<T,3>,3>& patch_input,
                                const OPENGL_MATERIAL& front_material_input,const OPENGL_MATERIAL& back_material_input);
    virtual ~OPENGL_B_SPLINE_PATCH();

    void Display() const override;
    virtual RANGE<TV> Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream &output_stream) const override;
    void Print_Selection_Info(std::ostream &output_stream,MATRIX<T,4>* transform) const;
    virtual RANGE<TV> Selection_Bounding_Box() const override;

    void Get_Vertex_Selection(int index);
    void Get_Element_Selection(int index);

    void Set_Two_Sided(bool two_sided_input=true) {two_sided=two_sided_input;}
    void Set_Front_Material(const OPENGL_MATERIAL &material_input);
    void Set_Back_Material(const OPENGL_MATERIAL &material_input);

    // Create_Display_Lists creates a display list which this object will own.
    // Use_Display_List makes this object use a display list owned by another object.
    int Create_Display_List();
    int Get_Display_List_Id() { return display_list_id; }
    void Use_Display_List(int input_display_list_id);

    void Highlight_Current_Node() const;

public:
    void Reinitialize_Display_List();
    void Draw() const;
    void Draw_Selected_Element(int index) const;
    void Draw_Vertices_For_Selection() const;
    void Draw_Elements_For_Selection() const;

    B_SPLINE_PATCH<VECTOR<T,3>,3>& patch;
    bool two_sided;
public:
    OPENGL_MATERIAL front_material;
    OPENGL_MATERIAL back_material;

protected:
    bool use_display_list;
    bool owns_display_list;
    int display_list_id;

public:
    int selected_vertex;
    int selected_element;
    int current_node;
    bool highlight_current_node;
    bool wireframe_only;
    bool draw_particles;
    const int substeps=5;

    friend class OPENGL_SELECTION_B_SPLINE_PATCH_VERTEX<T>;
    friend class OPENGL_SELECTION_B_SPLINE_PATCH_ELEMENT<T>;
};

}
#endif
