//#####################################################################
// Copyright 2002-2004, Eran Guendelman, Eilene Hao, Neil Molino, Robert Bridson, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_B_SPLINE_PATCH
//##################################################################### 
#ifndef __OPENGL_B_SPLINE_PATCH__
#define __OPENGL_B_SPLINE_PATCH__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
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
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Point;
    using OPENGL_OBJECT<T>::World_Space_Box;
    OPENGL_B_SPLINE_PATCH(STREAM_TYPE stream_type,B_SPLINE_PATCH<VECTOR<T,3>,3>& patch_input,
                                const OPENGL_MATERIAL& front_material_input,const OPENGL_MATERIAL& back_material_input);
    virtual ~OPENGL_B_SPLINE_PATCH();

    void Display() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection,MATRIX<T,4>* transform) const;

    OPENGL_SELECTION<T>* Get_Vertex_Selection(int index);
    OPENGL_SELECTION<T>* Get_Element_Selection(int index);

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
    OPENGL_SELECTION<T>* current_selection;
    int current_node;
    bool highlight_current_node;
    bool wireframe_only;
    bool draw_particles;
    const int substeps=5;

    friend class OPENGL_SELECTION_B_SPLINE_PATCH_VERTEX<T>;
    friend class OPENGL_SELECTION_B_SPLINE_PATCH_ELEMENT<T>;
};

template<class T>
class OPENGL_SELECTION_B_SPLINE_PATCH_VERTEX:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_B_SPLINE_PATCH_VERTEX(OPENGL_OBJECT<T>* object, int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::B_SPLINE_PATCH_VERTEX, object), index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_B_SPLINE_PATCH_ELEMENT:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_B_SPLINE_PATCH_ELEMENT(OPENGL_OBJECT<T>* object, int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::B_SPLINE_PATCH_ELEMENT, object), index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}
#endif
