//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D__
#define __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D.h>
namespace PhysBAM{

template<class T> class OPENGL_POINT_SIMPLICES_1D;
template<class T> class OPENGL_AXES;
template<class TV> class RIGID_BODY_COLLECTION;

template<class T>
class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D:public OPENGL_COMPONENT<T>
{
protected:
    typedef VECTOR<T,1> TV;
    typedef VECTOR<T,TV::m+(TV::m==1)> TV_BOX;

    bool valid;
    bool show_object_names;
    bool output_positions;
    bool draw_velocity_vectors;
    bool draw_node_velocity_vectors;
    bool draw_point_simplices;
    ARRAY<int> needs_init,needs_destroy;
    bool has_init_destroy_information;
public:
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
protected:
    ARRAY<OPENGL_POINT_SIMPLICES_1D<T>*,int> opengl_point_simplices;
    ARRAY<OPENGL_AXES<T>*,int> opengl_axes;
    ARRAY<bool,int> draw_object;
    ARRAY<bool,int> use_object_bounding_box;
    ARRAY<VECTOR<T,1> > positions;
    ARRAY<VECTOR<T,1> > node_positions;
//    void current_selection;
    bool need_destroy_rigid_body_collection;

public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_OBJECT<T>::viewer_callbacks;using OPENGL_COMPONENT<T>::viewer_dir;
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(const VIEWER_DIR& viewer_dir);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(const VIEWER_DIR& viewer_dir,RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    virtual ~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D();
    
//#####################################################################
    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    void Set_Draw_Object(int i, bool draw_it);  // Need to call Reinitialize after changing draw objects
    bool Get_Draw_Object(int i) const;
    void Set_Object_Color(int i, const OPENGL_COLOR &color);
    void Set_Use_Object_Bounding_Box(int i, bool use_it);

    void Toggle_Output_Positions();
    void Toggle_Show_Object_Names();
    void Toggle_Draw_Mode();

public:
    virtual void Reinitialize(const bool read_geometry=true); // Needs to be called after some state changes
protected:
    void Create_Geometry(const int id);
    void Update_Geometry(const int id);
    void Destroy_Geometry(const int id);
    virtual void Update_Object_Labels();
//#####################################################################
};

}
#endif
