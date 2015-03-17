//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D__
#define __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/PAIR.h>
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

    std::string basedir;
    int frame_loaded;
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
    OPENGL_SELECTION<T>* current_selection;
    bool need_destroy_rigid_body_collection;

public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(STREAM_TYPE stream_type,const std::string& basedir);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(STREAM_TYPE stream_type,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir);
    virtual ~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D();
    
//#####################################################################
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display() const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Set_Draw_Object(int i, bool draw_it);  // Need to call Reinitialize after changing draw objects
    bool Get_Draw_Object(int i) const;
    void Set_Object_Color(int i, const OPENGL_COLOR &color);
    void Set_Use_Object_Bounding_Box(int i, bool use_it);

    void Toggle_Output_Positions();
    void Toggle_Show_Object_Names();
    void Toggle_Draw_Mode();

public:
    virtual void Reinitialize(const bool force=false,const bool read_geometry=true); // Needs to be called after some state changes
protected:
    void Create_Geometry(const int id);
    void Update_Geometry(const int id);
    void Destroy_Geometry(const int id);
    virtual void Update_Object_Labels();
//#####################################################################
};

}
#endif
