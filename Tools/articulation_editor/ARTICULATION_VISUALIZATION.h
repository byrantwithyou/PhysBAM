//#####################################################################
// Copyright 2005, Jared Go, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATION_VISUALIZATION
//##################################################################### 
#ifndef __ARTICULATION_VISUALIZATION__
#define __ARTICULATION_VISUALIZATION__

#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
#include "OPENGL_COMPONENT_MUSCLE_3D.h"


namespace PhysBAM
{

template<class T> class OPENGL_CUSTOM_MOUSE_HANDLER;


template<class T> 
class ARTICULATION_VISUALIZATION : public ANIMATED_VISUALIZATION
{
public:
    ARTICULATION_VISUALIZATION();
    ~ARTICULATION_VISUALIZATION();

    // selection call (bad abstraction, but ok for now)
    void Handle_Object_Selected(OPENGL_COMPONENT_MUSCLE_3D<T>*m,int viapt);

protected:
    virtual void Add_Arguments(PARSE_ARGS &parse_args);
    virtual void Parse_Arguments(PARSE_ARGS &parse_args);
    virtual void Initialize_Components();
    virtual void Add_OpenGL_Initialization();
    virtual void Add_Key_Bindings();
    virtual void Update_OpenGL_Strings();

    // commands
    void Save_Muscles_Command();
    DEFINE_CALLBACK_CREATOR(ARTICULATION_VISUALIZATION, Save_Muscles_Command);
    void Change_Segment_Type_Command();
    DEFINE_CALLBACK_CREATOR(ARTICULATION_VISUALIZATION, Change_Segment_Type_Command);
    void Change_Segment_Type_Prompt();
    DEFINE_CALLBACK_CREATOR(ARTICULATION_VISUALIZATION, Change_Segment_Type_Prompt);
    bool Valid_Segment_Type(int segment_type);
    void Set_Segment_Type(int segment_type);

    void Change_Curve_Thickness_Prompt();
    DEFINE_CALLBACK_CREATOR(ARTICULATION_VISUALIZATION, Change_Curve_Thickness_Prompt);
    void Change_Curve_Offset_Thickness_Prompt();
    DEFINE_CALLBACK_CREATOR(ARTICULATION_VISUALIZATION, Change_Curve_Offset_Thickness_Prompt);
    void Change_Tendon_Fraction_1_Prompt();
    DEFINE_CALLBACK_CREATOR(ARTICULATION_VISUALIZATION, Change_Tendon_Fraction_1_Prompt);
    void Change_Tendon_Fraction_2_Prompt();
    DEFINE_CALLBACK_CREATOR(ARTICULATION_VISUALIZATION, Change_Tendon_Fraction_2_Prompt);

    // list loading
    void Load_List(std::string filename, ARRAY<std::string> &list);

    // data loading
    void Generate_Bones(std::string listfile);
    void Generate_Muscles(std::string listfile);    
    void Generate_Joints(std::string listfile); 

    // save back out
    void Save_Bones(std::string listfile);
    void Save_Muscles();
    void Save_Joints(std::string listfile);

    // parsing a single joint
    JOINT_3D<T>* Parse_Joint(std::istream &s, std::string &parent, std::string &child);

private:
   
    // our mouse handler
    OPENGL_CUSTOM_MOUSE_HANDLER<T> *mouse_handler;

    // currently highlighted muscle
    OPENGL_COMPONENT_MUSCLE_3D<T>* currmuscle;

    // data path of project
    std::string bone_path;
    std::string muscle_path;

    // skeleton data - named lists of muscles and bones we are showing
    ARRAY<std::string> bone_list;
    ARRAY<std::string> muscle_list;

    // map of loaded objects by name
    HASHTABLE<std::string,OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>*> bone_map;
    HASHTABLE<std::string,OPENGL_COMPONENT_MUSCLE_3D<T>*> muscle_map;

    // components
    ARRAY<OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>*> triangulated_surface_component;    
    ARRAY<std::string> filenames;
    ARRAY<FRAME_3D<T> > frames;    
};

}

#endif
