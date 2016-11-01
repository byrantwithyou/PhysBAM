//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_BASIC
//#####################################################################
#ifndef __OPENGL_COMPONENT_BASIC__
#define __OPENGL_COMPONENT_BASIC__

#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class T,class T2>
class OPENGL_COMPONENT_BASIC:public OPENGL_COMPONENT<T>
{
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;
    T2  &object;
    bool own_object;
    bool use_clip_planes;

    OPENGL_COMPONENT_BASIC(STREAM_TYPE stream_type,T2 &object)
        :OPENGL_COMPONENT<T>(stream_type),object(object),own_object(true)
    {
        Use_Clip_Planes(false);
    }

    virtual ~OPENGL_COMPONENT_BASIC()
    {
        if(own_object) delete &object;
    }

    void Use_Clip_Planes(bool use=true) 
    {use_clip_planes=use;}

    void Display() const override
    {if(draw){
        if(use_clip_planes){glEnable(GL_CLIP_PLANE0);glEnable(GL_CLIP_PLANE1);}
        object.Display();
        if(use_clip_planes){glDisable(GL_CLIP_PLANE0);glDisable(GL_CLIP_PLANE1);}}}

    bool Use_Bounding_Box() const  override
    {if(draw) return object.Use_Bounding_Box();else return false;}

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override
    {if(draw) return object.Bounding_Box();return RANGE<VECTOR<T,3> >::Centered_Box();}

    bool Is_Transparent() const override
    {return object.Is_Transparent();}

    virtual void Turn_Smooth_Shading_Off() override
    {object.Turn_Smooth_Shading_Off();}

    virtual void Turn_Smooth_Shading_On() override
    {object.Turn_Smooth_Shading_On();}

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override
    {return object.Get_Selection_Priority(indices);}

    bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override
    {return object.Set_Selection(indices,modifiers);}

    virtual void Clear_Selection() override
    {object.Clear_Selection();}

    virtual void Set_Slice(OPENGL_SLICE *slice_input) override
    {slice=slice_input;object.Set_Slice(slice_input);}

    virtual void Slice_Has_Changed() override
    {object.Slice_Has_Changed();}

    void Print_Selection_Info(std::ostream& ostream) const override
    {object.Print_Selection_Info(ostream);}
};
}

#endif
