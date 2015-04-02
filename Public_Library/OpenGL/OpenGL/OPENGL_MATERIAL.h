//#####################################################################
// Copyright 2002-2008, Robert Bridson, Geoffrey Irving, Neil Molino, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_MATERIAL
//##################################################################### 
#ifndef __OPENGL_MATERIAL__
#define __OPENGL_MATERIAL__

#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Tools/Utilities/EXCEPTIONS.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
namespace PhysBAM{

class OPENGL_MATERIAL
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    OPENGL_COLOR ambient, diffuse, specular;
    float shininess;

    OPENGL_MATERIAL()
        : ambient(OPENGL_COLOR::Gray()),diffuse(OPENGL_COLOR::Gray()),specular(OPENGL_COLOR::Gray()),shininess(45)
    {}

    OPENGL_MATERIAL(const OPENGL_COLOR& color) // default is a diffuse metallic
        : ambient(color),diffuse(color),specular(color),shininess(45)
    {}

    OPENGL_MATERIAL(const OPENGL_COLOR& ambient_input,const OPENGL_COLOR& diffuse_input,
                    const OPENGL_COLOR& specular_input,const float shininess_input=45)
        : ambient(ambient_input),diffuse(diffuse_input),specular(specular_input),shininess(shininess_input)
    {}
    
    OPENGL_MATERIAL Grayscale() const
    {return OPENGL_MATERIAL(ambient.Grayscale(),diffuse.Grayscale(),specular.Grayscale(),shininess);}

    void Send_To_GL_Pipeline(const GLenum face=GL_FRONT_AND_BACK) const
    {
        glMaterialfv(face,GL_AMBIENT,ambient.rgba);
        glMaterialfv(face,GL_DIFFUSE,diffuse.rgba);
        glMaterialfv(face,GL_SPECULAR,specular.rgba);
        glMaterialf(face,GL_SHININESS,shininess);
    }

    static OPENGL_MATERIAL Matte(const OPENGL_COLOR& color)
    {return OPENGL_MATERIAL(color,color,OPENGL_COLOR::Black(),60);}

    static OPENGL_MATERIAL Metal(const OPENGL_COLOR& color)
    {OPENGL_COLOR diff=OPENGL_COLOR(.3f*color.rgba[0],.3f*color.rgba[1],.3f*color.rgba[2],color.rgba[3]);
    return OPENGL_MATERIAL(color,diff,color,30);}
    
    static OPENGL_MATERIAL Plastic(const OPENGL_COLOR& color)
    {return OPENGL_MATERIAL(color,color,OPENGL_COLOR(1,1,1,color.rgba[3]),60);}

    static OPENGL_MATERIAL Random()
    {return OPENGL_MATERIAL(OPENGL_COLOR::Random(),OPENGL_COLOR::Random(),OPENGL_COLOR::Random());}

    static void Set_Global_Ambient(const float y)//by default (.2,.2,.2,1) in OpenGL
    {float global_ambient[4]={y,y,y,1}; glLightModelfv(GL_LIGHT_MODEL_AMBIENT,global_ambient);}
    
    static OPENGL_MATERIAL Affine_Combination(float t,const OPENGL_MATERIAL& material1,const OPENGL_MATERIAL& material2)
    {return OPENGL_MATERIAL(OPENGL_COLOR::Affine_Combination(t,material1.ambient,material2.ambient),OPENGL_COLOR::Affine_Combination(t,material1.diffuse,material2.diffuse),
        OPENGL_COLOR::Affine_Combination(t,material1.specular,material2.specular),material1.shininess+t*(material2.shininess-material1.shininess));}

    template<class RW> void Read(std::istream& input)
    {char version;Read_Binary<RW>(input,version);
    if(version!=1) throw READ_ERROR(LOG::sprintf("Unrecognized OPENGL_MATERIAL version %d",version));
    Read_Binary<RW>(input,ambient,diffuse,specular,shininess);}

    template<class RW> void Write(std::ostream& output) const
    {char version=1;Write_Binary<RW>(output,version,ambient,diffuse,specular,shininess);}
//########################################################################
};
}
#endif

