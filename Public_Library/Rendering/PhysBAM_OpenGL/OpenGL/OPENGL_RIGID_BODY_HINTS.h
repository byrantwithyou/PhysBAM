//#####################################################################
// Copyright 2002, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_RIGID_BODY_HINTS
//##################################################################### 
//
//#####################################################################
// Guendelman - November 2, 2002
//#####################################################################
#ifndef __OPENGL_RIGID_BODY_HINTS__
#define __OPENGL_RIGID_BODY_HINTS__

#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
namespace PhysBAM
{

class OPENGL_RIGID_BODY_HINTS
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    OPENGL_MATERIAL material;
    bool            include_bounding_box;

    OPENGL_RIGID_BODY_HINTS()
        : material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR::White())),
          include_bounding_box(true)
    {}

    OPENGL_RIGID_BODY_HINTS(const OPENGL_MATERIAL &material,
                            bool include_bounding_box)
        : material(material), include_bounding_box(include_bounding_box)
    {}

    template<class RW> void Read(std::istream& input)
    {char version;Read_Binary<RW>(input,version);
    if(version!=1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized OPENGL_RIGID_BODY_HINTS version %d",version));
    Read_Binary<RW>(input,material,include_bounding_box);}

    template<class RW> void Write(std::ostream& output) const
    {char version=1;Write_Binary<RW>(output,version,material,include_bounding_box);}
};
}
#endif
