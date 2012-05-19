//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEBUG_PARTICLES
//#####################################################################
#ifndef __DEBUG_PARTICLES__
#define __DEBUG_PARTICLES__

#include <PhysBAM_Tools/Arrays/ATTRIBUTE_ID.h>
#include <PhysBAM_Tools/Read_Write/TYPED_STREAM.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;

template<class TV>
struct DEBUG_OBJECT
{
    typedef typename TV::SCALAR T;
    enum TYPE {segment=2,triangle=3} type;
    VECTOR<TV,3> X;
    VECTOR<T,3> color,bgcolor;
    bool draw_vertices;

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,X,type,color,bgcolor,draw_vertices);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,X,type,color,bgcolor,draw_vertices);}
};

template<class TV>
class DEBUG_PARTICLES
{
public:
    typedef typename TV::SCALAR T;
    DEBUG_PARTICLES();
    ~DEBUG_PARTICLES();

    GEOMETRY_PARTICLES<TV>& debug_particles;
    mutable ARRAY<DEBUG_OBJECT<TV> > debug_objects;

    static DEBUG_PARTICLES<TV>* Store_Debug_Particles(DEBUG_PARTICLES<TV>* particle=0);
    void Write_Debug_Particles(STREAM_TYPE stream_type,const std::string& output_directory,int frame) const;
};
template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);
template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
template<class TV> void Add_Debug_Object(const VECTOR<TV,2>& object,const VECTOR<typename TV::SCALAR,3>& color);
template<class TV> void Add_Debug_Object(const VECTOR<TV,3>& object,const VECTOR<typename TV::SCALAR,3>& color,const VECTOR<typename TV::SCALAR,3>& bgcolor);
}
#endif
