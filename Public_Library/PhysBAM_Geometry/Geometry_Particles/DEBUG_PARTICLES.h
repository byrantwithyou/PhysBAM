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
class DEBUG_PARTICLES
{
public:
    typedef typename TV::SCALAR T;
    DEBUG_PARTICLES();
    ~DEBUG_PARTICLES();

    GEOMETRY_PARTICLES<TV>& debug_particles;

    static GEOMETRY_PARTICLES<TV>* Store_Debug_Particles(GEOMETRY_PARTICLES<TV>* particle=0);
    void Write_Debug_Particles(STREAM_TYPE stream_type,const std::string& output_directory,int frame) const;
};
template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);
template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
}
#endif
