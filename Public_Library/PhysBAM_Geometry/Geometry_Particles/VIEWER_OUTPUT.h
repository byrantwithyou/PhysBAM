//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIEWER_OUTPUT
//#####################################################################
#ifndef __VIEWER_OUTPUT__
#define __VIEWER_OUTPUT__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;
template<class TV> class GRID;

template<class TV>
class VIEWER_OUTPUT
{
public:
    typedef typename TV::SCALAR T;

    int frame;
    std::string output_directory;
    STREAM_TYPE stream_type;

    const GRID<TV>& grid;
    DEBUG_PARTICLES<TV>& debug_particles;

    VIEWER_OUTPUT(STREAM_TYPE stream_type,const GRID<TV>& grid,const std::string& output_directory);
    ~VIEWER_OUTPUT();

    static VIEWER_OUTPUT* Singleton(VIEWER_OUTPUT* vo=0);

    void Flush_Frame(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const char* title);
    void Flush_Frame(const char* title);
};
template<class TV> inline void Flush_Frame(const char* title)
{
    VIEWER_OUTPUT<TV>::Singleton()->Flush_Frame(title);
}
template<class T,int d> inline void Flush_Frame(const ARRAY<T,FACE_INDEX<d> >& face_velocities,const char* title)
{
    VIEWER_OUTPUT<VECTOR<T,d> >::Singleton()->Flush_Frame(face_velocities,title);
}
}
#endif
