//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_RAY_DEBUG
//#####################################################################
#ifndef __RENDERING_RAY_DEBUG__
#define __RENDERING_RAY_DEBUG__

#include <Core/Arrays/ARRAY.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Ray_Tracing/Rendering/RENDERING_RAY.h>
namespace PhysBAM{

template<class T> class PHOTON;
template<class T> class IRRADIANCE_SAMPLE;

template<class T>
class RENDERING_RAY_DEBUG:public NONCOPYABLE
{
public:
    RENDERING_RAY<T> ray;
    RENDERING_RAY_DEBUG* parent;
    ARRAY<RENDERING_RAY_DEBUG*> children;
    ARRAY<std::string> comments;
    ARRAY<PHOTON<T>*> photons_used;
    ARRAY<int> irradiance_cache_samples_used;
    bool hit_object;
    VECTOR<T,3> same_side_normal;

    RENDERING_RAY_DEBUG();
    RENDERING_RAY_DEBUG(const RENDERING_RAY<T>& ray_input);
    ~RENDERING_RAY_DEBUG();
    void Add_Child(RENDERING_RAY<T>& ray_to_add);
    void Add_Comment(const std::string& comment_string);
    void Print(std::ostream& output,const int spaces=0);
//#####################################################################
};
}
#endif
