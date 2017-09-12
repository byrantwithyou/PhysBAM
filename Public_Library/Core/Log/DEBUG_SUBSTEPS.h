//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEBUG_SUBSTEPS
//#####################################################################
#ifndef __DEBUG_SUBSTEPS__
#define __DEBUG_SUBSTEPS__

#include <Core/Log/LOG.h>
#include <functional>
#include <string>
namespace PhysBAM{
namespace DEBUG_SUBSTEPS
{
extern int write_substeps_level;
extern std::function<void(const std::string&)> writer;
template<class ...Args> inline void
PHYSBAM_DEBUG_WRITE_SUBSTEP(const std::string& title,int level,Args... args)
{
    if(level<=write_substeps_level && writer)
        writer(LOG::sprintf(title.c_str(),args...));
}
inline void PHYSBAM_DEBUG_WRITE_SUBSTEP(const std::string& title,int level)
{
    if(level<=write_substeps_level && writer)
        writer(title);
}
}
using DEBUG_SUBSTEPS::PHYSBAM_DEBUG_WRITE_SUBSTEP;
}

#endif
