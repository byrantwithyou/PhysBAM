//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace DEBUG_SUBSTEPS
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
namespace PhysBAM{
namespace DEBUG_SUBSTEPS
{
int write_substeps_level=-1;
std::function<void(const std::string&)> writer=0;
}
}
