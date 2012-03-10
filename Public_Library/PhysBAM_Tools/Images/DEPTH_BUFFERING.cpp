//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Images/DEPTH_BUFFERING.h>
using namespace PhysBAM;

//#####################################################################
template class DEPTH_BUFFERING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEPTH_BUFFERING<double>;
#endif
