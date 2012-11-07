//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Duc Nguyen, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RUNGEKUTTA
//#####################################################################
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
namespace PhysBAM{
//#####################################################################
// Function Start
//#####################################################################
//#####################################################################
template class RUNGEKUTTA_CORE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RUNGEKUTTA_CORE<double>;
#endif
}
