//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "MPLE_DRIVER.h"

using namespace PhysBAM;

//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPLE_DRIVER<TV>::
MPLE_DRIVER()
{
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPLE_DRIVER<TV>::
~MPLE_DRIVER()
{
}

template<class TV> void MPLE_DRIVER<TV>::
Initialize()
{
    u.Resize(grid.Node_Indices(ghost),false);
}

