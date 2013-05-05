//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include "MPLE_DRIVER.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    MPLE_DRIVER<VECTOR<double,3> > a;


    return 0;
}
