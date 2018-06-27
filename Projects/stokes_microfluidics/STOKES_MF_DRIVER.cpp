//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/FINE_TIMER.h>
#include <Core/Log/SCOPE.h>
#include "STOKES_MF_EXAMPLE.h"
#include "STOKES_MF_DRIVER.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STOKES_MF_DRIVER<TV>::
STOKES_MF_DRIVER(STOKES_MF_EXAMPLE<TV>& example)
    :example(example)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STOKES_MF_DRIVER<TV>::
~STOKES_MF_DRIVER()
{
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void STOKES_MF_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void STOKES_MF_DRIVER<TV>::
Initialize()
{
    example.Initialize();
    example.Write_Output_Files(0);
}
//#####################################################################
namespace PhysBAM{
template class STOKES_MF_DRIVER<VECTOR<float,2> >;
template class STOKES_MF_DRIVER<VECTOR<float,3> >;
template class STOKES_MF_DRIVER<VECTOR<double,2> >;
template class STOKES_MF_DRIVER<VECTOR<double,3> >;
}
