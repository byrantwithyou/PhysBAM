//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "MPM_PARTICLES.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int IN> MPM_PARTICLES<TV,IN>::
MPM_PARTICLES()
{
    Store_Velocity();
    Store_Mass();
    Add_Array(ATTRIBUTE_ID_VOLUME,&volume);
    Add_Array(ATTRIBUTE_ID_FE,&Fe);
    Add_Array(ATTRIBUTE_ID_FP,&Fp);
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV,int IN> MPM_PARTICLES<TV,IN>::
~MPM_PARTICLES()
{}
static int Initialize_MPM_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_VOLUME,"volume");
    Register_Attribute_Name(ATTRIBUTE_ID_FE,"Fe");
    Register_Attribute_Name(ATTRIBUTE_ID_FP,"Fp");
    return 0;
}
int initialize_mpm_particles=Initialize_MPM_Particles();
//#####################################################################
template class MPM_PARTICLES<VECTOR<float,1> >;
template class MPM_PARTICLES<VECTOR<float,2> >;
template class MPM_PARTICLES<VECTOR<float,3> >;
template class MPM_PARTICLES<VECTOR<double,1> >;
template class MPM_PARTICLES<VECTOR<double,2> >;
template class MPM_PARTICLES<VECTOR<double,3> >;
}
