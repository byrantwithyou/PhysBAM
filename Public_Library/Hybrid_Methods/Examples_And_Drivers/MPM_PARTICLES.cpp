//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
MPM_PARTICLES()
    :store_Fp(false),store_B(false),store_S(false),store_C(false),store_lame(false),store_lame0(false),
     store_phase(false)
{
    this->Store_Velocity();
    this->Store_Mass();
    Add_Array(ATTRIBUTE_ID_VOLUME,&volume);
    Add_Array(ATTRIBUTE_ID_F,&F);
    Add_Array(ATTRIBUTE_ID_VALID,&valid);
    this->template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
~MPM_PARTICLES()
{}
//#####################################################################
// Function Store_Fb
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_Fp(bool store)
{
    if(store_Fp==store) return;
    store_Fp=store;
    if(store) Add_Array(ATTRIBUTE_ID_FP,&Fp);
    else Remove_Array(ATTRIBUTE_ID_FP);
    Store_Lame0(store_lame && store_Fp);
}
//#####################################################################
// Function Store_B
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_B(bool store)
{
    if(store_B==store) return;
    store_B=store;
    if(store) Add_Array(ATTRIBUTE_ID_B,&B);
    else Remove_Array(ATTRIBUTE_ID_B);
}
//#####################################################################
// Function Store_S
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_S(bool store)
{
    if(store_S==store) return;
    store_S=store;
    if(store) Add_Array(ATTRIBUTE_ID_S,&S);
    else Remove_Array(ATTRIBUTE_ID_S);
}
//#####################################################################
// Function Store_C
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_C(bool store)
{
    if(store_C==store) return;
    store_C=store;
    if(store) Add_Array(ATTRIBUTE_ID_C,&C);
    else Remove_Array(ATTRIBUTE_ID_C);
}
//#####################################################################
// Function Store_Lame
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_Lame(bool store)
{
    if(store_lame==store) return;
    store_lame=store;
    if(store){
        Add_Array(ATTRIBUTE_ID_MU,&mu);
        Add_Array(ATTRIBUTE_ID_LAMBDA,&lambda);}
    else{
        Remove_Array(ATTRIBUTE_ID_MU);
        Remove_Array(ATTRIBUTE_ID_LAMBDA);}
    Store_Lame0(store_lame && store_Fp);
}
//#####################################################################
// Function Store_Lame0
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_Lame0(bool store)
{
    if(store_lame0==store) return;
    store_lame0=store;
    if(store){
        Add_Array(ATTRIBUTE_ID_MU0,&mu0);
        Add_Array(ATTRIBUTE_ID_LAMBDA0,&lambda0);}
    else{
        Remove_Array(ATTRIBUTE_ID_MU0);
        Remove_Array(ATTRIBUTE_ID_LAMBDA0);}
}
//#####################################################################
// Function Store_Phase
//#####################################################################
template <class TV> void MPM_PARTICLES<TV>::
Store_Phase(bool store)
{
    if(store_phase==store) return;
    store_phase=store;
    if(store) Add_Array(ATTRIBUTE_ID_PHASE,&phase);
    else Remove_Array(ATTRIBUTE_ID_PHASE);
}
//#####################################################################
// Function Initialize_MPM_Particles
//#####################################################################
static int Initialize_MPM_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_VOLUME,"volume");
    Register_Attribute_Name(ATTRIBUTE_ID_F,"F");
    Register_Attribute_Name(ATTRIBUTE_ID_FP,"Fp");
    Register_Attribute_Name(ATTRIBUTE_ID_B,"B");
    Register_Attribute_Name(ATTRIBUTE_ID_C,"C");
    Register_Attribute_Name(ATTRIBUTE_ID_VALID,"valid");
    Register_Attribute_Name(ATTRIBUTE_ID_S,"S");
    Register_Attribute_Name(ATTRIBUTE_ID_MU,"mu");
    Register_Attribute_Name(ATTRIBUTE_ID_LAMBDA,"lambda");
    Register_Attribute_Name(ATTRIBUTE_ID_MU0,"mu0");
    Register_Attribute_Name(ATTRIBUTE_ID_LAMBDA0,"lambda0");
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
