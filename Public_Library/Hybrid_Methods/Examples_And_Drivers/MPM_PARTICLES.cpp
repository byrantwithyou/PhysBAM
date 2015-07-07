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
    :store_B(false),store_S(false),store_C(false),volume(0,0),F(0,0),B(0,0),valid(0,0),one_over_lambda(0,0)
{
    this->Store_Velocity();
    this->Store_Mass();
    Add_Array(ATTRIBUTE_ID_VOLUME,&volume);
    Add_Array(ATTRIBUTE_ID_F,&F);
    Add_Array(ATTRIBUTE_ID_VALID,&valid);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
~MPM_PARTICLES()
{}
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
// Function Store_One_Over_Lambda
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_One_Over_Lambda(bool store)
{
    if(store_one_over_lambda==store) return;
    store_one_over_lambda=store;
    if(store) Add_Array(ATTRIBUTE_ID_ONE_OVER_LAMBDA,&one_over_lambda);
    else Remove_Array(ATTRIBUTE_ID_ONE_OVER_LAMBDA);
}
//#####################################################################
// Function Initialize_MPM_Particles
//#####################################################################
static int Initialize_MPM_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_VOLUME,"volume");
    Register_Attribute_Name(ATTRIBUTE_ID_F,"F");
    Register_Attribute_Name(ATTRIBUTE_ID_B,"B");
    Register_Attribute_Name(ATTRIBUTE_ID_C,"C");
    Register_Attribute_Name(ATTRIBUTE_ID_ONE_OVER_LAMBDA,"one_over_lambda");
    Register_Attribute_Name(ATTRIBUTE_ID_VALID,"valid");
    Register_Attribute_Name(ATTRIBUTE_ID_S,"S");
    return 0;
}
int initialize_deformables_particles=Initialize_MPM_Particles();
//#####################################################################
template class MPM_PARTICLES<VECTOR<float,2> >;
template class MPM_PARTICLES<VECTOR<float,3> >;
template class MPM_PARTICLES<VECTOR<double,2> >;
template class MPM_PARTICLES<VECTOR<double,3> >;
}
