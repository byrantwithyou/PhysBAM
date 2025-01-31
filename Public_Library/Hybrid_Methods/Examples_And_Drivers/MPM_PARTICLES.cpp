//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
MPM_PARTICLES()
    :store_Fp(false),store_B(false),store_C(false),store_S(false),store_lame(false),
    store_lame0(false),store_vort(false)
{
    this->Store_Velocity();
    this->Store_Mass();
    Add_Array("volume",&volume);
    Add_Array("F",&F);
    Add_Array("valid",&valid);
    this->template Add_Array<VECTOR<T,3> >("color");
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
    if(store) Add_Array("Fp",&Fp);
    else Remove_Array("Fp");
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
    if(store) Add_Array("B",&B);
    else Remove_Array("B");
}
//#####################################################################
// Function Store_B
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_C(bool store)
{
    if(store_C==store) return;
    store_C=store;
    if(store) Add_Array("C",&C);
    else Remove_Array("C");
}
//#####################################################################
// Function Store_S
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Store_S(bool store)
{
    if(store_S==store) return;
    store_S=store;
    if(store) Add_Array("S",&S);
    else Remove_Array("S");
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
        Add_Array("mu",&mu);
        Add_Array("lambda",&lambda);}
    else{
        Remove_Array("mu");
        Remove_Array("lambda");}
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
        Add_Array("mu0",&mu0);
        Add_Array("lambda0",&lambda0);}
    else{
        Remove_Array("mu0");
        Remove_Array("lambda0");}
}
//#####################################################################
// Function Store_Vort
//#####################################################################
template <class TV> void MPM_PARTICLES<TV>::
Store_Vort(bool store)
{
    if(store_vort==store) return;
    store_vort=store;
    if(store) Add_Array("vort",&vort);
    else Remove_Array("vort");
}
//#####################################################################
// Function Initialize_MPM_Particles
//#####################################################################
static int Initialize_MPM_Particles()
{
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
