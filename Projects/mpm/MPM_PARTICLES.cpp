//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "MPM_PARTICLES.h"
#include "MPM_PARTICLES_FORWARD.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
MPM_PARTICLES()
{
    Store_Velocity();
    Store_Mass();
    Add_Array(ATTRIBUTE_ID_MATERIAL_COORDINATE,&Xm);
    Add_Array(ATTRIBUTE_ID_VOLUME,&volume);
    Add_Array(ATTRIBUTE_ID_FE,&Fe);
    Add_Array(ATTRIBUTE_ID_FP,&Fp);
    Add_Array(ATTRIBUTE_ID_PARTICLE_DOMAIN,&particle_domain);
    Add_Array(ATTRIBUTE_ID_MU,&mu);
    Add_Array(ATTRIBUTE_ID_LAMBDA,&lambda);
    Add_Array(ATTRIBUTE_ID_MU0,&mu0);
    Add_Array(ATTRIBUTE_ID_LAMBDA0,&lambda0);
    Add_Array(ATTRIBUTE_ID_PRESSURE,&pressure);
    Add_Array(ATTRIBUTE_ID_USE_PLASTICITY_YIELD,&use_plasticity_yield);
    Add_Array(ATTRIBUTE_ID_USE_PLASTICITY_CLAMP,&use_plasticity_clamp);
    Add_Array(ATTRIBUTE_ID_YIELD_MIN,&yield_min);
    Add_Array(ATTRIBUTE_ID_YIELD_MAX,&yield_max);
    Add_Array(ATTRIBUTE_ID_CLAMP_MIN,&clamp_min);
    Add_Array(ATTRIBUTE_ID_CLAMP_MAX,&clamp_max);
    Add_Array(ATTRIBUTE_ID_USE_VISCO_PLASTICITY,&use_visco_plasticity);
    Add_Array(ATTRIBUTE_ID_VISCO_NU,&visco_nu);
    Add_Array(ATTRIBUTE_ID_VISCO_TAU,&visco_tau);
    Add_Array(ATTRIBUTE_ID_VISCO_KAPPA,&visco_kappa);
    Add_Array(ATTRIBUTE_ID_COMPRESS,&compress);
    rand_generator.Set_Seed(0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
~MPM_PARTICLES()
{}
//#####################################################################
// Function Add_Randomly_Sampled_Object
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Add_Randomly_Sampled_Implicit_Object(const IMPLICIT_OBJECT<TV>& object,const T exclude_radius)
{
    const_cast<IMPLICIT_OBJECT<TV>&>(object).Update_Box();
    RANGE<TV> bounding_box(object.box);
    TV_INT cells(ceil(bounding_box.Edge_Lengths()/exclude_radius));
    bounding_box.Scale_About_Center(exclude_radius*TV(cells)/bounding_box.Edge_Lengths());
    GRID<TV> grid(cells,bounding_box,true);
    ARRAY<bool,TV_INT> ok(grid.Domain_Indices(1));
    ARRAY<TV_INT> todo;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next())
        if(object.Extended_Phi(it.Location())<=0){
            ok(it.index)=true;
            todo.Append(it.index);}
    rand_generator.Random_Shuffle(todo);
    RANGE<TV_INT> neigh_index(TV_INT()-1,TV_INT()+2);
    for(int i=0;i<todo.m;i++){
        TV_INT index=todo(i);
        if(!ok(index)) continue;
        TV new_X;
        for(int j=0;j<1000;j++){
            new_X=rand_generator.Get_Uniform_Vector(grid.Cell_Domain(index));
            if(object.Extended_Phi(new_X)<=0) break;}
        for(RANGE_ITERATOR<TV::m> it(neigh_index);it.Valid();it.Next())
            ok(it.index+index)=false;
        int p=this->Add_Element();
        X(p)=new_X;
        Xm(p)=new_X;}
}
//#####################################################################
// Function Set_Material_Properties
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Set_Material_Properties(int start_index,int count,T mass_in,T mu_in,T lambda_in,bool compress_in,bool pressure_in)
{
    for(int i=start_index;i<start_index+count;i++){
        mass(i)=mass_in;
        mu(i)=mu_in;mu0(i)=mu_in;
        lambda(i)=lambda_in;lambda0(i)=lambda_in;
        compress(i)=compress_in;
        pressure(i)=pressure_in;}
}
//#####################################################################
// Function Set_Initial_State
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Set_Initial_State(int start_index,int count,MATRIX<T,TV::m> Fe_in,MATRIX<T,TV::m> Fp_in,TV V_in)
{
    for(int i=start_index;i<start_index+count;i++){
        Fe(i)=Fe_in;
        Fp(i)=Fp_in;
        V(i)=V_in;}
}
//#####################################################################
// Function Set_Plasticity
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Set_Plasticity(int start_index,int count,bool use_plasticity_yield_in,T yield_min_in,T yield_max_in,bool use_plasticity_clamp_in,T clamp_min_in,T clamp_max_in)
{
    for(int i=start_index;i<start_index+count;i++){
        use_plasticity_yield(i)=use_plasticity_yield_in;
        yield_min(i)=yield_min_in;
        yield_max(i)=yield_max_in;
        use_plasticity_clamp(i)=use_plasticity_clamp_in;
        clamp_min(i)=clamp_min_in;
        clamp_max(i)=clamp_max_in;}
}
//#####################################################################
// Function Set_Visco_Plasticity
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Set_Visco_Plasticity(int start_index,int count,bool use_visco_plasticity_in,T visco_nu_in,T visco_tau_in,T visco_kappa_in)
{
    for(int i=start_index;i<start_index+count;i++){
        use_visco_plasticity(i)=use_visco_plasticity_in;
        visco_nu(i)=visco_nu_in;
        visco_tau(i)=visco_tau_in;
        visco_kappa(i)=visco_kappa_in;}
}
static int Initialize_MPM_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_MATERIAL_COORDINATE,"Xm");
    Register_Attribute_Name(ATTRIBUTE_ID_VOLUME,"volume");
    Register_Attribute_Name(ATTRIBUTE_ID_FE,"Fe");
    Register_Attribute_Name(ATTRIBUTE_ID_FP,"Fp");
    Register_Attribute_Name(ATTRIBUTE_ID_PARTICLE_DOMAIN,"particle_domain");
    Register_Attribute_Name(ATTRIBUTE_ID_MU,"mu");
    Register_Attribute_Name(ATTRIBUTE_ID_LAMBDA,"lambda");
    Register_Attribute_Name(ATTRIBUTE_ID_MU0,"mu0");
    Register_Attribute_Name(ATTRIBUTE_ID_LAMBDA0,"lambda0");
    Register_Attribute_Name(ATTRIBUTE_ID_PRESSURE,"pressure");
    Register_Attribute_Name(ATTRIBUTE_ID_USE_PLASTICITY_YIELD,"use_plasticity_yield");
    Register_Attribute_Name(ATTRIBUTE_ID_USE_PLASTICITY_CLAMP,"use_plasticity_clamp");
    Register_Attribute_Name(ATTRIBUTE_ID_YIELD_MIN,"yield_min");
    Register_Attribute_Name(ATTRIBUTE_ID_YIELD_MAX,"yield_max");
    Register_Attribute_Name(ATTRIBUTE_ID_CLAMP_MIN,"clamp_min");
    Register_Attribute_Name(ATTRIBUTE_ID_CLAMP_MAX,"clamp_max");
    Register_Attribute_Name(ATTRIBUTE_ID_USE_VISCO_PLASTICITY,"use_visco_plasticity");
    Register_Attribute_Name(ATTRIBUTE_ID_VISCO_NU,"visco_nu");
    Register_Attribute_Name(ATTRIBUTE_ID_VISCO_TAU,"visco_tau");
    Register_Attribute_Name(ATTRIBUTE_ID_VISCO_KAPPA,"visco_kappa");
    Register_Attribute_Name(ATTRIBUTE_ID_COMPRESS,"compress");
    return 0;
}
int initialize_mpm_particles=Initialize_MPM_Particles();
//#####################################################################
template class MPM_PARTICLES<VECTOR<float,2> >;
template class MPM_PARTICLES<VECTOR<float,3> >;
template class MPM_PARTICLES<VECTOR<double,2> >;
template class MPM_PARTICLES<VECTOR<double,3> >;
}
