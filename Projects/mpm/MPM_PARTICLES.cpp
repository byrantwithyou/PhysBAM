//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
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
    Add_Array(ATTRIBUTE_ID_ISDIRICHLET,&is_dirichlet);
    Add_Array(ATTRIBUTE_ID_DENSITY,&density);
    Add_Array(ATTRIBUTE_ID_VOLUME,&volume);
    Add_Array(ATTRIBUTE_ID_FE,&Fe);
    Add_Array(ATTRIBUTE_ID_FP,&Fp);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
~MPM_PARTICLES()
{}
//#####################################################################
// Function Initialize_X_As_A_Grid
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Initialize_X_As_A_Grid(const VECTOR<int,TV::m>& count,const RANGE<TV>& box)
{
    typedef VECTOR<int,TV::m> TV_INT;
    GRID<TV> grid(count,box);
    RANGE<TV_INT> range(TV_INT(),TV_INT()+count);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        TV x=grid.X(it.index);
        sample_X.Append(x);}
    this->Resize(sample_X.m);
    X=sample_X;
}
//#####################################################################
// Function Initialize_X_As_A_Ball
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Initialize_X_As_A_Ball(const VECTOR<int,TV::m>& count,const RANGE<TV>& square_box)
{
    typedef VECTOR<int,TV::m> TV_INT;
    GRID<TV> grid(count,square_box);
    TV center=square_box.Center();
    T r=0.5*((-square_box.min_corner+square_box.max_corner)(0));
    RANGE<TV_INT> range(TV_INT(),TV_INT()+count);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        TV x=grid.X(it.index);
        if((x-center).Magnitude()<=r) sample_X.Append(x);}
    this->Resize(sample_X.m);
    X=sample_X;
}
//#####################################################################
// Function Identify_Dirichlet_Particles_With_A_Box
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Identify_Dirichlet_Particles_With_A_Box(const RANGE<TV>& box)
{
    for(int p=0;p<this->number;p++) if(box.Lazy_Inside(X(p))) is_dirichlet(p)=true;
}
static int Initialize_MPM_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_DENSITY,"density");
    Register_Attribute_Name(ATTRIBUTE_ID_VOLUME,"volume");
    Register_Attribute_Name(ATTRIBUTE_ID_FE,"Fe");
    Register_Attribute_Name(ATTRIBUTE_ID_FP,"Fp");
    return 0;
}
int initialize_mpm_particles=Initialize_MPM_Particles();
//#####################################################################
template class MPM_PARTICLES<VECTOR<float,2> >;
template class MPM_PARTICLES<VECTOR<float,3> >;
template class MPM_PARTICLES<VECTOR<double,2> >;
template class MPM_PARTICLES<VECTOR<double,3> >;
}
