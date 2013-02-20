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
    GRID<TV> grid(count,box);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
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
    GRID<TV> grid(count,square_box);
    TV center=square_box.Center();
    T r=0.5*((-square_box.min_corner+square_box.max_corner)(0));
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
        TV x=grid.X(it.index);
        if((x-center).Magnitude()<=r) sample_X.Append(x);}
    this->Resize(sample_X.m);
    X=sample_X;
}
//#####################################################################
// Function Add_X_As_A_Grid
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Add_X_As_A_Grid(const VECTOR<int,TV::m>& count,const RANGE<TV>& box)
{
    GRID<TV> grid(count,box);
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
        TV x=grid.X(it.index);
        sample_X.Append(x);}
    old_X.Append_Elements(sample_X);
    this->Resize(old_X.m);
    X=old_X;
}
//#####################################################################
// Function Add_X_As_A_Ball
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Add_X_As_A_Ball(const VECTOR<int,TV::m>& count,const RANGE<TV>& square_box)
{
    GRID<TV> grid(count,square_box);
    TV center=square_box.Center();
    T r=0.5*((-square_box.min_corner+square_box.max_corner)(0));
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
        TV x=grid.X(it.index);
        if((x-center).Magnitude()<=r) sample_X.Append(x);}
    old_X.Append_Elements(sample_X);
    this->Resize(old_X.m);
    X=old_X;
}
//#####################################################################
// Function Reduce_X_As_A_Ball
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Reduce_X_As_A_Ball(const RANGE<TV>& square_box)
{
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    TV center=square_box.Center();
    T r=0.5*((-square_box.min_corner+square_box.max_corner)(0));
    ARRAY<TV> sample_X;
    for(int i=0;i<old_X.m;i++)
        if((old_X(i)-center).Magnitude()>r)
            sample_X.Append(old_X(i));
    this->Resize(sample_X.m);
    X=sample_X;
}
static int Initialize_MPM_Particles()
{
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
