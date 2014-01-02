//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Finite_Elements/INTERFACE_STOKES_MULTIGRID.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_MULTIGRID<TV>::
INTERFACE_STOKES_MULTIGRID(int num_levels,INTERFACE_STOKES_SYSTEM_COLOR<TV>* iss)
    :levels(num_levels)
{
    levels(0).iss=iss;
    for(int i=1;i<levels.m;i++) Construct_Level(i);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_STOKES_MULTIGRID<TV>::
~INTERFACE_STOKES_MULTIGRID()
{
    for(int i=1;i<levels.m;i++) delete levels(i).iss;
}
//#####################################################################
// Function Construct_Level
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Construct_Level(int l)
{
    PHYSBAM_FATAL_ERROR("TODO");
}
//#####################################################################
// Function Update
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Update()
{
    for(int i=0;i<levels.m;i++){
        LEVEL& l=levels(i);
        l.iss->Resize_Vector(l.tmp0);
        l.iss->Resize_Vector(l.tmp1);
        l.iss->Resize_Vector(l.tmp2);}
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Apply_Preconditioner(T_VECTOR& z,const T_VECTOR& x,bool initial_guess)
{
    if(initial_guess){
        levels(0).iss->Multiply(z,levels(0).tmp1);
        levels(0).tmp2.Copy(-1,levels(0).tmp1,x);}
    else levels(0).tmp2=x;

    for(int i=0;i<levels.m-1;i++){
        levels(i).tmp0*=0;
        levels(i).Interior_Smoother(levels(i).tmp0,levels(i).tmp2);
        levels(i).Boundary_Smoother(levels(i).tmp0,levels(i).tmp2);
        levels(i).Interior_Smoother(levels(i).tmp0,levels(i).tmp2);
        levels(i).iss->Multiply(levels(i).tmp0,levels(i).tmp1);
        levels(i).tmp0.Copy(-1,levels(i).tmp1,levels(i).tmp2);
        Coarsen(levels(i+1).tmp2,levels(i).tmp2,i);}

    Exact_Solve(levels.Last().tmp0,levels.Last().tmp2);

    for(int i=levels.m-2;i>=0;i--){
        Prolongation(levels(i).tmp1,levels(i+1).tmp0,i);
        levels(i).tmp0.Copy(1,levels(i).tmp1,levels(i).tmp0);
        levels(i).Interior_Smoother(levels(i).tmp0,levels(i).tmp2);
        levels(i).Boundary_Smoother(levels(i).tmp0,levels(i).tmp2);
        levels(i).Interior_Smoother(levels(i).tmp0,levels(i).tmp2);}

    z=levels(0).tmp0;
}
//#####################################################################
// Function Coarsen
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Coarsen(T_VECTOR& z,const T_VECTOR& x,int fine_level) const
{
    PHYSBAM_FATAL_ERROR("TODO");
}
//#####################################################################
// Function Prolongation
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Prolongation(T_VECTOR& z,const T_VECTOR& x,int fine_level) const
{
    PHYSBAM_FATAL_ERROR("TODO");
}
//#####################################################################
// Function Prolongation
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Exact_Solve(T_VECTOR& z,const T_VECTOR& x) const
{
    PHYSBAM_FATAL_ERROR("TODO");
}
//#####################################################################
// Function Interior_Smoother
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Interior_Smoother(T_VECTOR& z,const T_VECTOR& x) const
{
    PHYSBAM_FATAL_ERROR("TODO");
}
//#####################################################################
// Function Boundary_Smoother
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Boundary_Smoother(T_VECTOR& z,const T_VECTOR& x) const
{
    PHYSBAM_FATAL_ERROR("TODO");
}
namespace PhysBAM{
template class INTERFACE_STOKES_MULTIGRID<VECTOR<float,2> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<float,3> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<double,2> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<double,3> >;
}
