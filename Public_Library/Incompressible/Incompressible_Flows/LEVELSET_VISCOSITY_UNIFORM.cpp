//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM_SYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_VISCOSITY_UNIFORM<TV>::
LEVELSET_VISCOSITY_UNIFORM(BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input,const GRID<TV>& grid_input,T dt,T density,T viscosity)
    :index_map(grid_input,callback_input),system(index_map),scale(dt*viscosity/density),print_matrix(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_VISCOSITY_UNIFORM<TV>::
~LEVELSET_VISCOSITY_UNIFORM()
{
    vectors.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Apply_Full_Viscosity
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Apply_Full_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,bool fully_explicit,bool fully_implicit,int axis)
{
    T local_scale=(fully_explicit || fully_implicit)?scale:scale/2;
    index_map.Compute(axis,periodic_boundary);
    LOG::cout<<"Solve size "<<index_map.index_to_face.m<<std::endl;
    system.Compute(axis,local_scale);
    Resize_Vectors(fully_explicit);
    index_map.Gather(u,b.v);

    static int solve_id=-1;solve_id++;
    if(print_matrix){
        LOG::cout<<"viscosity solve id "<<solve_id<<std::endl;
        OCTAVE_OUTPUT<T>(LOG::sprintf("visc-M-%i.txt",solve_id).c_str()).Write("M",system,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>(LOG::sprintf("visc-b-%i.txt",solve_id).c_str()).Write("b",b);}

    if(!fully_implicit) Apply_Explicit_Viscosity(u,axis);
    if(!fully_implicit && !fully_explicit) x.v.Exchange(b.v);
    if(!fully_explicit) Apply_Implicit_Viscosity(u,axis);
    if(print_matrix) OCTAVE_OUTPUT<T>(LOG::sprintf("visc-x-%i.txt",solve_id).c_str()).Write("x",x);
    index_map.Scatter(x.v,u);
}
//#####################################################################
// Function Apply_Viscosity
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Apply_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,bool fully_explicit,bool fully_implicit,bool coupled)
{
    if(coupled) Apply_Full_Viscosity(u,fully_explicit,fully_implicit,0);
    else for(int a=0;a<d;a++) Apply_Full_Viscosity(u,fully_explicit,fully_implicit,a);
}
//#####################################################################
// Function Apply_Implicit_Viscosity
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Apply_Implicit_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,int axis)
{
    system.Add_Constant_Part(b);
    x.v=b.v;
    CONJUGATE_GRADIENT<T> cg;
    cg.print_diagnostics=true;
    bool result=cg.Solve(system,x,b,vectors,(T)1e-4,1,1000);
    PHYSBAM_ASSERT(result);
}
//#####################################################################
// Function Apply_Explicit_Viscosity
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Apply_Explicit_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,int axis)
{
    system.poisson.P.Times(x.v,b.v);
    system.Add_Constant_Part(x);
    b.v=x.v+scale*b.v;
    x.v=b.v;
}
//#####################################################################
// Function Resize_Vectors
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Resize_Vectors(bool minimal)
{
    x.v.Resize(index_map.index_to_face.m);
    b.v.Resize(index_map.index_to_face.m);
    if(print_matrix) KRYLOV_SOLVER<T>::Ensure_Size(vectors,x,2);
}
namespace PhysBAM{
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<float,1> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<float,2> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<float,3> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<double,1> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<double,2> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<double,3> >;
}
