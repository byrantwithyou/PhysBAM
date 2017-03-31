//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER  
//##################################################################### 
//
// Inherited by EULER_1D, EULER_2D, and EULER_3D.
// Initializes the equation of state information.
//
//#####################################################################
#ifndef __EULER__
#define __EULER__    

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Compressible/Euler_Equations/BOUNDARY_OBJECT_EULER.h>
#include <Compressible/Euler_Equations/BOUNDARY_OBJECT_SOLID_VELOCITY.h>
namespace PhysBAM{

template<class TV>
class EULER
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    enum {d=TV::m+2};
    typedef VECTOR<T,d> TV_DIMENSION;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;typedef ARRAYS_ND_BASE<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_BASE;
public:
    EOS<T>* eos;
    BOUNDARY<TV,TV_DIMENSION>* boundary;
    CONSERVATION<TV,d>* conservation;
    T cfl_number;
    VECTOR<bool,2*TV::m> open_boundaries;
    bool use_solid_velocity_in_ghost_cells;

    bool cut_out_grid; // (1) cut out grid, (0) no cut out grid 
    T max_time_step;
    bool use_force;
    TV gravity;
    BOUNDARY<TV,TV_DIMENSION> boundary_default;
    CONSERVATION_ENO_LLF<TV,d> conservation_default;
    EOS_GAMMA<T> eos_default;

protected:
    EULER();
    virtual ~EULER();

public:
    static TV Get_Velocity(const TV_DIMENSION& u)
    {return u.template Slice<1,TV::m>()/u(0);}

    static T Get_Velocity_Component(const TV_DIMENSION& U,const int axis)
    {assert((unsigned)axis<(unsigned)TV::m);return U(axis+1)/U(0);}

    static T Get_Density(const T_ARRAYS_DIMENSION_BASE& U,const TV_INT& cell)
    {return U(cell)(0);}

    static T Get_Total_Energy(const T_ARRAYS_DIMENSION_BASE& U,const TV_INT& cell)
    {return U(cell)(TV::m+1);}

    static void Set_Euler_State_From_rho_velocity_And_internal_energy(T_ARRAYS_DIMENSION_BASE& U,const TV_INT& cell,const T rho,const TV& velocity,const T e)
    {U(cell)(0)=rho;
    for(int k=0;k<TV::m;k++) U(cell)(k+1)=rho*velocity[k];
    U(cell)(TV::m+1)=rho*(e+velocity.Magnitude_Squared()*(T).5);}

    static VECTOR<T,TV::m+2> Get_Euler_State_From_rho_velocity_And_internal_energy(const T rho,const TV& velocity,const T e)
    {VECTOR<T,TV::m+2> U;
     U(0)=rho;
     for(int k=0;k<TV::m;k++) U(k+1)=rho*velocity[k];
     U(TV::m+1)=rho*(e+(T).5*velocity.Magnitude_Squared());
     return U;}

    static T e(const TV_DIMENSION& u)
    {return u(d-1)/u(0)-Get_Velocity(u).Magnitude_Squared()/2;}

    static T p(EOS<T>* eos_input,const TV_DIMENSION& u)
    {return eos_input->p(u(0),e(u));}

    static T enthalpy(const EOS<T>& eos,const VECTOR<T,TV::m+2>& u)
    {T internal_energy=e(u);T rho=u(0);T p=eos.p(rho,internal_energy);return internal_energy+p/rho;}

    T enthalpy(const VECTOR<T,TV::m+2>& u) const
    {T internal_energy=e(u);T rho=u(0);T p=eos->p(rho,internal_energy);return internal_energy+p/rho;}

    void Set_Custom_Boundary(BOUNDARY<TV,TV_DIMENSION>& boundary_input)
    {boundary=&boundary_input;
    for(int axis=0;axis<TV::m;axis++){
        open_boundaries(2*axis)=boundary->Constant_Extrapolation(2*axis);
        open_boundaries(2*axis+1)=boundary->Constant_Extrapolation(2*axis+1);}}
    
    void Set_Custom_Conservation(CONSERVATION<TV,d>& conservation_input)
    {conservation=&conservation_input;
        if(use_solid_velocity_in_ghost_cells) conservation->Set_Custom_Object_Boundary(*new BOUNDARY_OBJECT_SOLID_VELOCITY<TV>);
        else conservation->Set_Custom_Object_Boundary(*new BOUNDARY_OBJECT_EULER<TV>);}

    void Set_Custom_Equation_Of_State(EOS<T>& eos_input)
    {eos=&eos_input;}

    void Set_Max_Time_Step(const T max_time_step_input=1e8)
    {max_time_step=max_time_step_input;}

    void Set_CFL_Number(const T cfl_number_input=.5)
    {cfl_number=cfl_number_input;}

//#####################################################################
    virtual void Log_Parameters() const;
//#####################################################################
};
//#####################################################################
}
#endif
