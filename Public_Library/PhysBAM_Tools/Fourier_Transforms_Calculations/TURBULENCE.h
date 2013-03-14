//#####################################################################
// Copyright 2002, Ron Fedkiw, Eran Guendelman, Robert Bridson
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TURBULENCE
//#####################################################################
//
// Generates a periodic turbulent field from a random complex sequence.
//
//#####################################################################
#ifndef __TURBULENCE__
#define __TURBULENCE__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class TURBULENCE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    RANDOM_NUMBERS<T> *random,random_default;
    int incompressible; // (1) incompressible, (0) compressible
    T k_inertial; // sets the scale of the turbulence
    T epsilon;    // dissipation
    T constant;  // Kolmogorov constant
    T rescaled_average_velocity; // rescales the computed velocity to have this as the average velocity

public:
    T time_start,time_end;
    GRID<TV> grid;
    VECTOR<ARRAY<T,TV_INT>,TV::m> u_old,u_new;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;

    TURBULENCE()
    {
        random=&random_default;
        Enforce_Incompressibility();
        Set_Lowest_Angular_Frequency();
        Set_Dissipation();
        Set_Constant();
        Set_Rescaled_Average_Velocity();
    }

    TURBULENCE(const GRID<TV>& grid_input)
    {
        random=&random_default;
        Enforce_Incompressibility();
        Set_Lowest_Angular_Frequency();
        Set_Dissipation();
        Set_Constant();
        Set_Rescaled_Average_Velocity();
        Initialize_Grid(grid_input);
    }

    void Set_Custom_Random(RANDOM_NUMBERS<T>& random_input)
    {random=&random_input;}

    void Enforce_Incompressibility()
    {incompressible=1;}

    void Do_Not_Enforce_Incompressibility()
    {incompressible=0;}

    void Set_Lowest_Angular_Frequency(const T k_inertial_input=4) // k_inertial=2*pi*|frequency|
    {k_inertial=k_inertial_input;}

    void Set_Dissipation(const T epsilon_input=1)
    {epsilon=epsilon_input;}

    void Set_Constant(const T constant_input=1.5)
    {constant=constant_input;}

    void Set_Rescaled_Average_Velocity(const T velocity_magnitude=1)
    {rescaled_average_velocity=velocity_magnitude;}

    void Initialize_Grid(const GRID<TV>& grid_input)
    {grid=grid_input;for(int i=0;i<TV::m;i++){u_old(i).Resize(grid.Domain_Indices());u_new(i).Resize(grid.Domain_Indices());}}

    void Generate_Initial_Turbulence(const T time_start_input=0,const T time_end_input=1)
    {time_start=time_start_input;time_end=time_end_input;Generate_Random_Turbulence(grid,u_old);Generate_Random_Turbulence(grid,u_new);}

    void Advance_Turbulence()
    {T increment=time_end-time_start;time_start=time_end;time_end+=increment;for(int i=0;i<TV::m;i++) u_old(i).Copy(u_new(i));Generate_Random_Turbulence(grid,u_new);}

    TV Turbulent_Velocity(const TV& X,const T fraction) const
    {TV X_new=wrap(X,grid.domain.min_corner,grid.domain.max_corner),u1,u2;
    for(int i=0;i<TV::m;i++){u1(i)=interpolation.Clamped_To_Array(grid,u_old(i),X_new);u2(i)=interpolation.Clamped_To_Array(grid,u_new(i),X_new);}
    return (1-fraction)*u1+fraction*u2;}

    T Turbulent_Face_Velocity(const int axis,const TV& X,const T fraction) const
    {TV X_new=wrap(X,grid.domain.min_corner,grid.domain.max_corner);
    T w1=interpolation.Clamped_To_Array(grid,u_old(axis),X_new),w2=interpolation.Clamped_To_Array(grid,u_new(axis),X_new);
    return (1-fraction)*w1+fraction*w2;}

//#####################################################################
    void Generate_Random_Turbulence(const GRID<TV>& grid,VECTOR<ARRAY<T,TV_INT>,TV::m>& u) const;
//#####################################################################
};
}
#endif
