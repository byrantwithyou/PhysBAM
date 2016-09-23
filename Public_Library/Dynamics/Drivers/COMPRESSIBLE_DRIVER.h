//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPRESSIBLE_DRIVER__
#define __COMPRESSIBLE_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Tools/Ordinary_Differential_Equations/DRIVER.h>
namespace PhysBAM{

template<class TV> class COMPRESSIBLE_EXAMPLE;

template<class TV>
class COMPRESSIBLE_DRIVER:public DRIVER<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef DRIVER<TV> BASE;

    using BASE::time;using BASE::Write_Substep;using BASE::Write_Output_Files;
protected:
    COMPRESSIBLE_EXAMPLE<TV>& example;
public:
    COMPRESSIBLE_DRIVER(COMPRESSIBLE_EXAMPLE<TV>& example);
    virtual ~COMPRESSIBLE_DRIVER();
    
    void Initialize();
    void Advance_To_Target_Time(const T target_time);
private:
    void Setup_Fluids(const T time);
    void Restore_Solids_To_Time_N(const T time_n_plus_one);
    void Restore_Solids_To_Time_N_Plus_One();
    void Advect_Fluid(const T dt,const int substep);
    void Advance_Fluid_One_Time_Step_Implicit_Part(const T dt_projection,const T time_projection,const int substep);
    T Compute_Dt(const T time,const T target_time,bool& done);

//#####################################################################
};
}
#endif
