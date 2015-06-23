//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_REFLECTION_ATTENUATION
//#####################################################################
// Attenuates to far field value for open walls (constant_extrapolation=true) and reflects for solid walls
//#####################################################################
#ifndef __BOUNDARY_REFLECTION_ATTENUATION__
#define __BOUNDARY_REFLECTION_ATTENUATION__

#include <Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_REFLECTION_ATTENUATION:public BOUNDARY_REFLECTION_UNIFORM<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,TV::m> TV_SIDES;typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
public:
    typedef BOUNDARY_REFLECTION_UNIFORM<TV,T2> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;
    using BASE::fixed_boundary_value;using BASE::Find_Ghost_Regions;using BASE::Boundary;

    T linear_attenuation;

    BOUNDARY_REFLECTION_ATTENUATION(const TV_SIDES& constant_extrapolation,const T2 far_field_value,const T linear_attenuation_input)
    {
        linear_attenuation=linear_attenuation_input;
        Set_Constant_Extrapolation(constant_extrapolation);
        fixed_boundary_value=far_field_value;
    }

//#####################################################################
    virtual void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    virtual void Fill_Single_Ghost_Region(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells=3) const override ;
    virtual void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const override {} // do nothing
private:
    T2 Attenuate_To_Far_Field_Value(const T2 boundary_value,const T dt) const;
//#####################################################################
};
}
#endif
