//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_REFLECTION_WATER
//#####################################################################
#ifndef __BOUNDARY_REFLECTION_WATER__
#define __BOUNDARY_REFLECTION_WATER__

#include <Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_REFLECTION_WATER:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef BOUNDARY<TV,T2> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;

    bool use_extrapolation_mode;
    T tolerance;
private:
    ARRAY<T,TV_INT>* phi;
    ARRAY<T,FACE_INDEX<TV::m> >* V;
public:

    BOUNDARY_REFLECTION_WATER(const bool left_constant_extrapolation_input=false,const bool right_constant_extrapolation_input=false,const bool bottom_constant_extrapolation_input=false,
        const bool top_constant_extrapolation_input=false,const bool front_constant_extrapolation_input=false,const bool back_constant_extrapolation_input=false)
        :phi(0),V(0)
    {
        Set_Constant_Extrapolation(left_constant_extrapolation_input,right_constant_extrapolation_input,bottom_constant_extrapolation_input,top_constant_extrapolation_input,
            front_constant_extrapolation_input,back_constant_extrapolation_input);
        Use_Extrapolation_Mode(false);
        Set_Tolerance();
    }

    void Use_Extrapolation_Mode(const bool use=true)
    {use_extrapolation_mode=use;}

    void Set_Phi_And_Velocity_Pointers(ARRAY<T,TV_INT>& phi_input,ARRAY<T,FACE_INDEX<TV::m> >& V_input)
    {phi=&phi_input;V=&V_input;}

    void Set_Tolerance(const T tolerance_input=(T)9.8/24)  // dt*gravity where dt=1/24 is based on the length of a frame
    {tolerance=tolerance_input;}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<T2,TV_INT>& u,ARRAY<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
//#####################################################################
};
}
#endif
