//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_ETHER_DRAG
//#####################################################################
#ifndef __RIGID_ETHER_DRAG__
#define __RIGID_ETHER_DRAG__

#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Rigids/Forces_And_Torques/RIGID_POINTWISE_FORCE.h>
namespace PhysBAM{

template<class TV>
class RIGID_ETHER_DRAG:public RIGID_POINTWISE_FORCE<TV>
{
private:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef RIGID_POINTWISE_FORCE<TV> BASE;
    using BASE::rigid_body_collection;using BASE::force_rigid_body_particles;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;

    bool use_constant_wind;
    T constant_wind_viscosity,constant_wind_angular_viscosity;
    TV constant_wind;
    bool use_spatially_varying_wind;
    T spatially_varying_wind_viscosity;
    RANGE<TV> spatially_varying_wind_domain;
    GRID<TV> V_grid;
    ARRAY<TV,TV_INT>* spatially_varying_wind;
    LINEAR_INTERPOLATION_UNIFORM<TV,TV> vector_interpolation;

    RIGID_ETHER_DRAG(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_rigid_body_particles_input,T dynamic_ether_viscosity=0,T angular_viscosity=0);
    RIGID_ETHER_DRAG(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_rigid_body_particles_input,T dynamic_ether_viscosity=0,T angular_viscosity=0);
    virtual ~RIGID_ETHER_DRAG();

    void Use_No_Drag()
    {use_constant_wind=use_spatially_varying_wind=false;}

    void Use_Constant_Wind(const T viscosity_input,const TV& wind_input=TV())
    {if(!viscosity_input && wind_input==TV()) Use_No_Drag();
    else{use_constant_wind=true;constant_wind_viscosity=viscosity_input;constant_wind=wind_input;}}

    void Use_Spatially_Varying_Wind(const T viscosity_input,const RANGE<TV>& domain_input,const GRID<TV>& grid_input,ARRAY<TV,TV_INT>& V_input)
    {use_spatially_varying_wind=true;spatially_varying_wind_viscosity=viscosity_input;
    spatially_varying_wind_domain=domain_input;V_grid=grid_input;spatially_varying_wind=&V_input;}

    TV Spatially_Varying_Wind_Velocity(const TV& X) const
    {return vector_interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind,X);}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
