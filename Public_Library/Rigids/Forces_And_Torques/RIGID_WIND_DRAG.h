//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_WIND_DRAG
//#####################################################################
#ifndef __RIGID_WIND_DRAG__
#define __RIGID_WIND_DRAG__

#include <Tools/Grids_Uniform/GRID.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class GRID;

template<class TV>
class RIGID_WIND_DRAG:public RIGIDS_FORCES<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::SCALAR T;
public:
    typedef RIGIDS_FORCES<TV> BASE;
    typedef typename BASE::FREQUENCY_DATA RIGID_FREQUENCY_DATA;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension-1>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_SIMPLEX;

    RIGID_BODY<TV>* rigid_body;
    bool use_constant_wind;
    T constant_wind_viscosity;
    TV constant_wind;
    bool use_spatially_varying_wind;
    T spatially_varying_wind_viscosity;
    GRID<TV> V_grid;
    ARRAY<TV,TV_INT> *spatially_varying_wind;
    T wind_density;
    ARRAY<T,TV_INT> *spatially_varying_wind_density,*spatially_varying_wind_pressure;
    T linear_normal_viscosity; // uses vertex normals
    mutable T surface_area;
private:
    static const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<TV,T> interpolation;
    static const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<TV,TV> vector_interpolation;
    struct OPTIMIZATION{
        OPTIMIZATION()
            :area_over_m((T)0)
        {}

        T area_over_m;TV inward_normal,center,wind_velocity;
    };
    ARRAY<OPTIMIZATION> optimization;
    ARRAY<TV> vertex_normals; // TODO: use triangulated_surface.vertex_normals
public:

    RIGID_WIND_DRAG(RIGID_BODY<TV>& rigid_body_input);
    virtual ~RIGID_WIND_DRAG();

    void Use_Constant_Wind(const T viscosity_input,const TV& wind_input=TV())
    {use_constant_wind=true;constant_wind_viscosity=viscosity_input;constant_wind=wind_input;}

    void Use_Spatially_Varying_Wind(const T viscosity_input,const GRID<TV>& grid_input,ARRAY<TV,TV_INT>& V_input)
    {use_spatially_varying_wind=true;spatially_varying_wind_viscosity=viscosity_input;
    V_grid=grid_input;spatially_varying_wind=&V_input;}

    void Set_Wind_Density(const T wind_density_input)
    {wind_density=wind_density_input;}

    void Set_Wind_Density(ARRAY<T,TV_INT>& density_input)
    {spatially_varying_wind_density=&density_input;}

    void Set_Wind_Pressure(ARRAY<T,TV_INT>& pressure_input) // only valid for volumetric objects
    {spatially_varying_wind_pressure=&pressure_input;}

    void Use_Linear_Normal_Viscosity(const T viscosity_input)
    {linear_normal_viscosity=viscosity_input;}

private:
    TV Spatially_Varying_Wind_Velocity(const TV& X) const
    {return vector_interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind,X);}

    T Spatially_Varying_Wind_Density(const TV& X) const
    {return interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind_density,X);}

    T Spatially_Varying_Wind_Pressure(const TV& X) const
    {return interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind_pressure,X);}

    void Enforce_Definiteness(const bool enforce_definiteness_input) override
    {}

    T CFL_Strain_Rate() const override
    {return FLT_MAX;}

    void Initialize_CFL(ARRAY_VIEW<RIGID_FREQUENCY_DATA> rigid_frequency) override
    {}

    void Update_Mpi(const ARRAY<bool>& particle_is_simulated) override
    {}

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override
    {}

    TV Add_Velocity_Independent_Forces_Helper(TV relative_velocity,int t) const;

public:

//#####################################################################
    void Update_Position_Based_State(const T time) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
//#####################################################################
};
}
#endif
