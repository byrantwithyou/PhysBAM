//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WIND_DRAG_3D
//#####################################################################
#ifndef __WIND_DRAG_3D__
#define __WIND_DRAG_3D__

#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Solids/Forces_And_Torques/SOLIDS_FORCES.h>
namespace PhysBAM{
template<class TV> class RIGID_BODY;

template<class T_input>
class WIND_DRAG_3D:public SOLIDS_FORCES<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename BASE::RIGID_FREQUENCY_DATA RIGID_FREQUENCY_DATA;
    typedef typename BASE::DEFORMABLE_FREQUENCY_DATA DEFORMABLE_FREQUENCY_DATA;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SIMPLEX;

    TRIANGULATED_SURFACE<T>* triangulated_surface;
    RIGID_BODY<TV>* rigid_body;
    bool use_constant_wind;
    T constant_wind_viscosity;
    TV constant_wind;
    bool use_spatially_varying_wind;
    T spatially_varying_wind_viscosity;
    GRID<TV> V_grid;
    ARRAY<TV,VECTOR<int,3> >* spatially_varying_wind;
    T wind_density;
    ARRAY<T,VECTOR<int,3> > *spatially_varying_wind_density,*spatially_varying_wind_pressure;
    T linear_normal_viscosity; // uses vertex normals
private:
    static const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<TV,T> interpolation;
    static const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<TV,TV> vector_interpolation;
    struct OPTIMIZATION{
        OPTIMIZATION()
            :one_third_area((T)0)
        {}

        T one_third_area;TV inward_normal,center,wind_velocity;
    };
    ARRAY<OPTIMIZATION> optimization;
    ARRAY<int> force_elements;
    ARRAY<int> force_particles;
    ARRAY<TV> vertex_normals; // TODO: use triangulated_surface.vertex_normals
    
    MPI_SOLIDS<TV>* mpi_solids;
public:

    template<class T_OBJECT> WIND_DRAG_3D(T_OBJECT& object,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    WIND_DRAG_3D(RIGID_BODY<TV>& rigid_body_input,DEFORMABLE_PARTICLES<TV>& deformable_body_particles_input);
    virtual ~WIND_DRAG_3D();

    void Use_Constant_Wind(const T viscosity_input,const TV& wind_input=TV())
    {use_constant_wind=true;constant_wind_viscosity=viscosity_input;constant_wind=wind_input;}

    void Use_Spatially_Varying_Wind(const T viscosity_input,const GRID<TV>& grid_input,ARRAY<TV,VECTOR<int,3> >& V_input)
    {use_spatially_varying_wind=true;spatially_varying_wind_viscosity=viscosity_input;
    V_grid=grid_input;spatially_varying_wind=&V_input;}

    void Set_Wind_Density(const T wind_density_input)
    {wind_density=wind_density_input;}

    void Set_Wind_Density(ARRAY<T,VECTOR<int,3> >& density_input)
    {spatially_varying_wind_density=&density_input;}

    void Set_Wind_Pressure(ARRAY<T,VECTOR<int,3> >& pressure_input) // only valid for volumetric objects
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

    void Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose=false) const override
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) override
    {}

    T CFL_Strain_Rate() const override
    {return FLT_MAX;}

    void Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency,ARRAY_VIEW<RIGID_FREQUENCY_DATA> rigid_frequency) override
    {}

    VECTOR<T,3> Add_Velocity_Independent_Forces_Helper(TV relative_velocity,int t) const;

public:

//#####################################################################
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Position_Based_State(const T time) override;
    void Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const override;
    void Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time) const override;
//#####################################################################
};
}
#endif
