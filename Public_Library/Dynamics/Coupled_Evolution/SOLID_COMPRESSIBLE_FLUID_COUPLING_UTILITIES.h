//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES
//#####################################################################
#ifndef __SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES__
#define __SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES__
#include <Geometry/Basic_Geometry/POLYGON.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <Incompressible/Collisions_And_Interactions/CUT_CELL.h>
namespace PhysBAM{
template<class TV> class GRID;
template<class T,int d> class VECTOR;
template<class TV> class MPI_UNIFORM_GRID;
template<class TV> class EULER_UNIFORM;
template<class TV> class GRID_BASED_COLLISION_GEOMETRY;
template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;
template<class TV> class EULER_FLUID_FORCES;

template<class TV>
class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef ARRAY<CUT_CELL<T,TV::dimension>*,TV_INT> T_ARRAYS_CUT_CELL;
    typedef ARRAY<int,FACE_INDEX<TV::m> > T_FACE_ARRAYS_INT;
    typedef ARRAY<TV_DIMENSION,FACE_INDEX<TV::m> > T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef LINEAR_INTERPOLATION_UNIFORM<TV,TV_DIMENSION> T_LINEAR_INTERPOLATION_DIMENSION;

public:
    EULER_UNIFORM<TV>& euler;
    MPI_UNIFORM_GRID<TV>* mpi_grid;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_bodies_affecting_fluid;

    ARRAY<bool,TV_INT> uncovered_cells;
    bool thinshell,use_fast_marching,use_higher_order_solid_extrapolation,fluid_affects_solid;
    int number_of_cells_to_extrapolate;
    TV_DIMENSION solid_state;

    T_ARRAYS_DIMENSION_SCALAR U_n;
    ARRAY<bool,FACE_INDEX<TV::m> > solid_fluid_face_time_n;
    ARRAY<bool,TV_INT> cells_inside_fluid_time_n,outside_fluid;
    EULER_FLUID_FORCES<TV>* euler_fluid_forces;
    ARRAY<T,FACE_INDEX<TV::m> > pressure_at_faces;
    ARRAY<T,TV_INT> phi_all_solids_negated;

    ARRAY<bool,TV_INT> near_interface;
    T_ARRAYS_CUT_CELL cut_cells_n,cut_cells_n_p_half,cut_cells_np1;
    ARRAY<T,TV_INT> cell_volumes_n,cell_volumes_n_p_half,cell_volumes_np1;

    ARRAY<bool,TV_INT> psi_n,psi_n_p_half,psi_np1,uncovered_cells_n_p_half;
    ARRAY<TV,TV_INT> advection_velocities_n;

    T_FACE_ARRAYS_DIMENSION_SCALAR accumulated_flux;

    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES(EULER_UNIFORM<TV>& euler_input,MPI_UNIFORM_GRID<TV>* mpi_grid_input=0);
    ~SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES();

//#####################################################################
    void Initialize_Solid_Fluid_Coupling(GRID_BASED_COLLISION_GEOMETRY<TV>* collision_bodies_affecting_fluid_input);
    void Update_Cut_Out_Grid();
    void Fill_Uncovered_Cells();
    void Fill_Solid_Cells(bool fill_pressure_only=false);
    void Project_Fluid_Pressure_At_Neumann_Faces(const ARRAY<T,TV_INT>& p_ghost,ARRAY<T,FACE_INDEX<TV::m> >& p_face) const;
    void Apply_Isobaric_Fix(const T dt,const T time);
    void Extract_Time_N_Data_For_Explicit_Fluid_Forces();
private:
    void Get_Neumann_Data(const TV& location,const T max_distance,TV& normal_direction,T& object_velocity_normal_component,TV& reflected_point) const;
    void Get_Neumann_Data(const TV& location,const T max_distance,TV& object_velocity,TV& reflected_point) const;
    void Extrapolate_State_Into_Solids(ARRAY<T,TV_INT>& phi_all_solids_negated,const int number_of_ghost_cells,const int number_of_cells_to_extrapolate);
    void Compute_Phi_Solids(const int number_of_ghost_cells);
//#####################################################################
public:
    void Snapshot_State(const T_ARRAYS_DIMENSION_SCALAR& U_ghost);
    void Initialize_Collision_Data();
    void Update_Np1_Collision_Data(const T dt);
    void Compute_Intermediate_Solid_Position_Data(const T dt);
    void Revert_Cells_Near_Interface(const int iteration_number);
    void Update_Cells_Near_Interface(const T dt,const int rk_order,const int rk_substep);
    void Compute_Post_Advected_Variables();
//#####################################################################
};
}
#endif
