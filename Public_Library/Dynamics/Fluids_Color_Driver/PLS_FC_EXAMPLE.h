//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_FC_EXAMPLE__
#define __PLS_FC_EXAMPLE__
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Finite_Elements/LEVELSET_COLOR.h>
#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class LEVELSET_MULTIPLE;
template<class TV> class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM;
template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;
template<class TV> class DEBUG_PARTICLES;

template<class TV_input>
class PLS_FC_EXAMPLE:public LEVELSET_CALLBACKS<TV_input>
{
    typedef TV_input TV;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    enum WORKAROUND1{num_bc=3};
    STREAM_TYPE stream_type;
    T initial_time;
    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    int substeps_delay_frame;
    bool write_output_files;
    std::string output_directory;
    int restart;
    int number_of_ghost_cells;

    T dt,time;
    int time_steps_per_frame;
    bool use_preconditioner;
    int max_iter;
    T solver_tolerance;
    bool dump_matrix;
    bool sparse_dump_matrix;
    bool use_advection;
    bool use_reduced_advection;
    bool omit_solve;
    ARRAY<T> mu,rho;
    int number_of_colors;
    bool use_discontinuous_velocity;
    bool use_p_null_mode;
    bool use_level_set_method;
    bool use_pls;
    bool dump_largest_eigenvector;
    bool save_pressure;
    bool test_system;
    bool use_polymer_stress;
    ARRAY<T> polymer_stress_coefficient,inv_Wi;

    int num_multigrid_levels;
    bool use_multigrid;

    GRID<TV> grid;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>& particle_levelset_evolution_multiple;
    VECTOR<ARRAY<T,TV_INT>,num_bc> bc_phis; // 0=Neumann, 1=Dirichlet, 2=Slip
    ARRAY<int,TV_INT> cell_color,prev_cell_color;
    ARRAY<int,FACE_INDEX<TV::dimension> > face_color,prev_face_color;
    ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > > face_velocities,prev_face_velocities;
    ARRAY<int,TV_INT> pressure_color;
    ARRAY<T,TV_INT> pressure;
    ADVECTION<TV,T>& advection_scalar;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T> boundary;
    BOUNDARY_MAC_GRID_PERIODIC<TV,int> boundary_int;
    BOUNDARY_MAC_GRID_PERIODIC<TV,SYMMETRIC_MATRIX<T,TV::m> > boundary_symmetric;
    LEVELSET_COLOR<TV> levelset_color;
    DEBUG_PARTICLES<TV>& debug_particles;
    ARRAY<ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> > polymer_stress,prev_polymer_stress;

    PLS_FC_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~PLS_FC_EXAMPLE();
    
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE;
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<TV>& levelset_multiple,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE;
    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Merge_Velocities(ARRAY<T,FACE_INDEX<TV::dimension> >& V,const ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > > u,const ARRAY<int,FACE_INDEX<TV::dimension> >& color) const;
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    virtual void Begin_Time_Step(const T time)=0;
    virtual void End_Time_Step(const T time)=0;
    virtual SYMMETRIC_MATRIX<T,TV::m> Polymer_Stress(const TV& X,int color,T time)=0;
    virtual SYMMETRIC_MATRIX<T,TV::m> Polymer_Stress_Forcing_Term(const TV& X,int color,T time)=0;
    virtual TV Jump_Interface_Condition(const TV& X,int color0,int color1,T time)=0;
    virtual TV Volume_Force(const TV& X,int color,T time)=0;
    virtual TV Velocity_Jump(const TV& X,int color0,int color1,T time)=0;
    virtual void Get_Initial_Velocities(T time)=0;
    virtual void Get_Initial_Polymer_Stresses()=0;
    int Color_At_Cell(const TV_INT& index) const;
    int Color_At_Cell(const TV_INT& index,T& phi) const;
    void Rebuild_Levelset_Color();
    void Make_Levelsets_Consistent();
    void Fill_Levelsets_From_Levelset_Color();
    void Enforce_Phi_Boundary_Conditions();
//#####################################################################
};
}
#endif
