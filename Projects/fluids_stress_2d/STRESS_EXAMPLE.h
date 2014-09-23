//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STRESS_EXAMPLE__
#define __STRESS_EXAMPLE__
#include <Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;

template<class TV>
class STRESS_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    enum BC_TYPE {DIRICHLET=0,NEUMANN=1};
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
    bool use_advection;
    bool use_reduced_advection;
    BC_TYPE bc_type;
    T inv_Wi;
    bool use_du_terms;

    GRID<TV> grid;
    ARRAY<T,TV_INT> bc_phi;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities,prev_face_velocities;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T> boundary;
    BOUNDARY_MAC_GRID_PERIODIC<TV,int> boundary_int;
    BOUNDARY_MAC_GRID_PERIODIC<TV,SYMMETRIC_MATRIX<T,TV::m> > boundary_stress;
    LEVELSET<TV> levelset;
    DEBUG_PARTICLES<TV>& debug_particles;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> polymer_stress,prev_polymer_stress;

    STRESS_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~STRESS_EXAMPLE();
    
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    virtual void Begin_Time_Step(const T time)=0;
    virtual void End_Time_Step(const T time)=0;
    virtual SYMMETRIC_MATRIX<T,TV::m> Polymer_Stress(const TV& X,T time)=0;
    virtual SYMMETRIC_MATRIX<T,TV::m> Polymer_Stress_Forcing_Term(const TV& X,T time)=0;
    virtual void Get_Velocities(T time)=0;
    virtual void Get_Initial_Polymer_Stresses()=0;
//#####################################################################
};
}
#endif
