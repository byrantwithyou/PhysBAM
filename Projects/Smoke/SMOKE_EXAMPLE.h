//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_EXAMPLE__
#define __SMOKE_EXAMPLE__
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Interpolation/CUBIC_MONOTONIC_INTERPOLATION_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Incompressible/Projection/PROJECTION_UNIFORM.h>

namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class SMOKE_PARTICLES;


template<class TV>
class SMOKE_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    enum workaround1{d=TV::m};
    // typedef CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > T_INTERPOLATION;
    typedef LINEAR_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > T_INTERPOLATION;

public:
    STREAM_TYPE stream_type;
    DEBUG_PARTICLES<TV>& debug_particles;
    int ghost;

    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    int restart;
    std::string frame_title;
    int write_substeps_level;
    bool write_debug_data;
    std::string output_directory;
    bool N_boundary;
    bool debug_divergence;
    T alpha;
    T cfl;
    GRID<TV> mac_grid;
    MPI_UNIFORM_GRID<TV> *mpi_grid;
    THREAD_QUEUE* thread_queue;    
    PROJECTION_UNIFORM<TV> projection;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T, AVERAGING_UNIFORM<TV, FACE_LOOKUP_UNIFORM<TV> >,T_INTERPOLATION > advection_scalar;
    BOUNDARY<TV,T> boundary_scalar;
    BOUNDARY<TV,T> *boundary;
    ARRAY<T,TV_INT> density;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;    
    RANGE<TV> source;
    pthread_mutex_t lock;

    // EAPIC
    bool use_eapic;
    int eapic_order;
    PARTICLE_GRID_WEIGHTS<TV>* weights; // cell center weights
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> face_weights; // face weights
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> face_weights0; // face weights of X0
    SMOKE_PARTICLES<TV>& particles;
    ARRAY<T,TV_INT> mass;
    ARRAY<T,FACE_INDEX<TV::m> > face_mass;

    SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,int refine=0);
    virtual ~SMOKE_EXAMPLE();
    T CFL(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities);
    void CFL_Threaded(RANGE<TV_INT>& domain,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,T& dt);
    T Time_At_Frame(const int frame) const;
    void Initialize_Grid(TV_INT counts,RANGE<TV> domain);
    void Initialize_Fields();
    void Get_Scalar_Field_Sources(const T time);
    void Set_Weights(PARTICLE_GRID_WEIGHTS<TV>* weights_input);
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Set_Boundary_Conditions(const T time);

//#####################################################################
};
}
#endif
