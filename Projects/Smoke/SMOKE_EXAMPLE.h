//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_EXAMPLE__
#define __SMOKE_EXAMPLE__
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Interpolation/CUBIC_MONOTONIC_INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>

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
    T beta;
    T cfl;
    GRID<TV> grid;
    MPI_UNIFORM_GRID<TV> *mpi_grid;
    PROJECTION_UNIFORM<TV> projection;
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T, AVERAGING_UNIFORM<TV, FACE_LOOKUP_UNIFORM<TV> >,T_INTERPOLATION > advection_scalar;
    BOUNDARY<TV,T> boundary_scalar;
    BOUNDARY<TV,T> *boundary;
    ARRAY<T,TV_INT> density;
    ARRAY<T,TV_INT> temperature;  //add temperature
    VECTOR<VECTOR<bool,2>,TV::m> domain_boundary;    
    RANGE<TV> source1;
    RANGE<TV> source2;
    
    // EAPIC
    bool use_eapic;
    int eapic_order;
    bool nrs;
    int np;
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> weights; // face weights
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> weights0; // face weights of X0
    SMOKE_PARTICLES<TV>& particles;
    ARRAY<T,TV_INT> mass;
    ARRAY<T,FACE_INDEX<TV::m> > face_mass;

    SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,int refine=0);
    virtual ~SMOKE_EXAMPLE();
    T CFL(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void CFL_Threaded(RANGE<TV_INT>& domain,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T& dt);
    T Time_At_Frame(const int frame) const;
    void Initialize_Grid(TV_INT counts,RANGE<TV> domain);
    void Initialize_Fields();
    void Get_Scalar_Field_Sources(const T time);
    void Set_Weights(int order);
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Set_Boundary_Conditions(const T time, const T source_velocities = 0);

//#####################################################################
};
}
#endif
