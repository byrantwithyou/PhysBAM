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
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Incompressible/Projection/PROJECTION_UNIFORM.h>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;

template<class TV>
class SMOKE_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    enum workaround1{d=TV::m};

public:
    STREAM_TYPE stream_type;
    DEBUG_PARTICLES<TV>& debug_particles;

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

    T cfl;

    GRID<TV> mac_grid;
    MPI_UNIFORM_GRID<TV> *mpi_grid;
    THREAD_QUEUE* thread_queue;    
    PROJECTION_UNIFORM<TV> projection;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T, AVERAGING_UNIFORM<TV, FACE_LOOKUP_UNIFORM<TV> >,LINEAR_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > > advection_scalar;
    BOUNDARY<TV,T> boundary_scalar;
    BOUNDARY<TV,T> *boundary;
    ARRAY<T,TV_INT> density;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;    
    RANGE<TV> source;
    pthread_mutex_t lock;

    SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,int refine=0);
    virtual ~SMOKE_EXAMPLE();

    T CFL(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities);
    void CFL_Threaded(RANGE<TV_INT>& domain,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,T& dt);
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    void Initialize_Grid(TV_INT counts,RANGE<TV> domain)
    {mac_grid.Initialize(counts,domain,true);}
    
    void Initialize_Fields() 
    {for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
    for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=0;}
    
    void Get_Scalar_Field_Sources(const T time)
    {for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next())
        if(source.Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;}

    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Set_Boundary_Conditions(const T time);

//#####################################################################
};
}
#endif
