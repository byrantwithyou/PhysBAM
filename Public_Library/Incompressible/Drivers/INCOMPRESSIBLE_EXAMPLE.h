//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INCOMPRESSIBLE_EXAMPLE__
#define __INCOMPRESSIBLE_EXAMPLE__
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Incompressible/Projection/PROJECTION_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class INCOMPRESSIBLE_EXAMPLE:public RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    enum workaround1{d=TV::m};

public:
    STREAM_TYPE stream_type;
    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    int restart;
    std::string frame_title;
    int write_substeps_level;
    bool write_debug_data;
    bool analytic_test;
    int order;
    std::string output_directory;

    int number_of_ghost_cells;
    T cfl;
    bool use_viscosity;

    GRID<TV> mac_grid;
    MPI_UNIFORM_GRID<TV> *mpi_grid;
    PROJECTION_DYNAMICS_UNIFORM<TV> projection;
    INCOMPRESSIBLE_UNIFORM<TV> incompressible;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities,face_velocities_save;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T> advection_scalar;
    BOUNDARY<TV,T> boundary_scalar;
    BOUNDARY<TV,T> *boundary;
    ARRAY<T,TV_INT> density,temperature;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;    
    RIGID_BODY_COLLECTION<TV> rigid_body_collection;

    INCOMPRESSIBLE_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~INCOMPRESSIBLE_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize_Fields()=0;
    virtual void Initialize_Confinement() {}
    virtual void Get_Scalar_Field_Sources(const T time)=0;
    virtual void Set_Boundary_Conditions(const T time)=0;
    virtual void Initialize_Bodies() {}
    virtual void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const ARRAY<T,TV_INT>& density_ghost,const T dt,const T time){}

//#####################################################################
};
}
#endif
