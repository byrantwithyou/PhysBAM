//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EXAMPLE
//#####################################################################
#ifndef __SOLIDS_EXAMPLE__
#define __SOLIDS_EXAMPLE__

#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class SOLIDS_PARAMETERS;
template<class TV> class SOLIDS_EVOLUTION;
template<class TV> class SOLID_BODY_COLLECTION;

template<class TV>
class SOLIDS_EXAMPLE:public EXAMPLE<TV>,public EXAMPLE_FORCES_AND_VELOCITIES<TV>,public SOLIDS_EVOLUTION_CALLBACKS<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef EXAMPLE<TV> BASE;typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<int,TV::m> T_VECTOR_INT;
    typedef ARRAY<char,TV_INT> T_ARRAYS_CHAR;
    typedef typename MATRIX_POLICY<TV>::TRANSFORMATION_MATRIX T_TRANSFORMATION_MATRIX;
    typedef typename TV::SPIN T_ANGULAR_VELOCITY;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    typedef AVERAGING_UNIFORM<TV> T_AVERAGING;
public:
    using EXAMPLE_FORCES_AND_VELOCITIES<TV>::Set_External_Positions;
    using BASE::output_directory;using BASE::frame_title;using BASE::stream_type;using BASE::initial_time;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::write_last_frame;using BASE::write_time;using BASE::write_substeps_level;using BASE::Set_Write_Substeps_Level;//using BASE::data_directory;
    using BASE::restart;using BASE::Write_Frame_Title;

protected:
    T minimum_collision_thickness; // needed for ray tracing, etc.
public:
    bool use_melting;
    SOLIDS_PARAMETERS<TV>& solids_parameters;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    SOLIDS_EVOLUTION<TV>* solids_evolution; // defaults to newmark
    bool opt_solidssymmqmr,opt_solidscr,opt_solidscg;
    DEBUG_PARTICLES<TV>& debug_particles;
    bool opt_skip_debug_data;

    SOLIDS_EXAMPLE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args);
    SOLIDS_EXAMPLE(const SOLIDS_EXAMPLE&) = delete;
    void operator=(const SOLIDS_EXAMPLE&) = delete;
    virtual ~SOLIDS_EXAMPLE();

    void Set_Minimum_Collision_Thickness(const T minimum_collision_thickness_input=1e-6)
    {minimum_collision_thickness=minimum_collision_thickness_input;}

//#####################################################################
    virtual void Post_Initialization();
    virtual void Preprocess_Frame(const int frame);
    virtual void Postprocess_Frame(const int frame);
    virtual void Preprocess_Substep(const T dt,const T time);
    virtual void Postprocess_Substep(const T dt,const T time);
    void Log_Parameters() const override;
    // solids
    virtual void Initialize_Bodies();
    virtual void Read_Output_Files_Solids(const int frame);
    void After_Initialization() override;
    void Post_Velocity_Advection_Callback(const T dt,const T time){}
    virtual void Write_Output_Files(const int frame) const override;
    void Adjust_Output_Directory_For_MPI(const MPI_SOLIDS<TV>* mpi);
//#####################################################################
};
}
#endif
