//#####################################################################
// Copyright 2011, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_EXAMPLE
//#####################################################################
#ifndef __DEFORMABLE_EXAMPLE__
#define __DEFORMABLE_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
#include "libmain.h"
namespace PhysBAM{

template<class T_GRID> class SOLIDS_FLUIDS_DRIVER_UNIFORM;
template<class T> class DEFORMABLE_EXAMPLE;
template<class TV> class DEFORMABLE_GRAVITY;
template<class TV,int d> class FINITE_VOLUME;

struct BASE_WRAPPER
{
    int id;
    DEFORMABLE_EXAMPLE<float>& de;

    BASE_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input,int id_input);
    virtual ~BASE_WRAPPER();
};

struct OBJECT_WRAPPER: public BASE_WRAPPER
{
    OBJECT_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input,int id_input);
};

struct DEFORMABLE_BODY_WRAPPER: public OBJECT_WRAPPER
{
    int structure_index;
    int free_particles_index;
    int enclosing_structure_index;
    ARRAY<int> particle_map;

    static int fixed_id(int s = -1){static int i = s; return i;}
    DEFORMABLE_BODY_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input);
};

struct SCRIPTED_GEOMETRY_WRAPPER: public OBJECT_WRAPPER
{
    int rigid_index;
    static int fixed_id(int s = -1){static int i = s; return i;}
    SCRIPTED_GEOMETRY_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input);
};

struct FORCE_WRAPPER: public BASE_WRAPPER
{
    FORCE_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input,int id_input=0);
    virtual ~FORCE_WRAPPER();
};

struct GRAVITY_WRAPPER: public FORCE_WRAPPER
{
    typedef VECTOR<float,3> TV;
    float magnitude;
    VECTOR<float,3> direction;

    ARRAY<DEFORMABLE_GRAVITY<TV>*> force_instances;

    static int fixed_id(int s = -1){static int i = s; return i;}
    GRAVITY_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input);
    virtual ~GRAVITY_WRAPPER();
};

struct VOLUMETRIC_FORCE_WRAPPER: public FORCE_WRAPPER
{
    typedef VECTOR<float,3> TV;
    float stiffness;
    float poissons_ratio;
    float damping;

    ARRAY<FINITE_VOLUME<TV,3>*> force_instances;

    static int fixed_id(int s = -1){static int i = s; return i;}
    VOLUMETRIC_FORCE_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input);
    virtual ~VOLUMETRIC_FORCE_WRAPPER();

    void Propagate_Parameters();
};

void register_accessors(int last_id);
inline void Register_Wrapper_Ids()
{
    int next_id=0;
    DEFORMABLE_BODY_WRAPPER::fixed_id(++next_id);
    SCRIPTED_GEOMETRY_WRAPPER::fixed_id(++next_id);
    GRAVITY_WRAPPER::fixed_id(++next_id);
    VOLUMETRIC_FORCE_WRAPPER::fixed_id(++next_id);
    register_accessors(next_id);
}

template<class T_input>
class DEFORMABLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    bool fully_implicit;
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> >& driver;

    ARRAY<OBJECT_WRAPPER*> object_wrappers;
    ARRAY<FORCE_WRAPPER*> force_wrappers;
    ARRAY<PAIR<OBJECT_WRAPPER*,FORCE_WRAPPER*> > new_forces_relations;
    bool added_body;
    bool want_log;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;using BASE::parse_args;using BASE::stream_type;

    DEFORMABLE_EXAMPLE(const STREAM_TYPE stream_type);
    ~DEFORMABLE_EXAMPLE();

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_End_Callback(const T dt,const T time) {}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return true;}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) {}

    void Initialize_Bodies() PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
    OBJECT_WRAPPER* Add_Simulation_Object(const data_exchange::simulation_object& body);
    DEFORMABLE_BODY_WRAPPER* Add_Deformable_Body(const data_exchange::deformable_body& body);
    SCRIPTED_GEOMETRY_WRAPPER* Add_Rigid_Body(const data_exchange::ground_plane& body);
    SCRIPTED_GEOMETRY_WRAPPER* Add_Rigid_Body(const data_exchange::scripted_geometry& body);
    FORCE_WRAPPER* Add_Force(const data_exchange::force& f);
    bool Instantiate_Force(FORCE_WRAPPER& wrapper,DEFORMABLE_BODY_WRAPPER& body);
    bool Instantiate_Force(GRAVITY_WRAPPER& wrapper,DEFORMABLE_BODY_WRAPPER& body);
    bool Instantiate_Force(VOLUMETRIC_FORCE_WRAPPER& wrapper,DEFORMABLE_BODY_WRAPPER& body);
    void Initialize_Simulation();
    void Initialize_After_New_Bodies();
    void Simulate_Frame();
    void Add_New_Forces();
//#####################################################################
};
}
#endif

