//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_COLLECTION
//#####################################################################
#ifndef __RIGID_BODY_COLLECTION__
#define __RIGID_BODY_COLLECTION__

#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class RIGID_BODY_CLUSTER_BINDINGS;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class RIGIDS_FORCES;
template<class TV> class RIGID_BODY_EVOLUTION_PARAMETERS;
template<class TV> struct ALLOCATE_HELPER{virtual RIGID_BODY<TV>* Create(int index=0)=0;virtual ~ALLOCATE_HELPER(){}};
template<class TV,class ID> class STRUCTURE_LIST;
template<class TV> class STRUCTURE;

template<class TV>
class RIGID_BODY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_TYPED_READ_WRITE;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles;
    COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list;
    STRUCTURE_LIST<TV,int>& structure_list;
    bool always_create_structure;
    HASHTABLE<std::string,int>& structure_hash; // maps to id
    mutable bool is_stale_key,is_stale_active;
    mutable ARRAY<int> *frame_list_key,*frame_list_active;
    mutable ARRAY<std::string,int> rigid_body_names;
    bool check_stale;
    int last_read_key,last_read_active;
    ALLOCATE_HELPER<TV>* allocate_helper;
    bool owns_collision_body_list;

    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings;
    RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>* rigids_example_forces_and_velocities;
    ARRAY<int> simulated_rigid_body_particles;
    ARRAY<int> dynamic_rigid_body_particles;
    ARRAY<int> static_rigid_bodies,kinematic_rigid_bodies,static_and_kinematic_rigid_bodies;
    ARRAY<RIGIDS_FORCES<TV>*> rigids_forces;
    
    bool print_diagnostics;
    bool print_residuals;
    bool print_energy;
    int iterations_used_diagnostic;

    RIGID_BODY<TV>* New_Body(int index);

    RIGID_BODY_COLLECTION(COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input);
    virtual ~RIGID_BODY_COLLECTION();

    RIGID_BODY_STATE<TV> State(const int particle) const
    {return RIGID_BODY_STATE<TV>(rigid_body_particles.frame(particle),rigid_body_particles.twist(particle));}

    void Set_State(const int particle,const RIGID_BODY_STATE<TV>& state)
    {rigid_body_particles.frame(particle)=state.frame;rigid_body_particles.twist(particle)=state.twist;}

    bool Exists(const int particle) const
    {return particle>=0 && particle<rigid_body_particles.Size() && rigid_body_particles.rigid_body(particle);}

    bool Is_Active(const int particle) const
    {return Exists(particle) && Rigid_Body(particle).particle_index>=0;}

    void Deactivate_Body(const int p)
    {assert(Exists(p) && Rigid_Body(p).particle_index>=0);Rigid_Body(p).particle_index=~Rigid_Body(p).particle_index;}

    void Reactivate_Body(const int p)
    {assert(Exists(p) && Rigid_Body(p).particle_index<=0);Rigid_Body(p).particle_index=~Rigid_Body(p).particle_index;}

    RIGID_BODY<TV>& Rigid_Body(const int particle_index);
    const RIGID_BODY<TV>& Rigid_Body(const int particle_index) const;
    int Add_Rigid_Body(RIGID_BODY<TV>* rigid_body,const int simplicial_boundary_id,const int implicit_object_id,const int simplicial_interior_id);
    int Add_Rigid_Body_And_Geometry(RIGID_BODY<TV>* rigid_body);
    int Add_Rigid_Body(const STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor=1,const bool read_simplicial_boundary=true,const bool read_implicit_object=true,
        const bool read_simplicial_interior=false,const bool read_rgd_file=true);
    int Add_Rigid_Body(const STREAM_TYPE stream_type,const bool thin_shell,const std::string& basename,const T scaling_factor=1,const bool read_simplicial_boundary=true,
        const bool read_implicit_object=true,const bool read_simplicial_interior=false,const bool read_rgd_file=true);
    int Add_Rigid_Body(RIGID_BODY<TV>* rigid_body,STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,
        const bool read_simplicial_boundary,const bool read_implicit_object,const bool read_simplicial_interior,const bool read_rgd_file);
    void Update_Angular_Velocity();
    void Update_Angular_Momentum();
    void Update_Angular_Velocity(const ARRAY<int>& particle_indices);
    void Update_Angular_Momentum(const ARRAY<int>& particle_indices);
    void Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame,ARRAY<int>* needs_init=0,ARRAY<int>* needs_destroy=0);
    void Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const; // TODO: optionally skip certain kinds of structures in output
    void Update_Simulated_Particles();

    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T scale,const T time) const;

    void Update_Position_Based_State(const T time);
    void Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const;
    void Print_Energy(const T time,const int step) const;
    T CFL_Rigid(const RIGID_BODY_EVOLUTION_PARAMETERS<TV>& rigid_body_evolution_parameters,const bool verbose_dt);
    int Add_Force(RIGIDS_FORCES<TV>* force);

    void Update_Kinematic_Particles();
    bool Register_Analytic_Replacement_Structure(const std::string& filename,const T scaling_factor,STRUCTURE<TV>* structure); // passing in zero skips reading
    bool Find_Or_Read_Structure(const STREAM_TYPE stream_type,ARRAY<int>& structure_id,const std::string& filename,const T scaling_factor,const TV& center);
    void Destroy_Unreferenced_Geometry();
};
}
#endif
