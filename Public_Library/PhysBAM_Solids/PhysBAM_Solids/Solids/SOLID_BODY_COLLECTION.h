//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_BODY_COLLECTION
//#####################################################################
#ifndef __SOLID_BODY_COLLECTION__
#define __SOLID_BODY_COLLECTION__

#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_FORCE_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/SOLIDS_FORCES.h>
namespace PhysBAM{

template<class TV> class EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class DEFORMALBLE_OBJECT_COLLISIONS;
template<class TV> class SOLIDS_PARAMETERS;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class RIGID_BODY_CLUSTER_BINDINGS;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class RIGIDS_FORCES;
template<class TV> class SOLID_FORCE_COLLECTION;

template<class TV>
class SOLID_BODY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename RIGIDS_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_RIGID;
    typedef typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_DEFORMABLE;
public:
    COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    SOLID_FORCE_COLLECTION<TV>& solid_force_collection;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>* example_forces_and_velocities;
public:
    bool print_diagnostics;
    bool print_residuals;
    bool simulate;
    int iterations_used_diagnostic;

    SOLID_BODY_COLLECTION(EXAMPLE_FORCES_AND_VELOCITIES<TV>* example_forces_and_velocities_input);
    virtual ~SOLID_BODY_COLLECTION();

    void Print_Diagnostics(const bool print_diagnostics_input=true)
    {print_diagnostics=print_diagnostics_input;}

    void Print_Residuals(const bool print_residuals_input=true)
    {print_residuals=print_residuals_input;}

//#####################################################################
    void Update_Time_Varying_Material_Properties(const T time);
    void Update_Simulated_Particles();
    void Compute_Linear_Momentum(TV& linear_momentum) const;
    TV Compute_Momentum() const;
    void Adjust_Mesh_For_Self_Collision(const T time);
    void Read(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const int static_frame,const bool include_static_variables,const bool read_rigid_body,
        const bool read_deformable_body,const bool read_from_every_process,ARRAY<int>* needs_init=0,ARRAY<int>* needs_destroy=0);
    void Write(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const int first_frame,const bool include_static_variables,const bool write_rigid_body,
        const bool write_deformable_body,const bool write_from_every_process,const bool output_interaction_pairs) const;
//#####################################################################
};
}
#endif
