//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonthan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_BODY_COLLECTION
//#####################################################################
#ifndef __DEFORMABLE_BODY_COLLECTION__
#define __DEFORMABLE_BODY_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class BINDING_LIST;
template<class TV> class DEFORMABLE_GEOMETRY_COLLECTION;
template<class TV> class COLLISION_GEOMETRY_COLLECTION;
template<class TV> class DEFORMABLE_OBJECT_COLLISIONS;
template<class TV> class DEFORMABLE_PARTICLES;
template<class TV> class SOFT_BINDINGS;
template<class TV> class TRIANGLE_REPULSIONS;
template<class TV> class TRIANGLE_COLLISIONS;
template<class TV> class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV> class COLLISION_PENALTY_FORCES;
template<class TV> class DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class MPI_SOLIDS;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
template<class TV> class STRUCTURE;
template<class TV> class DEFORMABLE_FORCE_COLLECTION;

template<class TV>
class DEFORMABLE_BODY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_DEFORMABLE;
public:
    DEFORMABLE_PARTICLES<TV>& particles;
    ARRAY<STRUCTURE<TV>*> structures;
    DEFORMABLE_FORCE_COLLECTION<TV>& deformable_force_collection;

    int Add_Structure(STRUCTURE<TV>* structure)
    {return structures.Append(structure);}

    template<class T_STRUCTURE> T_STRUCTURE
    Find_Structure(const int index=0)
    {return Find_Type<T_STRUCTURE>(structures,index);}

    template<class T_STRUCTURE> const T_STRUCTURE
    Find_Structure(const int index=0) const
    {return Find_Type<T_STRUCTURE>(structures,index);}

    bool simulate; // TODO: use one of those per fragment

    // These are pruned for MPI.
    ARRAY<int> simulated_particles;
    ARRAY<int> dynamic_particles;

    BINDING_LIST<TV>& binding_list;
    SOFT_BINDINGS<TV>& soft_bindings;

    MPI_SOLIDS<TV>* mpi_solids;

    bool print_diagnostics;
    bool print_residuals;
    int iterations_used_diagnostic;

    DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV>* deformables_example_forces_and_velocities;
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& triangle_repulsions_and_collisions_geometry;
    TRIANGLE_REPULSIONS<TV>& triangle_repulsions;
    TRIANGLE_COLLISIONS<TV>& triangle_collisions;
    DEFORMABLE_OBJECT_COLLISIONS<TV>& collisions;
    ARRAY<COLLISION_PENALTY_FORCES<TV>*> collision_penalty_forces;
    bool use_embedded_collisions;
    bool use_nonembedded_self_collision; // TODO: have one of these per fragment
    bool check_stale;

    DEFORMABLE_BODY_COLLECTION(DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV>* deformables_example_forces_and_velocities_input,COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list);
    virtual ~DEFORMABLE_BODY_COLLECTION();

//#####################################################################
    void Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collisions_parameters);
    void Update_Collision_Penalty_Forces_And_Derivatives();
    void Read_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame);
    void Write_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const;
    void Read_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame);
    void Write_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const;
    void Update_Simulated_Particles();
    void Update_Simulated_Particles(DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities);
    void Set_Mpi_Solids(MPI_SOLIDS<TV>* mpi_solids);

    void Read(const STREAM_TYPE,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,const bool read_from_every_process);
    void Write(const STREAM_TYPE,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,const bool write_from_every_process) const;
//#####################################################################
};
}
#endif
