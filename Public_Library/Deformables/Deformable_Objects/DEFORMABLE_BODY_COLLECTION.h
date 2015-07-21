//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonthan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_BODY_COLLECTION
//#####################################################################
#ifndef __DEFORMABLE_BODY_COLLECTION__
#define __DEFORMABLE_BODY_COLLECTION__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Tools/Utilities/Find_Type.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class BINDING_LIST;
template<class TV> class DEFORMABLE_GEOMETRY_COLLECTION;
template<class TV> class COLLISION_BODY_COLLECTION;
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

template<class TV>
class DEFORMABLE_BODY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_DEFORMABLE;
public:
    DEFORMABLE_PARTICLES<TV>& particles;
    ARRAY<STRUCTURE<TV>*> structures;

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
    ARRAY<DEFORMABLES_FORCES<TV>*> deformables_forces;

    MPI_SOLIDS<TV>* mpi_solids;
    ARRAY<T_FREQUENCY_DEFORMABLE> frequency; // hertz for deformable CFL
    T cfl_number;
    T cfl_elastic,cfl_damping;
    bool implicit_damping;

    bool print_diagnostics;
    bool print_residuals;
    bool print_energy;
    int iterations_used_diagnostic;

    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& triangle_repulsions_and_collisions_geometry;
    TRIANGLE_REPULSIONS<TV>& triangle_repulsions;
    TRIANGLE_COLLISIONS<TV>& triangle_collisions;
    DEFORMABLE_OBJECT_COLLISIONS<TV>& collisions;
    ARRAY<COLLISION_PENALTY_FORCES<TV>*> collision_penalty_forces;
    bool use_embedded_collisions;
    bool use_nonembedded_self_collision; // TODO: have one of these per fragment
    bool check_stale;
    bool own_particles,own_collision_body_collection;

    DEFORMABLE_BODY_COLLECTION(DEFORMABLE_PARTICLES<TV>* particles,COLLISION_BODY_COLLECTION<TV>* collision_body_list);
    virtual ~DEFORMABLE_BODY_COLLECTION();

    template<class T_FORCE> T_FORCE
    Find_Force(const int index=0)
    {return Find_Type<T_FORCE>(deformables_forces,index);}

    template<class T_FORCE> const T_FORCE
    Find_Force(const int index=0) const
    {return Find_Type<T_FORCE>(deformables_forces,index);}

    void Test_Forces(const T time)
    {Test_Energy(time);Test_Force_Derivatives(time);}

//#####################################################################
    void Set_CFL_Number(const T cfl_number_input=.5);
    int Add_Force(DEFORMABLES_FORCES<TV>* force);
    void Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collisions_parameters);
    void Update_Collision_Penalty_Forces_And_Derivatives();
    void Read_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame);
    void Write_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const;
    void Read_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame);
    void Write_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const;
    void Update_Simulated_Particles();
    void Set_Mpi_Solids(MPI_SOLIDS<TV>* mpi_solids);
    void Update_CFL();
    T CFL(const bool verbose=false);
    T CFL_Elastic_And_Damping();
    T CFL_Elastic();
    T CFL_Damping();
    T CFL_Strain_Rate();

    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian);
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F_full,const T time) const;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T time) const;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T time) const;

    void Test_Energy(const T time);
    void Test_Force_Derivatives(const T time);
    void Read(const STREAM_TYPE,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,const bool read_from_every_process);
    void Write(const STREAM_TYPE,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,const bool write_from_every_process) const;
//#####################################################################
};
}
#endif
