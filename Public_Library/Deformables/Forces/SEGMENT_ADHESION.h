//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_ADHESION
//#####################################################################
#ifndef __SEGMENT_ADHESION__
#define __SEGMENT_ADHESION__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Deformables/Deformable_Objects/HAIR_ID.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{
class SEGMENT_MESH;
class VIEWER_DIR;

template<class TV>
class SEGMENT_ADHESION_SPRING_STATE
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    VECTOR<int,4> nodes;
    VECTOR<T,2> weights;
    TV normal;
    T distance;
    T damping;
    bool external;    

    SEGMENT_ADHESION_SPRING_STATE()
        :distance(0),damping(0),external(false)
    {}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,nodes,weights,normal,distance,damping);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,nodes,weights,normal,distance,damping);}

};

template<class TV>
class SEGMENT_ADHESION:public DEFORMABLES_FORCES<TV>,public SPRINGS_TAG
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT T_SEGMENTED_CURVE;
    enum WORKAROUND {d=TV::m};
    MPI_SOLIDS<TV>* mpi_solids;
public:
    typedef SEGMENT_ADHESION_SPRING_STATE<TV> SPRING_STATE;
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    SEGMENT_MESH& mesh;
    T_SEGMENTED_CURVE curve;
    T on_distance;
    T off_distance;
    T overdamping_fraction;
    T youngs_modulus;
    T restlength;
    int max_connections;
    HASHTABLE<VECTOR<int,2> > existing_pairs;
    ARRAY<HAIR_ID>& particle_to_spring_id;
    
    ARRAY<int> internal_segment_indices;
    ARRAY<int> external_segment_indices;
    SEGMENT_MESH internal_mesh,external_mesh;
    T_SEGMENTED_CURVE internal_curve,external_curve;

    typedef HASHTABLE<VECTOR<int,2>,SPRING_STATE> T_SPRING_HASH;
    T_SPRING_HASH *springs; // Internal springs and external springs owned by this processor
    ARRAY<VECTOR<int,2> > external_spring_segments;
    ARRAY<SPRING_STATE> internal_springs; // cache for Add_Velocity_Dependent_Forces
    ARRAY<SPRING_STATE> external_springs;
    
    ARRAY<PAIR<ARRAY<PAIR<T,int> >,bool> > segments_with_springs;
    HASHTABLE<VECTOR<int,4> >& intersecting_edge_edge_pairs;
    HASHTABLE<VECTOR<int,4> > default_intersecting_edge_edge_pairs;

    SEGMENT_ADHESION(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& mesh,ARRAY<HAIR_ID>& particle_to_spring_id,HASHTABLE<VECTOR<int,4> >& intersecting_edge_edge_pairs);
    SEGMENT_ADHESION(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& mesh,ARRAY<HAIR_ID>& particle_to_spring_id);

    virtual ~SEGMENT_ADHESION()
    {
        delete springs;
    }

//#####################################################################
    void Update_Hierarchy();
    void Update_Springs(const bool search_hierarchy);
    void Update_Collisions_List(); 
    void Update_Partitions(bool restart,MPI_SOLIDS<TV>* mpi_solids,const VIEWER_DIR& viewer_dir);
    void Set_Parameters(const T youngs_modulus_input,const T overdamping_fraction_input,const T on_distance,const T off_distance, const int max_connections_input);
    void Set_Restlength(const T restlength);
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override;
    void Write_State(STREAM_TYPE type,const std::string& filename);
    void Read_State(const std::string& filename);
//#####################################################################
};
}
#endif
