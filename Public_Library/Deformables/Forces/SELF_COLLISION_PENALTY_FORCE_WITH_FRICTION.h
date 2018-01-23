//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION
//#####################################################################
#ifndef __SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION__
#define __SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <functional>
namespace PhysBAM{

template<class TV>
class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    typedef int HAS_TYPED_READ_WRITE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    using BASE::particles;

    enum class STATE
    {
        triangle,edge,point
    };
    
    ARRAY<T_SURFACE*> surfaces;
    T stiffness_coefficient=0;
    T friction=0;
    T trial_distance;
    
    struct DIFF_ENTRY
    {
        int e; // element of relaxed attachment point
        TV Y; // relaxed attachment point
        VECTOR<MATRIX<T,TV::m>,5> dYdI; // Dependence on X, Z, A, B, C
    };

    struct COLLISION_PAIR
    {
        int p; // colliding particle
        int s; // surface

        TV w0; // barycentric coords of original attachment point
        int e0; // element of original attachment point
        TV Y0; // original attachment point

        ARRAY<DIFF_ENTRY> diff_entry;
        TV w; // barycentric coords of relaxed attachment point
        VECTOR<MATRIX<T,TV::m>,5> dwdI; // Dependence on X, Z, A, B, C
        bool active;

        template<class RW> void Write(std::ostream& output) const
        {Write_Binary<RW>(output,p,s,w0,e0);}

        template<class RW> void Read(std::istream& input)
        {Read_Binary<RW>(input,p,s,w0,e0);}
    };

    ARRAY<COLLISION_PAIR> collision_pairs;
    HASHTABLE<PAIR<int,int> > hash; // p,s
    HASHTABLE<TV_INT,VECTOR<int,2> > object_from_element; // face -> (s,e)
    std::function<void()> get_candidates=0; // Call Add_Pair on collision candidates.

    SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input);
    virtual ~SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION();

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const override;
    void Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    T Potential_Energy(const T time) const override;
    void Relax_Attachment(int cp);
    void Update_Attachments_And_Prune_Pairs();
    void Add_Pair(int p,int s,const TV& w0,int e0);
    void Test_Relax(int cp);
    void Add_Surface(T_SURFACE& surface);
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
//#####################################################################
};
}
#endif
