//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_LINEAR_SPRINGS
//#####################################################################
#ifndef __RIGID_LINEAR_SPRINGS__
#define __RIGID_LINEAR_SPRINGS__

#include <Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGID_BODY;
template<class TV>
class RIGID_LINEAR_SPRINGS:public RIGIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef RIGIDS_FORCES<TV> BASE;
    using BASE::Invalidate_CFL;using BASE::cfl_number;using BASE::rigid_body_collection;
    using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;
    typedef typename FORCE_ELEMENTS::ITERATOR SEGMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

public:
    SEGMENT_MESH segment_mesh;
    ARRAY<T> youngs_modulus; // units of force (i.e. force per unit strain)
    ARRAY<T> restlength,visual_restlength,current_lengths; // visual restlength corresponds to length between particles; restlength may be larger than this to avoid zero/small restlength
    ARRAY<T> damping; // units of force*time (i.e. force per unit strain rate)
    mutable ARRAY<VECTOR<T,2> > strains_of_segment; // VECTOR<T,2>(strain_rate, strain)
    ARRAY<VECTOR<TV,2> > attachment_radius;

    struct STATE{
        STATE()
            :coefficient(0)
        {}

        VECTOR<int,2> nodes; // copy of nodes to reduce associativity needed in cache
        T coefficient;
        TV direction;
        VECTOR<TV,2> r;
    };
    ARRAY<STATE> states;
public:
    FORCE_ELEMENTS force_segments;
    RIGID_LINEAR_SPRINGS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    virtual ~RIGID_LINEAR_SPRINGS();

    void Enforce_Definiteness(const bool enforce_definiteness_input) override
    {}

//#####################################################################
    void Add_Spring(int body0,int body1,const TV& r1,const TV& r2);
    TV Attachment_Location(int s,int b) const;
    TV Endpoint_Velocity(int s,int b) const;
    T Spring_Length(int s) const;
    void Set_Restlengths();
    void Set_Stiffness(int b,T stiffness);
    void Set_Damping(int b,T damp);
    void Set_Restlength(int b,T length=-1,T visual=-1);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated) override;
    void Set_Overdamping_Fraction(int b,const T overdamping_fraction=1); // 1 is critically damped
    void Update_Position_Based_State(const T time) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    void Add_Force(ARRAY_VIEW<TWIST<TV> > rigid_F,const STATE& state,const TV& force) const;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
    T Average_Restlength() const;
    void Print_Restlength_Statistics() const;
    void Print_Deformation_Statistics() const;
    T Maximum_Compression_Or_Expansion_Fraction(int* index=0) const;
    T Potential_Energy(int s,const T time) const;
    T Potential_Energy(const T time) const override;
    T Compute_Total_Energy(const T time) const;
    T Effective_Impulse_Factor(int s,int b) const;
    T Effective_Impulse_Factor(int s) const;
    const RIGID_BODY<TV>& Body(int s,int b) const;
    RIGID_BODY<TV>& Body(int s,int b);
//#####################################################################
};
}
#endif
