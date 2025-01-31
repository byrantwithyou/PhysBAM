//#####################################################################
// Copyright 2006-2007, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_ZERO_LENGTH_SPRINGS
//#####################################################################
#ifndef __IMPLICIT_ZERO_LENGTH_SPRINGS__
#define __IMPLICIT_ZERO_LENGTH_SPRINGS__

#include <Core/Arrays/ARRAY.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <climits>
namespace PhysBAM{

template<class TV>
class IMPLICIT_ZERO_LENGTH_SPRINGS:public DEFORMABLES_FORCES<TV>,public SPRINGS_TAG
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    SEGMENT_MESH& segment_mesh;
    ARRAY<T> stiffness;
    T constant_stiffness;
    ARRAY<T> damping;
    T constant_damping;
protected:
    ARRAY<int> force_segments;

public:
    IMPLICIT_ZERO_LENGTH_SPRINGS(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh_input);
    virtual ~IMPLICIT_ZERO_LENGTH_SPRINGS();

    virtual void Set_Stiffness(const T stiffness_input)
    {constant_stiffness=stiffness_input;stiffness.Clean_Memory();}

    virtual void Set_Stiffness(ARRAY_VIEW<const T> stiffness_input)
    {constant_stiffness=0;stiffness=stiffness_input;}

    virtual void Set_Damping(const T damping_input)
    {constant_damping=damping_input;damping.Clean_Memory();}

    virtual void Set_Damping(ARRAY_VIEW<const T> damping_input)
    {constant_damping=0;damping=damping_input;}

    void Enforce_Definiteness(const bool enforce_definiteness_input) override
    {}

    T CFL_Strain_Rate() const override
    {return FLT_MAX;}

    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override
    {}

    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override
    {}

//#####################################################################
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    virtual void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass is already defined
    void Set_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction); // 1 is critically damped
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    T Potential_Energy(const T time) const override;
//#####################################################################
};

template<class T,class TV> IMPLICIT_ZERO_LENGTH_SPRINGS<TV>*
Create_Edge_Zero_Length_Springs(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const T stiffness=2e3,const T overdamping_fraction=1,const bool verbose=true);
}
#endif
