//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_FORCE_COLLECTION
//#####################################################################
#ifndef __RIGID_FORCE_COLLECTION__
#define __RIGID_FORCE_COLLECTION__

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
template<class TV,class ID> class STRUCTURE_LIST;
template<class TV> class STRUCTURE;

template<class TV>
class RIGID_FORCE_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_TYPED_READ_WRITE;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    ARRAY<RIGIDS_FORCES<TV>*> rigids_forces;
    bool print_energy;
    
    RIGID_FORCE_COLLECTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    virtual ~RIGID_FORCE_COLLECTION();

    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T scale,const T time) const;

    void Update_Position_Based_State(const T time);
    void Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const;
    void Print_Energy(const T time,const int step) const;
    T CFL_Rigid(const RIGID_BODY_EVOLUTION_PARAMETERS<TV>& rigid_body_evolution_parameters,const bool verbose_dt);
    int Add_Force(RIGIDS_FORCES<TV>* force);
};
}
#endif
