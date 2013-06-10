//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FORCE_COLLECTION
//#####################################################################
#ifndef __SOLID_FORCE_COLLECTION__
#define __SOLID_FORCE_COLLECTION__

#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/SOLIDS_FORCES.h>
namespace PhysBAM{

template<class TV> class EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class DEFORMALBLE_OBJECT_COLLISIONS;
template<class TV> class SOLIDS_PARAMETERS;
template<class TV> class RIGID_FORCE_COLLECTION;
template<class TV> class DEFORMABLE_FORCE_COLLECTION;
template<class TV> class RIGID_FORCE_CLUSTER_BINDINGS;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class RIGIDS_FORCES;

template<class TV>
class SOLID_FORCE_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename RIGIDS_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_RIGID;
    typedef typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_DEFORMABLE;
public:
    DEFORMABLE_FORCE_COLLECTION<TV>& deformable_force_collection;
    RIGID_FORCE_COLLECTION<TV>& rigid_force_collection;
    ARRAY<SOLIDS_FORCES<TV>*> solids_forces;
    T cfl_number;
    T cfl_elastic,cfl_damping;
private:
    ARRAY<T_FREQUENCY_DEFORMABLE> frequency; // hertz for deformable CFL
    ARRAY<T_FREQUENCY_RIGID> rigid_frequency; // hertz for rigid CFL
public:
    bool implicit_damping;
    bool print_energy;

    SOLID_FORCE_COLLECTION(DEFORMABLE_FORCE_COLLECTION<TV>& deformable_force_collection,RIGID_FORCE_COLLECTION<TV>& rigid_force_collection);
    virtual ~SOLID_FORCE_COLLECTION();

    void Set_CFL_Number(const T cfl_number_input=.5)
    {cfl_number=cfl_number_input;for(int i=0;i<solids_forces.m;i++) solids_forces(i)->Set_CFL_Number(cfl_number_input);
    deformable_force_collection.Set_CFL_Number(cfl_number);}

    void Set_Implicit_Damping(const bool implicit_damping_input=true)
    {implicit_damping=implicit_damping_input;}

    template<class T_FORCE> T_FORCE
    Find_Force(const int index=0)
    {return Find_Type<T_FORCE>(solids_forces,index);}

    template<class T_FORCE> const T_FORCE
    Find_Force(const int index=0) const
    {return Find_Type<T_FORCE>(solids_forces,index);}

//#####################################################################
    void Add_All_Forces(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time,const bool damping_only=false);
    int Add_Force(SOLIDS_FORCES<TV>* force);
    int Add_Force(DEFORMABLES_FORCES<TV>* force);
    int Add_Force(RIGIDS_FORCES<TV>* force);
    void Update_CFL();
    T CFL(const bool verbose=false);
    T CFL_Elastic_And_Damping();
    T CFL_Elastic();
    T CFL_Damping();
    T CFL_Strain_Rate();
    void Disable_Finite_Volume_Damping();
    void Disable_Spring_Elasticity();
    void Delete_Forces();
    void Update_Position_Based_State(const T time,const bool is_position_update);
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T scale,
        const T time) const;
    void Force_Differential(ARRAY_VIEW<const TV> dX_full,ARRAY_VIEW<TV> dF_full,const T time) const;
    void Enforce_Definiteness(const bool enforce_definiteness_input=true);
    void Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const;
    void Print_Energy(const T time,const int step) const;
//#####################################################################
};
}
#endif
