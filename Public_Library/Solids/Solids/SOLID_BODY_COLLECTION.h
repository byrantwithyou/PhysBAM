//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_BODY_COLLECTION
//#####################################################################
#ifndef __SOLID_BODY_COLLECTION__
#define __SOLID_BODY_COLLECTION__

#include <Core/Utilities/Find_Type.h>
#include <Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <Solids/Forces_And_Torques/SOLIDS_FORCES.h>
namespace PhysBAM{

template<class TV> class DEFORMALBLE_OBJECT_COLLISIONS;
template<class TV> class SOLIDS_PARAMETERS;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGID_SLENDER_ROD_COLLECTION;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class RIGID_BODY_CLUSTER_BINDINGS;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class RIGIDS_FORCES;
template<class TV> class GENERALIZED_VELOCITY;
class VIEWER_DIR;

template<class TV>
class SOLID_BODY_COLLECTION
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename RIGIDS_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_RIGID;
    typedef typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_DEFORMABLE;
public:
    COLLISION_BODY_COLLECTION<TV>& collision_body_list;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    RIGID_SLENDER_ROD_COLLECTION<TV>& rigid_slender_rod_collection;
    ARRAY<SOLIDS_FORCES<TV>*> solids_forces;
    T cfl_number;
    T cfl_elastic,cfl_damping;
private:
    ARRAY<T_FREQUENCY_DEFORMABLE> frequency; // hertz for deformable CFL
    ARRAY<T_FREQUENCY_RIGID> rigid_frequency; // hertz for rigid CFL
public:
    bool implicit_damping;
    bool print_diagnostics;
    bool print_residuals;
    bool print_energy;
    bool simulate;
    int iterations_used_diagnostic;

    SOLID_BODY_COLLECTION(DEFORMABLE_PARTICLES<TV>* particles=0);
    SOLID_BODY_COLLECTION(const SOLID_BODY_COLLECTION&) = delete;
    void operator=(const SOLID_BODY_COLLECTION&) = delete;
    virtual ~SOLID_BODY_COLLECTION();

    void Print_Diagnostics(const bool print_diagnostics_input=true)
    {print_diagnostics=print_diagnostics_input;}

    void Print_Residuals(const bool print_residuals_input=true)
    {print_residuals=print_residuals_input;}

    void Set_CFL_Number(const T cfl_number_input=.5)
    {cfl_number=cfl_number_input;for(int i=0;i<solids_forces.m;i++) solids_forces(i)->Set_CFL_Number(cfl_number_input);
    deformable_body_collection.Set_CFL_Number(cfl_number);}

    void Set_Implicit_Damping(const bool implicit_damping_input=true)
    {implicit_damping=implicit_damping_input;}

    template<class T_FORCE> T_FORCE
    Find_Force(const int index=0)
    {return Find_Type<T_FORCE>(solids_forces,index);}

    template<class T_FORCE> const T_FORCE
    Find_Force(const int index=0) const
    {return Find_Type<T_FORCE>(solids_forces,index);}

//#####################################################################
    void Add_All_Forces(GENERALIZED_VELOCITY<TV>& F,const T time,const bool damping_only=false);
    int Add_Force(SOLIDS_FORCES<TV>* force);
    int Add_Force(DEFORMABLES_FORCES<TV>* force);
    int Add_Force(RIGIDS_FORCES<TV>* force);
    void Update_CFL();
    T CFL(const bool verbose=false);
    T CFL_Elastic_And_Damping();
    T CFL_Elastic();
    T CFL_Damping();
    T CFL_Strain_Rate();
    void Update_Simulated_Particles();
    void Delete_Forces();
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian);
    void Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const;
    void Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time) const;
    void Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,
        GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose=false) const;
    void Force_Differential(const GENERALIZED_VELOCITY<TV>& dX_full,GENERALIZED_VELOCITY<TV>& dF_full,const T time) const; 
    void Enforce_Definiteness(const bool enforce_definiteness_input=true);
    void Compute_Linear_Momentum(TV& linear_momentum) const;
    void Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const;
    TV Compute_Momentum() const;
    void Print_Energy(const T time,const int step) const;
    void Read(const VIEWER_DIR& viewer_dir,const bool static_variables_every_frame=false,
        ARRAY<int>* needs_init=0,ARRAY<int>* needs_destroy=0);
    void Write(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir,const bool static_variables_every_frame=false,const bool write_from_every_process=true) const;
    GENERALIZED_VELOCITY<TV>& New_Generalized_Velocity() const;
//#####################################################################
};
}
#endif
