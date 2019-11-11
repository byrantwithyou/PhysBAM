//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_SLENDER_ROD_COLLECTION
//#####################################################################
#ifndef __RIGID_SLENDER_ROD_COLLECTION__
#define __RIGID_SLENDER_ROD_COLLECTION__

#include <Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_FORWARD.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>

namespace PhysBAM{

template<class TV> class RIGID_SLENDER_ROD_PARTICLES;
template<class TV> class SLENDER_ROD_FORCES;
class VIEWER_DIR;

template<class TV>
class RIGID_SLENDER_ROD_COLLECTION
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_TYPED_READ_WRITE;

    RIGID_SLENDER_ROD_PARTICLES<TV>& rigid_slender_rod_particles;
    ARRAY<SLENDER_ROD_FORCES<TV>*> slender_rod_forces;

    ARRAY<int> simulated_rigid_rod_particles;
    ARRAY<int> dynamic_rigid_rod_particles;
    bool print_energy=false;

    RIGID_SLENDER_ROD_COLLECTION();
    RIGID_SLENDER_ROD_COLLECTION(const RIGID_SLENDER_ROD_COLLECTION&) = delete;
    void operator=(const RIGID_SLENDER_ROD_COLLECTION&) = delete;
    ~RIGID_SLENDER_ROD_COLLECTION();
    void Read(const VIEWER_DIR& viewer_dir);
    void Write(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const;
    void Update_Simulated_Particles();

    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time,bool transpose=false) const;

    void Update_Position_Based_State(const T time);
    void Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const;
    void Print_Energy(const T time,const int step) const;
    T CFL_Rigid_Rod(const bool verbose_dt);
    int Add_Force(SLENDER_ROD_FORCES<TV>* force);

    void Update_Kinematic_Particles();
    bool Register_Analytic_Replacement_Structure(const std::string& filename,const T scaling_factor,STRUCTURE<TV>* structure); // passing in zero skips reading
    bool Find_Or_Read_Structure(ARRAY<int>& structure_id,const std::string& filename,const T scaling_factor,const TV& center);
    void Destroy_Unreferenced_Geometry();
};
}
#endif
