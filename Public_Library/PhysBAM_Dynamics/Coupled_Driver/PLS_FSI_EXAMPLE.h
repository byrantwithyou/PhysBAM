//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLS_FSI_EXAMPLE
//#####################################################################
#ifndef __PLS_FSI_EXAMPLE__
#define __PLS_FSI_EXAMPLE__

#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION_CALLBACKS.h>
#include <PhysBAM_Fluids/PhysBAM_Fluids/Fluids/FLUID_COLLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/INCOMPRESSIBLE_COLLISIONS_FORWARD.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class SOLIDS_FLUIDS_PARAMETERS;
template<class TV> class SOLID_BODY_COLLECTION;
template<class T_GRID> class LAPLACE_UNIFORM;
template<class TV> class KANG_POISSON_VISCOSITY;

template<class TV_input>
class PLS_FSI_EXAMPLE:public EXAMPLE<TV_input>,public EXAMPLE_FORCES_AND_VELOCITIES<TV_input>,public SOLIDS_EVOLUTION_CALLBACKS<TV_input>,public SOLIDS_FLUIDS_CALLBACKS<TV_input>,
                      public LEVELSET_CALLBACKS<GRID<TV_input> >,public FLUIDS_PARAMETERS_CALLBACKS<GRID<TV_input> >,public NONCOPYABLE
{
    typedef TV_input TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MATRIX_POLICY<TV>::TRANSFORMATION_MATRIX T_TRANSFORMATION_MATRIX;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP T_FACE_LOOKUP;typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<TV> > T_FACE_LOOKUP_COLLIDABLE;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<GRID<TV> >::AVERAGING T_AVERAGING;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
public:
    typedef EXAMPLE<TV> BASE;using BASE::parse_args;using BASE::Set_Write_Substeps_Level;
    using BASE::output_directory;using BASE::first_frame;using BASE::restart;using BASE::Write_Frame_Title;using BASE::stream_type;
    using FLUIDS_PARAMETERS_CALLBACKS<GRID<TV> >::Get_Source_Reseed_Mask;
    using FLUIDS_PARAMETERS_CALLBACKS<GRID<TV> >::Get_Source_Velocities;using FLUIDS_PARAMETERS_CALLBACKS<GRID<TV> >::Get_Object_Velocities; // silence -Woverloaded-virtual

protected:
    T minimum_collision_thickness; // needed for ray tracing, etc.
public:
    SOLIDS_PARAMETERS<TV>& solids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    SOLIDS_EVOLUTION<TV>* solids_evolution; // defaults to newmark
    FLUIDS_PARAMETERS_UNIFORM<GRID<TV> > fluids_parameters;
    FLUID_COLLECTION<TV> fluid_collection;
    int resolution;
    int convection_order;
    bool use_pls_evolution_for_structure;
    bool two_phase;
    bool use_kang;
    bool print_matrix;
    bool test_system;
    KANG_POISSON_VISCOSITY<TV>* kang_poisson_viscosity;
    T m,s,kg;
    bool opt_skip_debug_data,opt_solidscg,opt_solidscr,opt_solidssymmqmr;

    PLS_FSI_EXAMPLE(const STREAM_TYPE stream_type,const int number_of_regions);
    virtual ~PLS_FSI_EXAMPLE();

    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::m> >& V_levelset,const T time) const PHYSBAM_OVERRIDE
    {V_levelset=fluid_collection.incompressible_fluid_collection.face_velocities;}

    void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time) PHYSBAM_OVERRIDE
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    void Set_Minimum_Collision_Thickness(const T minimum_collision_thickness_input=1e-6)
    {minimum_collision_thickness=minimum_collision_thickness_input;}

//#####################################################################
    void Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision=true,bool add_coupling=true);
    void Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision=true,bool add_coupling=true);
    void Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,bool add_collision=true,bool add_coupling=true);
    template<class GEOMETRY> void Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity);
    template<class GEOMETRY> void Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity,
        const ARRAY<bool,FACE_INDEX<TV::m> >& invalid_mask);
    template<class GEOMETRY> void Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source);
    template<class GEOMETRY> void Adjust_Phi_With_Source(const GEOMETRY& source,const int region,const T_TRANSFORMATION_MATRIX& world_to_source);
    void Revalidate_Fluid_Scalars();
    void Revalidate_Phi_After_Modify_Levelset();
    void Revalidate_Fluid_Velocity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Get_Object_Velocities(LAPLACE_UNIFORM<GRID<TV> >* elliptic_solver,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<GRID<TV> >& levelset_multiple,ARRAY<T,FACE_INDEX<TV::m> >& V_levelset,const T time) const PHYSBAM_OVERRIDE;
    void Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Read_Output_Files_Fluids(const int frame);
    void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
    void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time) PHYSBAM_OVERRIDE;
    void Log_Parameters() const PHYSBAM_OVERRIDE;
    void Register_Options() PHYSBAM_OVERRIDE;
    void Parse_Options() PHYSBAM_OVERRIDE;
    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,
        const T time) PHYSBAM_OVERRIDE;
    virtual void Update_Fluid_Parameters(const T dt,const T time);
//#####################################################################
    virtual void Post_Initialization(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Frame(const int frame){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Frame(const int frame){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Substep(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Substep(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();} // time at start of substep
    // solids
    virtual void Initialize_Bodies() {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Read_Output_Files_Solids(const int frame);
    // fluids
    virtual void Extrapolate_Phi_Into_Objects(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Phi(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual bool Adjust_Phi_With_Sources(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN(); return false; }
    virtual void Adjust_Phi_With_Objects(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Initialize_Velocities(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Setup_Initial_Refinement(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Initialize_Advection(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Clamp_Velocities(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    // melting
    virtual void Melting_Substep(const T dt,const T time){}
    virtual void Modify_Fluid_For_Melting(const T dt,const T time){}
    virtual void Update_Melting_Substep_Parameters(const T dt,const T time){}
    void Parse_Late_Options() PHYSBAM_OVERRIDE;
    template<class T_MPI> void Adjust_Output_Directory_For_MPI(const T_MPI mpi);
    virtual void Set_Boundary_Conditions_Callback(ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& psi_D_value,
        ARRAY<T,FACE_INDEX<TV::dimension> >& psi_N_value) const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    void Set_Boundary_Conditions(ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& psi_D_value,
        ARRAY<T,FACE_INDEX<TV::dimension> >& psi_N_value) const;
};
}
#endif
