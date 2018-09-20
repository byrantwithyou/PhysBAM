//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_EXAMPLE_RB__
#define __MPM_EXAMPLE_RB__
#include <Core/Data_Structures/CHAINED_ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Matrices/MATRIX.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Rigids/Forces_And_Torques/MOVE_RIGID_BODY_DIFF.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <functional>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class DEBUG_PARTICLES;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class SOLIDS_FORCES;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class GATHER_SCATTER;
template<class TV> class IMPLICIT_OBJECT;
template<class T>  class KRYLOV_VECTOR_BASE;
template<class TV> class MPM_FORCE_HELPER;
template<class T>  class MPM_KRYLOV_VECTOR_RB;
template<class TV> class MPM_PARTICLES;
template<class TV> class PARTICLE_GRID_FORCES;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class MPM_PLASTICITY_MODEL;
template<class TV> class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION;
template<class T> class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION;
template<class T> class RIGID_PENALTY_WITH_FRICTION;
template<class TV> class PENALTY_FORCE_COLLECTION;

template<class TV>
class MPM_EXAMPLE_RB
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
    GRID<TV> grid;
    STREAM_TYPE stream_type;
    MPM_PARTICLES<TV>& particles;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    DEBUG_PARTICLES<TV>& debug_particles;
    ARRAY<int> simulated_particles;
    ARRAY<bool> particle_is_simulated;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > Dp_inv;
    ARRAY<T,TV_INT> mass,volume;
    ARRAY<TV,TV_INT> location;
    ARRAY<TV,TV_INT> velocity,velocity_save,velocity_friction_save;
    ARRAY<int> valid_grid_indices;
    ARRAY<TV_INT> valid_grid_cell_indices;
    ARRAY<PARTICLE_GRID_FORCES<TV>*> forces;
    ARRAY<MPM_PLASTICITY_MODEL<TV>*> plasticity_models;
    ARRAY<DEFORMABLES_FORCES<TV>*>& lagrangian_forces;
    PARTICLE_GRID_WEIGHTS<TV>* weights=0;
    GATHER_SCATTER<TV>& gather_scatter;
    ARRAY<MPM_COLLISION_OBJECT<TV>*> collision_objects;
    mutable ARRAY<TV> lagrangian_forces_V,lagrangian_forces_F;
    MPM_FORCE_HELPER<TV>& force_helper;

    T initial_time=0;
    int last_frame=100;
    std::string frame_title;
    int write_substeps_level=-1;
    int substeps_delay_frame=-1;
    std::string output_directory="output";
    std::string data_directory="../../Public_Data";
    std::string test_output_prefix;
    bool use_test_output=false;
    T mass_contour=-1;
    int restart=0;
    T dt=0,time=0,frame_dt=(T)1/24;
    T min_dt=0,max_dt=frame_dt;
    int ghost=3;
    bool use_affine=true;
    bool use_midpoint=false;
    bool use_symplectic_euler=false;
    bool lag_Dp=false; // Bp will actually store Cp
    bool use_oldroyd=false;
    bool print_stats=false;
    bool only_write_particles=false;
    T flip=0;
    T cfl=1;
    T cfl_F=(T).5;
    T cfl_sound=(T).9;
    T inv_Wi=0;
    bool use_strong_cfl=false;
    bool use_sound_speed_cfl=false;
    int reflection_bc_flags=0;
    
    T newton_tolerance=1;
    int newton_iterations=100;
    T solver_tolerance=(T).5;
    int solver_iterations=1000;
    T angle_tol=(T).01;
    bool test_diff=false;
    bool test_system=false;
    int threads=1;
    bool use_gradient_magnitude_objective=false;
    bool debug_newton=false;
    
    TV last_linear_momentum;
    typename TV::SPIN last_angular_momentum;
    T last_te=0;
    T last_grid_ke=0;
    bool output_structures_each_frame=false;
    T quad_F_coeff=0;
    bool asymmetric_system=false;

    ARRAY<bool> rigid_body_is_simulated;
    HASHTABLE<PAIR<int,int>,TV> stored_contacts_rr;
    HASHTABLE<PAIR<int,TV_INT>,TV> stored_contacts_rm;
    T contact_factor=1; // omega
    T impulse_interpolation=1; // lambda
    T min_impulse_change=1e-4;
    
    bool pairwise_collisions=false;
    bool projected_collisions=false;
    int collision_iterations=5;

    T rd_penalty_stiffness=0;
    T rd_penalty_friction=0.3;
    bool use_rd=false,use_rr=false,use_di=false;
    bool use_rd_k=false,use_rr_k=false,use_di_k=false;
    T rd_k=0,rr_k=0,di_k=0;
    bool use_rd_mu=false,use_rr_mu=false,use_di_mu=false;
    T rd_mu=0,rr_mu=0,di_mu=0;
    PENALTY_FORCE_COLLECTION<TV>* pfd=0;
    ARRAY<MOVE_RIGID_BODY_DIFF<TV> > move_rb_diff;
    bool use_bisection=false;
    bool use_rd_ccd=false,use_rr_ccd=false,use_di_ccd=false;
    
    MPM_EXAMPLE_RB(const STREAM_TYPE stream_type_input);
    MPM_EXAMPLE_RB(const MPM_EXAMPLE_RB&) = delete;
    void operator=(const MPM_EXAMPLE_RB&) = delete;
    virtual ~MPM_EXAMPLE_RB();

    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    std::function<void(int frame)> begin_frame;
    std::function<void(int frame)> end_frame;
    std::function<void(T time)> begin_time_step;
    std::function<void(T time)> end_time_step;

    void Capture_Stress();
    void Precompute_Forces(const T time,const T dt,const bool update_hessian);
    T Potential_Energy(const T time) const;
    void Add_Forces(ARRAY<TV,TV_INT>& F,ARRAY<TWIST<TV> >& RF,const T time) const;
    void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,
        ARRAY<TWIST<TV> >& RF,const ARRAY<TWIST<TV> >& RV,
        const T time,bool transpose=false) const;
    void Update_Lagged_Forces(const T time) const;
    int Add_Force(PARTICLE_GRID_FORCES<TV>& force);
    int Add_Force(DEFORMABLES_FORCES<TV>& force);
    int Add_Force(SOLIDS_FORCES<TV>& force);
    void Set_Weights(int order);

    TV Total_Particle_Linear_Momentum() const;
    TV Total_Grid_Linear_Momentum(const ARRAY<TV,TV_INT>& u) const;
    typename TV::SPIN Total_Grid_Angular_Momentum(T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0) const;
    typename TV::SPIN Total_Particle_Angular_Momentum() const;
    T Total_Grid_Kinetic_Energy(const ARRAY<TV,TV_INT>& u) const;
    T Total_Particle_Kinetic_Energy() const;
    T Average_Particle_Mass() const;
    template<class S> void Reflection_Boundary_Condition(ARRAY<S,TV_INT>& u,bool flip_sign) const;
//#####################################################################
};
}
#endif
