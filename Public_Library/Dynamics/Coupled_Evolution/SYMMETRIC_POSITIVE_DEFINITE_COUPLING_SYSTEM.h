//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM
//#####################################################################
#ifndef __SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM__
#define __SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM__
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>

namespace PhysBAM{
template<class TV> class BACKWARD_EULER_SYSTEM;
template<class TV> class GENERALIZED_MASS;
template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class MATRIX_FLUID_GRADIENT_BASE;
template<class TV> class MATRIX_FLUID_POISSON;
template<class TV> class MATRIX_SOLID_INTERPOLATION_BASE;
template<class TV> class MATRIX_SOLID_FORCES;
template<class TV> class FLUID_TO_SOLID_INTERPOLATION_BASE;
template<class TV> class MATRIX_VISCOUS_FORCES;
template<class TV> class MATRIX_FLUID_INTERPOLATION_BASE;
template<class TV> class COUPLED_SYSTEM_VECTOR;
template<class TV> class GENERALIZED_FLUID_MASS;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class BOUNDARY_CONDITION_COLLECTION;
template<class TV> class FLUID_COLLECTION;
template<class TV> class IMPLICIT_BOUNDARY_CONDITION_COLLECTION;
template<class TV> class MPI_UNIFORM_GRID;
template<class TV> class MPI_SOLID_FLUID;

template<class TV>
class SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef COUPLED_SYSTEM_VECTOR<TV> VECTOR_T;
    typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
protected:
    using BASE::use_preconditioner;using BASE::Test_System;

    GENERALIZED_MASS<TV>* solid_mass; // M
    const ARRAY<bool,TV_INT>& outside_fluid;
    const ARRAY<bool,FACE_INDEX<TV::m> >& inactive_faces;
    const ARRAY<T,TV_INT>& density;
    const ARRAY<TV,TV_INT>& centered_velocity;
    ARRAY<T> constrained_fluid_velocity;
    ARRAY<T,FACE_INDEX<TV::m> > beta_face;
    ARRAY<T> constrained_beta_face;
    const ARRAY<T,TV_INT>& one_over_rho_c_squared;
    ARRAY<T> one_over_rho_c_squared_flat;
    mutable ARRAY<T> full_precondition_in,full_precondition_out;

public:
    BACKWARD_EULER_SYSTEM<TV>* solid_system;
    COLLISION_AWARE_INDEX_MAP<TV>& index_map;
    MATRIX_SOLID_INTERPOLATION_BASE<TV>* solid_interpolation; // J
    MATRIX_FLUID_GRADIENT_BASE<TV>* fluid_gradient; // G
    MATRIX_FLUID_POISSON<TV>* fluid_poisson; // top left block
    MATRIX_SOLID_FORCES<TV>* solid_forces; // C
    MATRIX_FLUID_INTERPOLATION_BASE<TV>* fluid_interpolation; // W
    GENERALIZED_FLUID_MASS<TV>* fluid_mass; // beta
    MATRIX_VISCOUS_FORCES<TV>* fluid_viscous_forces;
    FLUID_TO_SOLID_INTERPOLATION_BASE<TV>* fluid_to_solid_interpolation; // H
protected:

    const bool leakproof_solve;
    bool fully_implicit;
    T dt;

    mutable ARRAY<T,COUPLING_CONSTRAINT_ID> coupling_faces;
    mutable ARRAY<TV> temporary_velocities;
    mutable ARRAY<TWIST<TV> > temporary_twists;
    mutable ARRAY<TWIST<TV> > temporary_rod_twists;
    mutable ARRAY<T> pressure;
    mutable ARRAY<T> temporary_faces;
    mutable ARRAY<T,COUPLING_CONSTRAINT_ID> temporary_lambdas;
    mutable ARRAY<T> temporary_viscous_velocities;
    VECTOR_T tolerances;
    const LEVELSET<TV>* levelset;
    ARRAY<T,COUPLING_CONSTRAINT_ID> lambda_diagonal_preconditioner;

public:
    mutable ARRAY<TV> pressure_impulses;
    mutable ARRAY<TWIST<TV> > pressure_impulses_twist;
    mutable ARRAY<TWIST<TV> > pressure_impulses_rod_twist;

    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_matrix;
    bool print_rhs;
    bool print_each_matrix;
    bool print_index_map;
    static int solve_id;

    bool use_viscous_forces;
    T surface_tension_coefficient;

    bool solid_node;
    bool fluid_node;
    bool use_full_ic;
    SPARSE_MATRIX_FLAT_MXN<T> full_matrix;
    MPI_SOLID_FLUID<TV>* mpi_solid_fluid;
    MPI_UNIFORM_GRID<TV>* mpi_grid;
    mutable ARRAY<T,TV_INT> pressure_on_grid;
    ARRAY<T,FACE_INDEX<TV::m> >* debug_velocity;
    GENERALIZED_VELOCITY<TV>* debug_generalized_velocity;

    SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM(const bool use_preconditioner_input,BACKWARD_EULER_SYSTEM<TV>* solid_system_input,
        IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>& boundary_condition_collection,UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,
        const SOLID_BODY_COLLECTION<TV>& solid_body_collection,
        const FLUID_COLLECTION<TV>& fluid_collection,const ARRAY<T,TV_INT>& density_input,const ARRAY<TV,TV_INT>& centered_velocity_input,
        const ARRAY<T,TV_INT>& one_over_rho_c_squared_input,const bool leakproof_solve_input,const bool using_slip,const bool using_viscosity,bool fully_implicit);

    virtual ~SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM();

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& bV) const override
    {if(run_self_tests) Test_Matrix();}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& bV) const override// TODO
    {}

//#####################################################################
    void Project(KRYLOV_VECTOR_BASE<T>& bV) const override;
    //void Project() const;
    void Get_Pressure(const VECTOR_T& V,ARRAY<T,TV_INT>& fluid_pressures) const;
    void Zero_Coupling_Faces_Values(ARRAY<T,FACE_INDEX<TV::m> >& face_array) const;
    void Apply_Lambda_To_Euler_State(const VECTOR_T& V,const ARRAY<T,COUPLING_CONSTRAINT_ID>& coupled_faces_solid_interpolated_velocity_n,
        const ARRAY<T,COUPLING_CONSTRAINT_ID>& coupled_faces_solid_interpolated_velocity_np1,ARRAY<TV_DIMENSION,TV_INT>& U) const;
    void Interpolate_Solid_Velocity_To_Coupled_Faces(const GENERALIZED_VELOCITY<TV>& solids_velocity,ARRAY<T,COUPLING_CONSTRAINT_ID>& coupled_faces_solid_interpolated_velocity);
    void Add_Dirichlet_Pressures_To_Velocity(const ARRAY<T,TV_INT>& pressure,ARRAY<T>& fluid_velocity_vector) const;
    void Add_Surface_Tension(ARRAY<T>& fluid_velocity_vector) const;
    void Print_Matrix(const VECTOR_T& vec) const;
    void Print_Each_Matrix(int n) const;
    T Residual_Linf_Norm(const VECTOR_T& x,const VECTOR_T& rhs) const;
    //void Apply(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF) const override;
    void Test_Matrix() const;
    void Test_Incompressibility(const ARRAY<T,FACE_INDEX<TV::m> >& fluid_velocity,const ARRAY<T>& constrained_fluid_velocity) const;
    void Test_Viscosity(const VECTOR_T& V,const VECTOR_T& B) const;
    void Set_Coupling_Faces(const int ghost_cells,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N) const;
    void Show_Constraints(ARRAY<bool,FACE_INDEX<TV::m> >& psi_N) const;
    void Check_Constraints(const GENERALIZED_VELOCITY<TV>& solids_rhs,const ARRAY<T,FACE_INDEX<TV::m> >& fluids_rhs,const ARRAY<T>& constrained_rhs) const;
    void Compute_Lambda_Diagonal_Preconditioner();
    void Compute_Full_Preconditioner();
    void Compute_Scatter_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& scatter_matrix);
    void Compute_Inverse_Mass_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& inverse_mass);
    void Gather(const KRYLOV_VECTOR_BASE<T>& bV,ARRAY<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& structure_velocity) const;
    void Scatter(const ARRAY<T>& fluid_velocity,const GENERALIZED_VELOCITY<TV>& structure_velocity,KRYLOV_VECTOR_BASE<T>& bF) const;
    void Inverse_Mass(ARRAY<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& structure_velocity) const;
    void Massless_Gather(const VECTOR_T& V,ARRAY<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& structure_velocity) const;
    void Massless_Scatter(const ARRAY<T>& fluid_velocity,const GENERALIZED_VELOCITY<TV>& structure_velocity,VECTOR_T& F) const;
    void Setup_Tolerances(const VECTOR_T& F,const ARRAY<T>& fluid_velocity,const GENERALIZED_VELOCITY<TV>& structure_velocity);
    void Setup_Initial_Guess(const VECTOR_T& F,VECTOR_T& V,const ARRAY<T,TV_INT>& p_advected) const;
    void Exchange_Pressure(ARRAY<T>& pressure) const;
    void Exchange_Coupled_Vector(VECTOR_T& V) const;
    void Exchange_Velocities(ARRAY<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& structure_velocity) const;
    void Set_MPI(MPI_SOLID_FLUID<TV>& mpi_solid_fluid_input,MPI_UNIFORM_GRID<TV>& mpi_grid_input);
    void Apply_Massless_Structure_Force_To_Fluid(ARRAY<T>& fluid_velocity,T time) const;
    void Dump_Substep(const ARRAY<T>& fluid_velocity,const char* name,int level=1) const;
    void Dump_Substep(const ARRAY<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity,const char* name,int level=1) const;
    void Fill_Extra_Velocities(ARRAY<T>& fluid_velocity_vector) const;
//#####################################################################
    void Resize_Coupled_System_Vector(VECTOR_T& b) const;
    void Apply_One_Sided_Interpolation_At_Coupling_Faces(const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N_domain_boundary,
        const bool use_one_sided_face_velocty_interpolation,ARRAY<T,FACE_INDEX<TV::m> >& fluids_velocity);
    void Compute(int ghost_cells,const T dt_input,const T current_velocity_time,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N_domain_boundary,const bool disable_thinshell,
        const bool use_one_sided_face_velocty_interpolation,ARRAY<T,FACE_INDEX<TV::m> >& fluids_velocity,T mu,bool use_second_order_cut_cell,const LEVELSET<TV>* levelset);
    void Set_Up_RHS(VECTOR_T& V,VECTOR_T& F,const GENERALIZED_VELOCITY<TV>& solids_velocity_star,const ARRAY<T,FACE_INDEX<TV::m> >& fluids_velocity_star,
        const ARRAY<T,TV_INT>& p_advected_over_rho_c_squared_dt,const ARRAY<T,TV_INT>& p_advected,const ARRAY<T,TV_INT>& fluid_pressures);
    virtual void Multiply(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF,bool transpose=false) const override;
    virtual double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bV0,const KRYLOV_VECTOR_BASE<T>& bV1) const override;
    virtual void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bR) const override;  // solve MR=V
    T Linf_Norm(const VECTOR_T& bR) const;
    virtual T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bR) const override;
    void Apply_Velocity_Update(const VECTOR_T& V,ARRAY<T,FACE_INDEX<TV::m> >& fluid_velocity,ARRAY<T,TV_INT>& fluid_pressures,GENERALIZED_VELOCITY<TV>& solid_velocity,
        GENERALIZED_VELOCITY<TV>& force_on_solid,bool want_solid,bool want_fluid) const;
    void Mark_Valid_Faces(ARRAY<bool,FACE_INDEX<TV::m> >& valid) const;
//#####################################################################
};
}
#endif
