//#####################################################################
// Copyright 2007, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_PROJECTION_UNIFORM  
//#####################################################################
#ifndef __EULER_PROJECTION_UNIFORM__
#define __EULER_PROJECTION_UNIFORM__

#include <Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <Compressible/Euler_Equations/BOUNDARY_OBJECT_REFLECTION.h>
#include <Compressible/Euler_Equations/EULER_PROJECTION.h>
namespace PhysBAM{

template<class TV> class EULER_UNIFORM;
template<class TV> class EULER_LAPLACE;
template<class TV> class POISSON_COLLIDABLE_UNIFORM;
template<class TV> class INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS;

template<class TV>
class EULER_PROJECTION_UNIFORM:public EULER_PROJECTION<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef typename ARRAY<T,FACE_INDEX<TV::m> >::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef typename ARRAY<T,TV_INT>::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    typedef CELL_ITERATOR<VECTOR<T,TV::m-1> > CELL_ITERATOR_LOWER_DIM;
    typedef CELL_ITERATOR<VECTOR<T,1> > CELL_ITERATOR_1D;
    
public:

    ARRAY<T,TV_INT> p; // p should always store the real pressure.
    ARRAY<T,TV_INT> p_save_for_projection;
    ARRAY<T,TV_INT> p_advected;
    ARRAY<T,TV_INT> one_over_rho_c_squared;
    bool is_pressure_scaled;
    T dt_scale_pressure;
    ARRAY<T,TV_INT> p_dirichlet;
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities,face_velocities_save;
    T_FACE_ARRAYS_DIMENSION_SCALAR *fluxes;
    EULER_UNIFORM<TV>* euler;
    EULER_LAPLACE<POISSON_COLLIDABLE_UNIFORM<TV> >* elliptic_solver;
    POISSON_COLLIDABLE_UNIFORM<TV>* poisson;
    ARRAY<T,TV_INT> divergence; // use this to set up a non-zero divergence
    BOUNDARY<TV,T>* pressure_boundary;
    BOUNDARY_REFLECTION_UNIFORM<TV,T> pressure_boundary_default;
    bool save_fluxes,use_exact_neumann_face_location,use_neumann_condition_for_outflow_boundaries;
#if 1
    ADVECTION_HAMILTON_JACOBI_ENO<VECTOR<T,1>,T> pressure_advection_HJ;
#else
    ADVECTION_HAMILTON_JACOBI_ENO<TV,T> pressure_advection_HJ;
#endif
    BOUNDARY_OBJECT_REFLECTION<TV,T> pressure_object_boundary;
    int hj_eno_order;

private:
    INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS<TV>* incompressible_coupling_callbacks;
    bool transition_to_using_implicit_pressure;
    
public:

    EULER_PROJECTION_UNIFORM(EULER_UNIFORM<TV>* euler_input);
    virtual ~EULER_PROJECTION_UNIFORM();

    void Save_State(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_s)
    {face_velocities_s.Copy(face_velocities);}

    void Restore_State(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_s)
    {face_velocities.Copy(face_velocities_s);}

    void Scale_Pressure_By_Dt(const T dt)
    {assert(!is_pressure_scaled);p*=dt;dt_scale_pressure=dt;is_pressure_scaled=true;}

    void Unscale_Pressure_By_Dt(const T dt)
    {assert(is_pressure_scaled && (dt==dt_scale_pressure));p*=(1/dt);is_pressure_scaled=false;}

    void Set_Transition_To_Using_Implicit_Pressure(const bool transition_to_using_implicit_pressure_input)
    {transition_to_using_implicit_pressure=transition_to_using_implicit_pressure_input;}

    void Exchange_Pressures_For_Projection()
    {ARRAY<T,TV_INT>::Exchange(p,p_save_for_projection);}

    void Set_Incompressible_Coupling_Callbacks(INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS<TV>* incompressible_coupling_callbacks_input)
    {incompressible_coupling_callbacks=incompressible_coupling_callbacks_input;}

    void Set_Use_Exact_Neumann_Face_Location(const bool use_exact_neumann_face_location_input)
    {use_exact_neumann_face_location=use_exact_neumann_face_location_input;
    euler->conservation->Set_Use_Exact_Neumann_Face_Location(use_exact_neumann_face_location);}

    void Set_Constant_Extrapolation_Pressure_Boundary()
    {pressure_boundary->Set_Constant_Extrapolation(euler->boundary->constant_extrapolation);}

    void Set_Custom_Pressure_Boundary(BOUNDARY<TV,T> &pressure_boundary_input)
    {pressure_boundary=&pressure_boundary_input;Set_Constant_Extrapolation_Pressure_Boundary();}

//#####################################################################
private:
    void Compute_Pressure(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    template<class FACE_LOOKUP> void Compute_Divergence(const FACE_LOOKUP &face_lookup);
public:
    virtual void Initialize_Grid();
    void Get_Pressure(ARRAY<T,TV_INT>& pressure) const;
    void Fill_Face_Weights_For_Projection(const T dt,const T time,ARRAY<T,FACE_INDEX<TV::m> >& beta_face);
    void Get_Ghost_Density(const T dt,const T time,const int number_of_ghost_cells,ARRAY<T,TV_INT>& density_ghost) const;
    void Get_Ghost_Centered_Velocity(const T dt,const T time,const int number_of_ghost_cells,ARRAY<TV,TV_INT>& centered_velocity_ghost) const;
    void Make_Boundary_Faces_Neumann(ARRAY<bool,FACE_INDEX<TV::m> >& psi_N);
    void Project(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    void Get_Dirichlet_Boundary_Conditions(const T_ARRAYS_DIMENSION_SCALAR& U_dirichlet);
    void Set_Dirichlet_Boundary_Conditions(const T time);
    void Compute_Advected_Pressure(const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<T,FACE_INDEX<TV::m> > face_velocities_for_solid_faces,const T dt);
    void Compute_Right_Hand_Side(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    void Compute_One_Over_rho_c_Squared();
    void Compute_Density_Weighted_Face_Velocities(const T dt,const T time,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N);
    static void Compute_Density_Weighted_Face_Velocities(const GRID<TV>& face_grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N);
    static void Compute_Face_Pressure_From_Cell_Pressures(const GRID<TV>& face_grid,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,ARRAY<T,FACE_INDEX<TV::m> >& p_face,const ARRAY<T,TV_INT>& p_cell);
    void Get_Ghost_Pressures(const T dt,const T time,const ARRAY<bool,TV_INT>& psi_D,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
        const ARRAY<T,TV_INT>& pressure,ARRAY<T,TV_INT>& p_ghost);
    void Get_Pressure_At_Faces(const T dt,const T time,const ARRAY<T,TV_INT>& p_ghost,ARRAY<T,FACE_INDEX<TV::m> >& p_face);
    void Apply_Pressure(const ARRAY<T,TV_INT>& p_ghost,const ARRAY<T,FACE_INDEX<TV::m> >& p_face,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_star,
        const ARRAY<bool,TV_INT>& psi_D,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T dt,const T time);
    static void Apply_Pressure(const ARRAY<T,TV_INT>& p_ghost,const ARRAY<T,FACE_INDEX<TV::m> >& p_face,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_star,
            const ARRAY<bool,TV_INT>& psi_D,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T dt,const T time,
            T_FACE_ARRAYS_DIMENSION_SCALAR *fluxes,EULER_UNIFORM<TV>* euler);
    bool Consistent_Boundary_Conditions() const;
    void Log_Parameters() const;
//#####################################################################
};
}
#endif
