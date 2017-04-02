//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION  
//##################################################################### 
//
// This class solves U_t + F(U)_x + G(U)_y + H(U)_z = 0 for one Euler step.
//
// Input vector U as (0,m) by (0,n) by (0,mn).
// Input vector U_ghost as (-2,m+3) by (-2,n+3) by (-2,mn+3). 
// Input psi as (0,m) by (0,n) by (0,mn). When psi=true, solve the equaitions. When psi=false, do NOT solve the equations.
//
//#####################################################################
#ifndef __CONSERVATION__
#define __CONSERVATION__   

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Utilities/STATIC_ASSERT.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Arrays/FLOOD_FILL.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class EULER_EIGENSYSTEM;
template<class T,int d> class EIGENSYSTEM;
template<class TV,class TV_DIMENSION> class BOUNDARY_OBJECT;
template<class TV,class TV_DIMENSION> class BOUNDARY_OBJECT_REFLECTION;
template<class TV> class GRID;

template<class TV,int d>
class CONSERVATION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<T,d> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef ARRAY<TV_DIMENSION,FACE_INDEX<TV::m> > T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef VECTOR<T,TV::m-1> TV_LOWER_DIM;
public: 
    int order; // 1,2, or 3
    int field_by_field_alpha; // (1) field by field, (0) max over all the fields 
    T amplification_factor; // for amplifying alpha to increase dissipation
    int save_fluxes; // saves the fluxes
    ARRAY<TV_DIMENSION,VECTOR<int,1> > flux_temp; // temporary array used for saving 1d fluxes 
    CONSERVATION_CALLBACKS<T> *callbacks;
    BOUNDARY_OBJECT<TV,TV_DIMENSION> *object_boundary;
    BOUNDARY_OBJECT_REFLECTION<TV,TV_DIMENSION>& object_boundary_default;
    T_FACE_ARRAYS_DIMENSION_SCALAR fluxes;
    bool use_exact_neumann_face_location,scale_outgoing_fluxes_to_clamp_variable;
    int clamped_variable_index;
    T clamped_value;
    // temp arrays used inside functions declared here for optimization
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_flux_1d_axis;
public: 

    CONSERVATION();

    virtual void Set_Order(const int order_input=3)
    {order=order_input;}
    
    virtual void Use_Field_By_Field_Alpha()
    {field_by_field_alpha=1;}

    virtual void Use_Maximum_Alpha()
    {field_by_field_alpha=0;}

    virtual void Set_Use_Exact_Neumann_Face_Location(const bool use_exact_neumann_face_location_input)
    {use_exact_neumann_face_location=use_exact_neumann_face_location_input;}

    virtual void Amplify_Alpha(const T amplification_factor_input=1)
    {amplification_factor=amplification_factor_input;}

    virtual void Save_Fluxes()
    {save_fluxes=1;}

    virtual void Scale_Outgoing_Fluxes_To_Clamp_Variable(bool scale_outgoing_fluxes_to_clamp_variable_input,int clamped_variable_index_input,T clamped_value_input)
    {scale_outgoing_fluxes_to_clamp_variable=scale_outgoing_fluxes_to_clamp_variable_input;clamped_variable_index=clamped_variable_index_input;clamped_value=clamped_value_input;
    save_fluxes=save_fluxes||scale_outgoing_fluxes_to_clamp_variable;}

    virtual void Set_Callbacks(CONSERVATION_CALLBACKS<T> *callbacks_input)
    {callbacks=callbacks_input;}

    virtual void Set_Custom_Object_Boundary(BOUNDARY_OBJECT<TV,TV_DIMENSION>& object_boundary_input)
    {object_boundary=&object_boundary_input;}

//#####################################################################
    virtual ~CONSERVATION();
    void Compute_Flux_Without_Clamping(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
        VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
        const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const VECTOR<bool,2*TV::m>& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,const T_ARRAYS_DIMENSION_SCALAR* U_ghost_clamped=0);
    void Compute_Flux_With_Clamping(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
        VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
        const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const VECTOR<bool,2*TV::m>& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,
        VECTOR<EIGENSYSTEM<T,d>*,TV::m>* eigensystems_auxiliary,T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary);
    void Compute_Flux(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
        VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
        const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const VECTOR<bool,2*TV::m>& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,
        VECTOR<EIGENSYSTEM<T,d>*,TV::m>* eigensystems_auxiliary,T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary);
    virtual void Update_Conservation_Law(GRID<TV>& grid,T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
        VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
        const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool thinshell=false,const VECTOR<bool,2*TV::m>& outflow_boundaries=(VECTOR<bool,2*TV::m>()),VECTOR<EIGENSYSTEM<T,d>*,TV::m>* eigensystems_auxiliary=0,
        T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary=0);
    template<class T_ARRAYS> void Update_Conservation_Law_For_Specialized_Shallow_Water_Equations(GRID<TV>& grid,T_ARRAYS& U,const T_ARRAYS& U_ghost,const ARRAY<bool,VECTOR<int,2> >& psi,const T dt,
        EIGENSYSTEM<T,2>& eigensystem_F,EIGENSYSTEM<T,2>& eigensystem_G,CONSERVATION<TV,2>& solver,const VECTOR<bool,2*TV::m>& outflow_boundaries=(VECTOR<bool,2*TV::m>()));
    T Alpha(const VECTOR<T,d>& lambda_left,const VECTOR<T,d>& lambda_right,const int k,const int length);
    void Compute_Delta_Flux_For_Clamping_Variable(const GRID<TV>& grid,const int number_of_ghost_cells,T dt,const int clamped_variable_index,const T clamped_value,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
        const T_ARRAYS_DIMENSION_SCALAR& U,const T_FACE_ARRAYS_DIMENSION_SCALAR& flux,T_FACE_ARRAYS_DIMENSION_SCALAR& delta_flux,T_ARRAYS_DIMENSION_SCALAR& rhs,ARRAY<T,TV_INT>& overshoot_percentages);
    virtual void Log_Parameters() const;
    virtual void Conservation_Solver(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,d>& eigensystem,EIGENSYSTEM<T,d>& eigensystem_explicit,
        const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux=0) {PHYSBAM_FATAL_ERROR();}
//#####################################################################
};
}
#endif
