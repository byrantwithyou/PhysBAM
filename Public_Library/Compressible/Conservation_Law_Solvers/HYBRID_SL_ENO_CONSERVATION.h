//#####################################################################
// Copyright 2010-2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYBRID_SL_ENO_CONSERVATION
//##################################################################### 
#ifndef __HYBRID_SL_ENO_CONSERVATION__
#define __HYBRID_SL_ENO_CONSERVATION__   

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class EULER_EIGENSYSTEM;
template<class T,class TV_DIMENSION> class EIGENSYSTEM;
template<class TV> class GRID;

template<class TV,int d>
class HYBRID_SL_ENO_CONSERVATION:public CONSERVATION<TV,d>
{
    typedef typename TV::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
    typedef VECTOR<int,TV::m> TV_INT;typedef TV_INT INDEX;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef ARRAY<TV_DIMENSION,FACE_INDEX<TV::m> > T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef CONSERVATION<TV,d> BASE;

    const ARRAY<bool,FACE_INDEX<TV::m> >& flux_face;
    CONSERVATION<TV,d> *conservation;

public:
    HYBRID_SL_ENO_CONSERVATION(const ARRAY<bool,FACE_INDEX<TV::m> >& flux_face_in, CONSERVATION<TV,d> *conservation_law_solver)
        : flux_face(flux_face_in),conservation(conservation_law_solver)
    {conservation->Save_Fluxes();}

    void Set_Order(const int order_input=3)
    {conservation->Set_Order(order_input);}
            
    void Use_Field_By_Field_Alpha()
    {conservation->Use_Field_By_Field_Alpha();}

    void Use_Maximum_Alpha()
    {conservation->Use_Maximum_Alpha();}

    void Set_Use_Exact_Neumann_Face_Location(const bool use_exact_neumann_face_location_input)
    {conservation->Set_Use_Exact_Neumann_Face_Location(use_exact_neumann_face_location_input);}

    void Amplify_Alpha(const T amplification_factor_input=1)
    {conservation->Amplify_Alpha(amplification_factor_input);}

    void Save_Fluxes()
    {conservation->Save_Fluxes();}

    void Scale_Outgoing_Fluxes_To_Clamp_Variable(bool scale_outgoing_fluxes_to_clamp_variable_input,int clamped_variable_index_input,T clamped_value_input)
    {conservation->Scale_Outgoing_Fluxes_To_Clamp_Variable(scale_outgoing_fluxes_to_clamp_variable_input, clamped_variable_index_input, clamped_value_input);}

    void Set_Callbacks(CONSERVATION_CALLBACKS<T> *callbacks_input)
    {conservation->Set_Callbacks(callbacks_input);}

    void Set_Custom_Object_Boundary(BOUNDARY_OBJECT<TV,TV_DIMENSION>& object_boundary_input)
    {conservation->Set_Custom_Object_Boundary(object_boundary_input);}

    virtual void Log_Parameters() const
    {BASE::Log_Parameters();conservation->Log_Parameters();}

    virtual void Update_Conservation_Law(GRID<TV>& grid,T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
        VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
        const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool thinshell=false,const VECTOR<bool,2*TV::m>& outflow_boundaries=(VECTOR<bool,2*TV::m>()),VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,TV::m>* eigensystems_auxiliary=0,
        T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary=0);
};
}
#endif
