//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES
//#####################################################################
#ifndef __COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES__
#define __COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES__
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS.h>
namespace PhysBAM{
template<class T> class EOS;
template<class T_GRID> class EULER;

template<class TV>
class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES:public INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS<TV>
{
    typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<VECTOR<T,T_GRID::dimension+2> >::TYPE T_ARRAYS_DIMENSION_SCALAR;

private:
    const T_FACE_ARRAYS_SCALAR& incompressible_face_velocities_;
    const T& incompressible_density_;
    const T_ARRAYS_SCALAR& incompressible_phi_;
public:
    T_ARRAYS_SCALAR p_dirichlet_incompressible;

    COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR* incompressible_face_velocities,const T* incompressible_density,
        const T_ARRAYS_SCALAR* incompressible_phi);

//#####################################################################
    static void Extrapolate_Compressible_State_Into_Incompressible_Region(const T dt,const T time,const T bandwidth,const int ghost_cells,const EOS<T>& eos,const T_GRID& grid,
        const T_ARRAYS_SCALAR& phi_ghost,const T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,T_ARRAYS_DIMENSION_SCALAR& U);
    void Get_Dirichlet_Boundary_Conditions_For_Incompressible_Region(const T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& U_dirichlet,const EOS<T>& euler_eos,const T incompressible_density,const T dt);
    static void Compute_Compressible_Incompressible_Face_Velocities(const T_GRID& face_grid,const T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T incompressible_density,
        const T_ARRAYS_SCALAR& incompressible_phi,const T_ARRAYS_DIMENSION_SCALAR& U,const ARRAY<bool,TV_INT>& euler_psi,T_FACE_ARRAYS_SCALAR& compressible_face_velocities);
    static void Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(const T_GRID& face_grid,
        const T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T incompressible_density,const T_ARRAYS_SCALAR& incompressible_phi,
        const T_ARRAYS_DIMENSION_SCALAR& U,const ARRAY<bool,TV_INT>& euler_psi,const T_ARRAYS_SCALAR& p_cell,T_FACE_ARRAYS_SCALAR& p_face);
    virtual void Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(const T_GRID& face_grid,const T_ARRAYS_DIMENSION_SCALAR& U,const ARRAY<bool,TV_INT>& euler_psi,const T_ARRAYS_SCALAR& p_cell,T_FACE_ARRAYS_SCALAR& p_face) const PHYSBAM_OVERRIDE;
    static void Fill_Incompressible_Beta_Face(const T_GRID& grid,const T incompressible_density,const T_ARRAYS_SCALAR& incompressible_phi,
        T_FACE_ARRAYS_SCALAR& beta_face);
    virtual void Fill_Incompressible_Beta_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& beta_face) const PHYSBAM_OVERRIDE;
    static void Apply_Pressure_At_Incompressible_Faces(const T_GRID& face_grid,const T incompressible_density,
        const T_ARRAYS_SCALAR& incompressible_phi,const T_FACE_ARRAYS_BOOL& psi_N,const T_ARRAYS_SCALAR& p_hat,
        T_FACE_ARRAYS_SCALAR& incompressible_face_velocities);
//#####################################################################

};
}
#endif
