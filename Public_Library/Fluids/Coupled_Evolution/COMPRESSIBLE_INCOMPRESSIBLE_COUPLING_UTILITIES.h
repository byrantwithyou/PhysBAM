//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES
//#####################################################################
#ifndef __COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES__
#define __COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES__
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Compressible/Euler_Equations/INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS.h>
namespace PhysBAM{
template<class T> class EOS;
template<class TV> class EULER;

template<class TV>
class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES:public INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<VECTOR<T,TV::m+2>,TV_INT> T_ARRAYS_DIMENSION_SCALAR;

private:
    const ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities_;
    const T& incompressible_density_;
    const ARRAY<T,TV_INT>& incompressible_phi_;
public:
    ARRAY<T,TV_INT> p_dirichlet_incompressible;

    COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >* incompressible_face_velocities,const T* incompressible_density,
        const ARRAY<T,TV_INT>* incompressible_phi);

//#####################################################################
    static void Extrapolate_Compressible_State_Into_Incompressible_Region(const T dt,const T time,const T bandwidth,const int ghost_cells,const EOS<T>& eos,const GRID<TV>& grid,
        const ARRAY<T,TV_INT>& phi_ghost,const ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,T_ARRAYS_DIMENSION_SCALAR& U);
    void Get_Dirichlet_Boundary_Conditions_For_Incompressible_Region(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_SCALAR& U_dirichlet,const EOS<T>& euler_eos,const T incompressible_density,const T dt);
    static void Compute_Compressible_Incompressible_Face_Velocities(const GRID<TV>& face_grid,const ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities,const T incompressible_density,
        const ARRAY<T,TV_INT>& incompressible_phi,const T_ARRAYS_DIMENSION_SCALAR& U,const ARRAY<bool,TV_INT>& euler_psi,ARRAY<T,FACE_INDEX<TV::m> >& compressible_face_velocities);
    static void Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(const GRID<TV>& face_grid,
        const ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities,const T incompressible_density,const ARRAY<T,TV_INT>& incompressible_phi,
        const T_ARRAYS_DIMENSION_SCALAR& U,const ARRAY<bool,TV_INT>& euler_psi,const ARRAY<T,TV_INT>& p_cell,ARRAY<T,FACE_INDEX<TV::m> >& p_face);
    virtual void Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(const GRID<TV>& face_grid,const T_ARRAYS_DIMENSION_SCALAR& U,const ARRAY<bool,TV_INT>& euler_psi,const ARRAY<T,TV_INT>& p_cell,ARRAY<T,FACE_INDEX<TV::m> >& p_face) const PHYSBAM_OVERRIDE;
    static void Fill_Incompressible_Beta_Face(const GRID<TV>& grid,const T incompressible_density,const ARRAY<T,TV_INT>& incompressible_phi,
        ARRAY<T,FACE_INDEX<TV::m> >& beta_face);
    virtual void Fill_Incompressible_Beta_Face(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& beta_face) const PHYSBAM_OVERRIDE;
    static void Apply_Pressure_At_Incompressible_Faces(const GRID<TV>& face_grid,const T incompressible_density,
        const ARRAY<T,TV_INT>& incompressible_phi,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const ARRAY<T,TV_INT>& p_hat,
        ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities);
//#####################################################################

};
}
#endif
