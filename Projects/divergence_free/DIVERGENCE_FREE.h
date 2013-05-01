//#####################################################################
// Copyright 2013
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIVERGENCE_FREE
//#####################################################################
#ifndef __DIVERGENCE_FREE__
#define __DIVERGENCE_FREE__

#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class DIVERGENCE_FREE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    const GRID<TV>& grid; // normal grid
    const GRID<TV>& mac_grid; // corresponding mac grid
    ARRAY<TV,TV_INT>& v; // velocities on grid nodes
    const ARRAY<T,TV_INT>& m; // masses on grid nodes
    T m_eps; // the minimum mass to identify Dirichlet pressure cells

    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    PROJECTION_UNIFORM<GRID<TV> > projection;

    DIVERGENCE_FREE(const GRID<TV>& grid_in,const GRID<TV>& mac_grid_in,ARRAY<TV,TV_INT>& v_in,const ARRAY<T,TV_INT>& m_in,const T m_eps_in);
    ~DIVERGENCE_FREE();

    void Interpolate_Velocities_Nodes_To_Faces();
    void Identify_Dirichlet_Cells(); // vacuum
    void Identify_Neumann_Cells(); // collision objects
    void Set_Up_Solver(); // prameters for cg solver
    void Make_Divergence_Free(T dt,T time); // poisson solve and projection
};
}
#endif
