//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_BOUNDARY_CONDITION_COLLECTION
//#####################################################################
#ifndef __IMPLICIT_BOUNDARY_CONDITION_COLLECTION__
#define __IMPLICIT_BOUNDARY_CONDITION_COLLECTION__
#include <Core/Arrays/ARRAY.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Incompressible/Solids_And_Fluids/BOUNDARY_CONDITIONS_CALLBACKS.h>
#include <Dynamics/Coupled_Evolution/BOUNDARY_CONDITION_INFO.h>
#include <Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{
template<class TV>
class IMPLICIT_BOUNDARY_CONDITION_COLLECTION
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
public:
    ARRAY<bool,TV_INT> psi_D;
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N;
    ARRAY<T,FACE_INDEX<TV::m> > psi_R;

    BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback;
    ARRAY<BOUNDARY_CONDITION_INFO<TV> > boundary_condition_info;

    ARRAY<IMPLICIT_BOUNDARY_CONDITION<TV>*> boundary_conditions;

    bool set_all_neumann_cells_to_dirichlet,zero_all_dirichlet_face_velocities,use_psi_R;
    bool use_boundary_condition_info;
    VECTOR<bool,TV::m> periodic_boundary;

    IMPLICIT_BOUNDARY_CONDITION_COLLECTION(BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input,bool set_all_neumann_cells_to_dirichlet_input,
        bool zero_all_dirichlet_face_velocities_input,bool use_psi_R_input,bool use_boundary_condition_info_input,VECTOR<bool,TV::m> periodic_boundary_input);
    ~IMPLICIT_BOUNDARY_CONDITION_COLLECTION();

    void Add_Boundary_Condition(IMPLICIT_BOUNDARY_CONDITION<TV>* boundary_condition)
    {boundary_conditions.Append(boundary_condition);}

    bool All_Cell_Faces_Neumann(const TV_INT& cell_index) const;
    void Compute(const GRID<TV>& grid,ARRAY<T,TV_INT>& p,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time);
    void Compute_Boundary_Condition_Info(const GRID<TV>& grid,const ARRAY<T,TV_INT>& p,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Compute_Boundary_Condition_Info(const ARRAY<T,TV_INT>& p,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const FACE_INDEX<TV::m>& f,int in_side);
};
}
#endif
