//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_VISCOSITY_UNIFORM
//#####################################################################
#ifndef __LEVELSET_VISCOSITY_UNIFORM__
#define __LEVELSET_VISCOSITY_UNIFORM__

#include <Incompressible/Incompressible_Flows/LEVELSET_INDEX_MAP_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM_SYSTEM.h>
namespace PhysBAM{

template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;

template<class TV>
class LEVELSET_VISCOSITY_UNIFORM
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
public:

    LEVELSET_INDEX_MAP_UNIFORM<TV> index_map;
    LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV> system;
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > x,b;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    T scale;
    bool print_matrix;
    VECTOR<bool,d> periodic_boundary;

    LEVELSET_VISCOSITY_UNIFORM(BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input,const GRID<TV>& grid_input,T dt,T density,T viscosity);
    LEVELSET_VISCOSITY_UNIFORM(const LEVELSET_VISCOSITY_UNIFORM&) = delete;
    void operator=(const LEVELSET_VISCOSITY_UNIFORM&) = delete;
    virtual ~LEVELSET_VISCOSITY_UNIFORM();

//#####################################################################
    void Apply_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,bool fully_explicit,bool fully_implicit,bool coupled);
    void Apply_Full_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,bool fully_explicit,bool fully_implicit,int axis);
    void Apply_Implicit_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,int axis);
    void Apply_Explicit_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,int axis);
    void Resize_Vectors(bool minimal);
//#####################################################################
};
}
#endif
