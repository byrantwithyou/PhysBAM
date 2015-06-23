//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DENSITY_CONTAINER
//#####################################################################
#ifndef __DENSITY_CONTAINER__    
#define __DENSITY_CONTAINER__

#include <Incompressible/Advection_Collidable/ADVECTION_COLLIDABLE_FORWARD.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <Incompressible/Grid_Based_Fields/GRID_AND_ARRAY_CONTAINER.h>
namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV>
class DENSITY_CONTAINER:public GRID_AND_ARRAY_CONTAINER<TV,typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM<TV,T> T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> T_FACE_LOOKUP_COLLIDABLE;
public:
    typedef GRID_AND_ARRAY_CONTAINER<TV,T> BASE;
    using BASE::array;using BASE::grid;
    using BASE::boundary_default;using BASE::boundary;using BASE::Set_To_Constant_Value;using BASE::Set_Custom_Advection;

    ARRAY<T,TV_INT>& density;
    T ambient_density;

private:
    ARRAY<bool,TV_INT> valid_mask_current;
    ARRAY<bool,TV_INT> valid_mask_next;
public:
    T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL* nested_semi_lagrangian_collidable;
private:
    ADVECTION_WRAPPER_COLLIDABLE_CELL<TV,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>* semi_lagrangian_collidable;
public:    

    DENSITY_CONTAINER(GRID<TV>& grid_input);
    ~DENSITY_CONTAINER();

    void Set_Ambient_Density(const T density_input=0)
    {ambient_density=density_input;}

    void Set_To_Ambient_Density()
    {Set_To_Constant_Value(ambient_density);}

//#####################################################################
    void Euler_Step(const T dt,const T time,const int number_of_ghost_cells) override;
    void Initialize_Array(const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true) override;
    void Use_Semi_Lagrangian_Collidable_Advection(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list,const ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid_mask_input);
    void Fill_Beta_At_Faces(const T dt,const T time,ARRAY<T,FACE_INDEX<TV::m> >& beta_face) const;
    void Get_Ghost_Density(const T dt,const T time,const int number_of_ghost_cells,ARRAY<T,TV_INT>& density_ghost) const;
//#####################################################################
};      
}
#endif
