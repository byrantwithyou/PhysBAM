//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <Tools/Log/LOG.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <Incompressible/Interpolation_Collidable/AVERAGING_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VORTICITY_CONFINEMENT<TV>::
VORTICITY_CONFINEMENT(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_body_list,ARRAY<bool,FACE_INDEX<TV::m> >* valid_mask,const bool use_variable_vorticity_confinement,const T vorticity_confinement)
    :collision_body_list(collision_body_list),valid_mask(valid_mask),use_variable_vorticity_confinement(use_variable_vorticity_confinement),vorticity_confinement(vorticity_confinement)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> VORTICITY_CONFINEMENT<TV>::
~VORTICITY_CONFINEMENT()
{
}
//#####################################################################
// Function Apply_Vorticity_Confinement_Force
//#####################################################################
template<class TV,class T_ARRAYS_TV,class T> static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T_ARRAYS_TV& F,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_body_list)
{
    typedef VECTOR<int,TV::m> TV_INT;
    // want cells to face averaging here
    if(collision_body_list){
        AVERAGING_COLLIDABLE_UNIFORM<TV,FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> > vorticity_averaging_collidable(*collision_body_list,T());
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities.Component(axis)(iterator.Face_Index())+=vorticity_averaging_collidable.Cell_To_Face(grid,axis,iterator.Face_Index(),F);}}
    else
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
            face_velocities.Component(axis)(iterator.Face_Index())+=F(iterator.First_Cell_Index())[axis]+F(iterator.Second_Cell_Index())[axis];}
}
static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<float,1> >&,ARRAY<float,FACE_INDEX<1> >&,ARRAY<VECTOR<float,1> ,VECTOR<int,1> >&,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<float,1> >*)
{PHYSBAM_NOT_IMPLEMENTED();}
static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<double,1> >&,ARRAY<double,FACE_INDEX<1> >&,ARRAY<VECTOR<double,1> ,VECTOR<int,1> >&,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<double,1> >*)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void VORTICITY_CONFINEMENT<TV>::
Apply_Vorticity_Confinement_Force(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<TV,TV_INT>& F)
{
    Apply_Vorticity_Confinement_Force_Helper(grid,face_velocities,F,collision_body_list);
}
//#####################################################################
// Function Compute_Vorticity_Confinement_Force
//#####################################################################
template<class TV,class T,class TV_INT,class T_VALID_MASK> static void 
Compute_Vorticity_Confinement_Force_Helper(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<TV,TV_INT>& F,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_body_list,const T_VALID_MASK* valid_mask)
{
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    ARRAY<typename TV::SPIN,TV_INT> vorticity(grid.Cell_Indices(2),false);
    ARRAY<T,TV_INT> vorticity_magnitude(grid.Cell_Indices(2));
    if(collision_body_list){
        FACE_LOOKUP_UNIFORM<TV> face_velocities_lookup_uniform(face_velocities_ghost);
        FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> face_velocities_lookup(face_velocities_lookup_uniform,*collision_body_list,valid_mask);
        VORTICITY_UNIFORM<TV>::Vorticity(grid,face_velocities_lookup,vorticity,vorticity_magnitude);}
    else VORTICITY_UNIFORM<TV>::Vorticity(grid,T_FACE_LOOKUP(face_velocities_ghost),vorticity,vorticity_magnitude);
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){ // do collision awareness when these are averaged to faces
        TV vortex_normal_vector=LEVELSET<TV>::Normal_At_Node(grid,vorticity_magnitude,iterator.Cell_Index());
        F(iterator.Cell_Index())=TV::Cross_Product(vortex_normal_vector,vorticity(iterator.Cell_Index()));}
}
template<class T,class TV_INT,class T_VALID_MASK> static void
Compute_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<T,1> >&,const ARRAY<T,FACE_INDEX<1> >&,ARRAY<VECTOR<T,1>,TV_INT>&,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<VECTOR<T,1> >*,const T_VALID_MASK*)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class TV> void VORTICITY_CONFINEMENT<TV>::
Compute_Vorticity_Confinement_Force(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const ARRAY<bool,FACE_INDEX<TV::m> >* valid_mask,ARRAY<TV,TV_INT>& F)
{
    Compute_Vorticity_Confinement_Force_Helper(grid,face_velocities_ghost,F,collision_body_list,valid_mask);
}
//#####################################################################
// Function Add_Explicit_Forces
//#####################################################################
template<class TV> void VORTICITY_CONFINEMENT<TV>::
Add_Explicit_Forces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    ARRAY<TV,TV_INT> F(grid.Cell_Indices(1),false);
    Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,valid_mask,F);
    if(collision_body_list){
        if(use_variable_vorticity_confinement){F*=dt;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement;}
    else{
        if(use_variable_vorticity_confinement){F*=dt*(T).5;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement*(T).5;}
    Apply_Vorticity_Confinement_Force(grid,face_velocities,F);
}
//#####################################################################
// Function Add_Implicit_Forces_Projection
//#####################################################################
template<class TV> void VORTICITY_CONFINEMENT<TV>::
Add_Implicit_Forces_Projection(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class TV> void VORTICITY_CONFINEMENT<TV>::
Initialize_Grids(const GRID<TV>& grid)
{
    if(use_variable_vorticity_confinement) variable_vorticity_confinement.Resize(grid.Cell_Indices(1));
    else variable_vorticity_confinement.Clean_Memory();
}
//#####################################################################
namespace PhysBAM{
template class VORTICITY_CONFINEMENT<VECTOR<float,1> >;
template class VORTICITY_CONFINEMENT<VECTOR<float,2> >;
template class VORTICITY_CONFINEMENT<VECTOR<float,3> >;
template class VORTICITY_CONFINEMENT<VECTOR<double,1> >;
template class VORTICITY_CONFINEMENT<VECTOR<double,2> >;
template class VORTICITY_CONFINEMENT<VECTOR<double,3> >;
}
