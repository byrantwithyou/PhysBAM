//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/PROJECTION_COLLIDABLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PROJECTION_COLLIDABLE_UNIFORM<TV>::
PROJECTION_COLLIDABLE_UNIFORM(const GRID<TV>& mac_grid,const bool multiphase,const bool use_poisson,const bool use_variable_beta,THREAD_QUEUE* thread_queue)
    :laplace_collidable(0),poisson_collidable(0)
{
    if(use_poisson){
        poisson=poisson_collidable=new POISSON_COLLIDABLE_UNIFORM<TV>(p_grid,p,true,multiphase,true);
        if(use_variable_beta) poisson->Set_Variable_beta();elliptic_solver=poisson;collidable_solver=poisson_collidable;}
    else{
        laplace=laplace_collidable=new LAPLACE_COLLIDABLE_UNIFORM<TV>(p_grid,p,true,false,true,thread_queue);elliptic_solver=laplace;collidable_solver=laplace_collidable;}
    Initialize_Grid(mac_grid);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PROJECTION_COLLIDABLE_UNIFORM<TV>::
PROJECTION_COLLIDABLE_UNIFORM(const GRID<TV>& mac_grid,LEVELSET<TV>& levelset_input)
    :laplace_collidable(0),poisson_collidable(0)
{
    poisson=poisson_collidable=new POISSON_COLLIDABLE_UNIFORM<TV>(p_grid,p,levelset_input,true,false,true);elliptic_solver=poisson;collidable_solver=poisson_collidable;
    Initialize_Grid(mac_grid);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PROJECTION_COLLIDABLE_UNIFORM<TV>::
~PROJECTION_COLLIDABLE_UNIFORM()
{
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class TV> void PROJECTION_COLLIDABLE_UNIFORM<TV>::
Initialize_Grid(const GRID<TV>& mac_grid)
{
    BASE::Initialize_Grid(mac_grid);
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class TV> void PROJECTION_COLLIDABLE_UNIFORM<TV>::
Apply_Pressure(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,bool scale_by_dt)
{
    // find divergence free u, v and w 
    if(collidable_solver->second_order_cut_cell_method){
        //Zero_Out_Neumann_Pocket_Velocities(face_velocities); TODO: Why is this here?
        ARRAY<bool,TV_INT>& psi_D=elliptic_solver->psi_D;
        ARRAY<bool,FACE_INDEX<TV::m> >& psi_N=elliptic_solver->psi_N;
        ARRAY<T,TV_INT>& phi=collidable_solver->levelset->phi;
        TV dx=p_grid.dX,one_over_dx=Inverse(dx);
        if(scale_by_dt) p*=dt;
        if(laplace){
            for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(!psi_N.Component(axis)(face_index) && !(psi_D(first_cell) && psi_D(second_cell))){
                    if(psi_D(first_cell) && !psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(second_cell),phi(first_cell)))
                        face_velocities.Component(axis)(face_index)-=(p(second_cell)-collidable_solver->u_interface.Component(axis)(face_index))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(second_cell),phi(first_cell)),collidable_solver->second_order_cut_cell_threshold)*dx[axis]);
                    else if(!psi_D(first_cell) && psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(first_cell),phi(second_cell)))
                        face_velocities.Component(axis)(face_index)-=(collidable_solver->u_interface.Component(axis)(face_index)-p(first_cell))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(first_cell),phi(second_cell)),collidable_solver->second_order_cut_cell_threshold)*dx[axis]);
                    else face_velocities.Component(axis)(face_index)-=(p(second_cell)-p(first_cell))*one_over_dx[axis];}}}
        else if(poisson){
            assert(!poisson->use_variable_beta); // assumes constant beta in each phase
            for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(!psi_N.Component(axis)(face_index) && !(psi_D(first_cell) && psi_D(second_cell))){
                    if(psi_D(first_cell) && !psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(second_cell),phi(first_cell)))
                        face_velocities.Component(axis)(face_index)-=poisson->beta_face.Component(axis)(face_index)*(p(second_cell)-collidable_solver->u_interface.Component(axis)(face_index))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(second_cell),phi(first_cell)),collidable_solver->second_order_cut_cell_threshold)*dx[axis]);
                    else if(!psi_D(first_cell) && psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(first_cell),phi(second_cell)))
                        face_velocities.Component(axis)(face_index)-=poisson->beta_face.Component(axis)(face_index)*(collidable_solver->u_interface.Component(axis)(face_index)-p(first_cell))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(first_cell),phi(second_cell)),collidable_solver->second_order_cut_cell_threshold)*dx[axis]);
                    else face_velocities.Component(axis)(face_index)-=poisson->beta_face.Component(axis)(face_index)*(p(second_cell)-p(first_cell))*one_over_dx[axis];}}}
            if(scale_by_dt) p*=1/dt;}
    else
        BASE::Apply_Pressure(face_velocities,dt,time,scale_by_dt);
}
//#####################################################################
namespace PhysBAM{
template class PROJECTION_COLLIDABLE_UNIFORM<VECTOR<float,1> >;
template class PROJECTION_COLLIDABLE_UNIFORM<VECTOR<float,2> >;
template class PROJECTION_COLLIDABLE_UNIFORM<VECTOR<float,3> >;
template class PROJECTION_COLLIDABLE_UNIFORM<VECTOR<double,1> >;
template class PROJECTION_COLLIDABLE_UNIFORM<VECTOR<double,2> >;
template class PROJECTION_COLLIDABLE_UNIFORM<VECTOR<double,3> >;
}
