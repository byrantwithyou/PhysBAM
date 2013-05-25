//#####################################################################
// Copyright 2013
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include "MPM_POISSON_VECTOR.h"
#include "MPM_POISSON_SYSTEM.h"
#include "MPM_PROJECTION.h"
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PROJECTION<TV>::
MPM_PROJECTION(MPM_SIMULATION<TV>& sim_in)
    :sim(sim_in)
{
    mac_grid.Initialize(sim.grid.numbers_of_cells,RANGE<TV>(sim.grid.domain.min_corner,sim.grid.domain.max_corner),true);
    face_velocities.Resize(mac_grid);face_masses.Resize(mac_grid);face_momenta.Resize(mac_grid);
    cell_dirichlet.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    cell_neumann.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    div_u.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    pressure.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PROJECTION<TV>::
~MPM_PROJECTION()
{}

//#####################################################################
// Function Reinitialize
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Reinitialize()
{
    face_velocities.Fill((T)0);
    face_masses.Fill((T)0);
    face_momenta.Fill((T)0);
    cell_dirichlet.Fill(false);
    cell_neumann.Fill(false);
    nodes_non_dirichlet_cells.Clean_Memory();
    div_u.Fill((T)0);
    pressure.Fill((T)0);
}

//#####################################################################
// Function Identify_Dirichlet_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Identify_Dirichlet_Cells()
{
    cell_dirichlet.Fill(true);
#pragma omp parallel for
    for(int p=0;p<sim.particles.number;p++){
        PHYSBAM_ASSERT(mac_grid.Cell(sim.particles.X(p),0)==sim.grid.Cell(sim.particles.X(p),0));
        TV_INT cell=mac_grid.Cell(sim.particles.X(p),0);
        cell_dirichlet(cell)=false;}
}

//#####################################################################
// Function Identify_Neumann_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Identify_Neumann_Cells()
{
    // HMandatory wall : the third and fourth outmost MAC grid layer
    cell_neumann.Fill(false);
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),mac_grid.counts));it.Valid();it.Next()){
        for(int d=0;d<TV::m;d++){
            if(it.index(d)==2 || it.index(d)==3 || it.index(d)==mac_grid.counts(d)-3 || it.index(d)==mac_grid.counts(d)-4){
                cell_neumann(it.index)=true;
                if(cell_dirichlet(it.index)) cell_dirichlet(it.index)=false;}}}
}

//#####################################################################
// Function Identify_Nodes_Of_Non_Dirichlet_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Identify_Nodes_Of_Non_Dirichlet_Cells()
{
    nodes_non_dirichlet_cells.Clean_Memory();
    TV_INT nodes[TV::m*4-4];
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),mac_grid.counts));it.Valid();it.Next()){
        if(!cell_dirichlet(it.index)){
            mac_grid.Nodes_In_Cell_From_Minimum_Corner_Node(it.index,nodes);
            for(int i=0;i<TV::m*4-4;i++)
                nodes_non_dirichlet_cells.Get_Or_Insert(nodes[i])=true;}}
}

//#####################################################################
// Function Velocities_Corners_To_Faces_MPM_Style
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Velocities_Corners_To_Faces_MPM_Style()
{
    face_velocities.Fill((T)0);
    face_masses.Fill((T)0);
    face_momenta.Fill((T)0);
    ARRAY<FACE_INDEX<TV::m> > faces_got_rasterized;
    for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(nodes_non_dirichlet_cells);it.Valid();it.Next()){
        TV_INT node=it.Key();
        for(int axis=0;axis<TV::m;axis++){
            for(int face=0;face<TV::m*2-2;face++){
                FACE_INDEX<TV::m> face_index(axis,mac_grid.Node_Face_Index(axis,node,face));
                faces_got_rasterized.Append(face_index);
                face_masses(face_index)+=0.5*sim.node_mass(node);
                face_momenta(face_index)+=0.5*sim.node_mass(node)*sim.node_V(node)(axis);}}}
    for(int i=0;i<faces_got_rasterized.m;i++){
        FACE_INDEX<TV::m> face_index=faces_got_rasterized(i);
        if(face_masses(face_index)>sim.min_mass)
            face_velocities(face_index)=face_momenta(face_index)/face_masses(face_index);}
}

//#####################################################################
// Function Build_Velocity_Divergence
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Build_Velocity_Divergence()
{
    div_u.Fill((T)0);
    T one_over_h=(T)1/mac_grid.dX.Min();
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),mac_grid.counts));it.Valid();it.Next()){
        if(!cell_dirichlet(it.index)){
            for(int axis=0;axis<TV::m;axis++)
                div_u(it.index)+=face_velocities(FACE_INDEX<TV::m>(axis,mac_grid.Second_Face_Index_In_Cell(axis,it.index)))
                    -face_velocities(FACE_INDEX<TV::m>(axis,mac_grid.First_Face_Index_In_Cell(axis,it.index)));
            div_u(it.index)*=one_over_h;}
        else div_u(it.index)=(T)0;} // dirichlet p cells
}

//#####################################################################
// Function Fix_RHS_Neumann_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Fix_RHS_Neumann_Cells(ARRAY<T,TV_INT>& rhs)
{
    T one_over_h=(T)1/mac_grid.dX.Min();
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+mac_grid.counts));it.Valid();it.Next()){
        if(!cell_dirichlet(it.index) && !cell_neumann(it.index)){ // cell is fluid
            for(int d=0;d<TV::m;d++){
                TV_INT left_index=it.index;left_index(d)--;
                TV_INT right_index=it.index;right_index(d)++;
                if(cell_neumann(left_index))
                    rhs(it.index)-=one_over_h*face_velocities(FACE_INDEX<TV::m>(d,mac_grid.First_Face_Index_In_Cell(d,it.index)));
                if(cell_neumann(right_index))
                    rhs(it.index)+=one_over_h*face_velocities(FACE_INDEX<TV::m>(d,mac_grid.Second_Face_Index_In_Cell(d,it.index)));}}}
}

//#####################################################################
// Function Solve_For_Pressure
// The equation is -div(dt/rho*grad(p))=-div(u)
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Solve_For_Pressure(const T dt,const T rho)
{
    MPM_POISSON_SYSTEM<TV> system(debug_cast<MPM_PROJECTION<TV>&>(*this));
    MPM_POISSON_VECTOR<TV> rhs,x;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    rhs.v.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    x.v.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    KRYLOV_SOLVER<T>::Ensure_Size(vectors,x,3);
    rhs.v.Copy(-1,div_u);
    Fix_RHS_Neumann_Cells(rhs.v);
    x.v=rhs.v; // x=dt/rho*p (for constant rho)
    system.Test_System(*vectors(0),*vectors(1),*vectors(2));
    CONJUGATE_GRADIENT<T> cg;
    CONJUGATE_RESIDUAL<T> cr;
    KRYLOV_SOLVER<T>* solver=&cg;
    solver->print_residuals=false;
    solver->Solve(system,x,rhs,vectors,(T)1e-7,0,1000);
    pressure.Copy(rho/dt,x.v);
    vectors.Delete_Pointers_And_Clean_Memory();
}

//#####################################################################
// Function Do_Projection
// The equation is u=u-dt/rho*grad(p)
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Do_Projection(const T dt,const T rho)
{
    LOG::cout<<"Maximum velocity divergence before projection: "<<div_u.Max_Abs()<<std::endl;        
    T one_over_h=(T)1/mac_grid.dX.Min();
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face_index=iterator.Full_Index();
        int axis=iterator.Axis(); // the axis of first_cell -> second_cell
        TV_INT first_cell=iterator.First_Cell_Index();
        TV_INT second_cell=iterator.Second_Cell_Index();        
        if(first_cell(axis)>=0&&second_cell(axis)<mac_grid.counts(axis)){ // only deal with non-boundary faces
            T grad_p=(pressure(second_cell)-pressure(first_cell))*one_over_h;
            face_velocities(face_index)-=dt/rho*grad_p;}}
    // Enforce face velocities for neumann cells
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),mac_grid.counts));it.Valid();it.Next()){
        if(cell_neumann(it.index)){
            for(int d=0;d<TV::m;d++){
                face_velocities(FACE_INDEX<TV::m>(d,mac_grid.First_Face_Index_In_Cell(d,it.index)))=T(0);
                face_velocities(FACE_INDEX<TV::m>(d,mac_grid.Second_Face_Index_In_Cell(d,it.index)))=T(0);}}}
    // check whether divergence free
    Build_Velocity_Divergence();
    LOG::cout<<"Maximum velocity divergence after projection: "<<div_u.Max_Abs()<<std::endl;
}

//#####################################################################
// Function Velocities_Faces_To_Corners_MPM_Style
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Velocities_Faces_To_Corners_MPM_Style()
{
    for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(nodes_non_dirichlet_cells);it.Valid();it.Next()){
        TV_INT node=it.Key();
        sim.node_V(node)=TV();
        for(int axis=0;axis<TV::m;axis++){
            for(int face=0;face<TV::m*2-2;face++){
                FACE_INDEX<TV::m> face_index(axis,mac_grid.Node_Face_Index(axis,node,face));
                sim.node_V(node)(axis)+=0.5*face_velocities(face_index);}}}
}

//#####################################################################
template class MPM_PROJECTION<VECTOR<float,2> >;
template class MPM_PROJECTION<VECTOR<float,3> >;
template class MPM_PROJECTION<VECTOR<double,2> >;
template class MPM_PROJECTION<VECTOR<double,3> >;
}
