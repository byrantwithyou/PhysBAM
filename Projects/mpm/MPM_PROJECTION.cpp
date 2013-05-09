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
    mac_grid.Initialize(sim.grid.numbers_of_cells+1,RANGE<TV>(sim.grid.domain.min_corner-sim.grid.dX*0.5,sim.grid.domain.max_corner+sim.grid.dX*0.5),true);
    LOG::cout<<"grid.h = "<<sim.grid.dX.Min()<<std::endl;
    LOG::cout<<"mac_grid.h = "<<mac_grid.dX.Min()<<std::endl;
    face_velocities.Resize(mac_grid);
    LOG::cout<<"mac_grid.counts = "<<mac_grid.counts<<std::endl; // count of cell centers
    cell_dirichlet.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    LOG::cout<<"number of cells via mac_grid = "<<mac_grid.numbers_of_cells<<std::endl;
    LOG::cout<<"number of cells via cell_dirichlet = "<<cell_dirichlet.Size()<<std::endl;
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
    cell_dirichlet.Fill(false);
    cell_neumann.Fill(false);
    div_u.Fill((T)0);
    pressure.Fill((T)0);
}

//#####################################################################
// Function Identify_Dirichlet_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Identify_Dirichlet_Cells()
{
    // Cells without any particles are dirichlet cells
    // (zero pressure, extrapolated face velocity).
    cell_dirichlet.Fill(true);
#pragma omp parallel for
    for(int p=0;p<sim.particles.number;p++)
        cell_dirichlet(mac_grid.Cell(sim.particles.X(p),0))=false;
}

//#####################################################################
// Function Identify_Neumann_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Identify_Neumann_Cells()
{
    // TODO
}

//#####################################################################
// Function Interpolate_Velocities_To_Faces
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Interpolate_Velocities_To_Faces()
{
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face_index=iterator.Full_Index();
        int axis=iterator.Axis(); // the axis of first_cell -> second_cell
        TV_INT first_cell=iterator.First_Cell_Index();
        TV_INT second_cell=iterator.Second_Cell_Index();        
        if(first_cell(axis)>=0&&second_cell(axis)<mac_grid.counts(axis)) // only deal with non-boundary faces
            face_velocities(face_index)=0.5*(sim.node_V(first_cell)(axis)+sim.node_V(second_cell)(axis));
        LOG::cout<<"axis: "<<axis<<"first_cell: "<<first_cell<<"second_cell: "<<second_cell<<"v: "<<face_velocities(face_index)<<std::endl;}
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
// Function Solve_For_Pressure
// The equation is -div(dt/rho*grad(p))=-div(u)
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Solve_For_Pressure(const T dt,const T rho)
{
    // static int solve_id=-1;solve_id++;
    // MPM_POISSON_SYSTEM<TV> system(debug_cast<MPM_PROJECTION<TV>&>(*this));
    // MPM_POISSON_VECTOR<TV> rhs,x;
    // ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    // rhs.v.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    // x.v.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    // KRYLOV_SOLVER<T>::Ensure_Size(vectors,x,3);
    
    // rhs.v.
    // x.v=rhs.v; // x=dt/rho*p (for constant rho)
    // if(test_system) system.Test_System(*vectors(0),*vectors(1),*vectors(2));
    // CONJUGATE_GRADIENT<T> cg;
    // CONJUGATE_RESIDUAL<T> cr;
    // KRYLOV_SOLVER<T>* solver=&cg;
    // solver->print_residuals=true;
    // if(dump_matrix){ // load M-1.txt;load m-1.txt;mm=reshape([m m]',rows(M),1);R=diag(mm)*M;max(max(abs(R-R')))
    //     LOG::cout<<"solve id "<<solve_id<<std::endl;
    //     OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",system,*vectors(0),*vectors(1));
    //     OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("m-%i.txt",solve_id).c_str()).Write("m",node_mass.array);
    //     OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",rhs);} 
    // solver->Solve(system,x,rhs,vectors,(T)1e-7,0,1000);
    // if(dump_matrix){
    //     LOG::cout<<"solve id "<<solve_id<<std::endl;
    //     OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",x);}
    // pressure=x.v*rho/dt;
    // vectors.Delete_Pointers_And_Clean_Memory();
}

//#####################################################################
// Function Do_Projection
// The equation is u=u-dt/rho*grad(p)
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Do_Projection(const T dt,const T rho)
{
    T one_over_h=(T)1/mac_grid.dX.Min();
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face_index=iterator.Full_Index();
        int axis=iterator.Axis(); // the axis of first_cell -> second_cell
        TV_INT first_cell=iterator.First_Cell_Index();
        TV_INT second_cell=iterator.Second_Cell_Index();        
        if(first_cell(axis)>=0&&second_cell(axis)<mac_grid.counts(axis)){ // only deal with non-boundary faces
            T grad_p=(pressure(second_cell)-pressure(first_cell))*one_over_h;
            face_velocities(face_index)-=dt/rho*grad_p;}}
}

//#####################################################################
template class MPM_PROJECTION<VECTOR<float,2> >;
template class MPM_PROJECTION<VECTOR<float,3> >;
template class MPM_PROJECTION<VECTOR<double,2> >;
template class MPM_PROJECTION<VECTOR<double,3> >;
}
