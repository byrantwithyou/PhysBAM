//#####################################################################
// Copyright 2013
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
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
    face_velocities.Resize(mac_grid);
    LOG::cout<<"mac_grid.counts = "<<mac_grid.counts<<std::endl; // count of cell centers
    cell_dirichlet.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    LOG::cout<<"number of cells via mac_grid = "<<mac_grid.numbers_of_cells<<std::endl;
    LOG::cout<<"number of cells via cell_dirichlet = "<<cell_dirichlet.Size()<<std::endl;
    cell_neumann.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
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
// Function Extrapolate_Velocities_For_Dirichlet_Cells
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Extrapolate_Velocities_For_Dirichlet_Cells()
{
    // Dirichlet p=0, div(u)=0.
    // This function makes oposite face center velocities equal to the
    // larger one.
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),mac_grid.counts));it.Valid();it.Next()){
        if(cell_dirichlet(it.index)){
            for(int axis=0;axis<TV::m;axis++){
                FACE_INDEX<TV::m> first_face_index(axis,mac_grid.First_Face_Index_In_Cell(axis,it.index));
                FACE_INDEX<TV::m> second_face_index(axis,mac_grid.Second_Face_Index_In_Cell(axis,it.index));
                T larger=max(face_velocities(first_face_index),face_velocities(second_face_index));
                face_velocities(first_face_index)=larger;
                face_velocities(second_face_index)=larger;}}}
}

//#####################################################################
// Function Solve_For_Pressure
//#####################################################################
template<class TV> void MPM_PROJECTION<TV>::
Solve_For_Pressure(const T dt,const T rho)
{
    // static int solve_id=-1;solve_id++;
    // MPM_SYSTEM<TV> system(debug_cast<MPM_SIMULATION<TV>&>(*this));
    // MPM_VECTOR<TV> rhs,x;
    // ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    // rhs.v.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    // x.v.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    // KRYLOV_SOLVER<T>::Ensure_Size(vectors,x,3);
    // rhs.v=node_V_star;
    // x.v=rhs.v;
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
 //    node_V_old=node_V;
//     node_V=x.v;
//     vectors.Delete_Pointers_And_Clean_Memory();
}

//#####################################################################
template class MPM_PROJECTION<VECTOR<float,2> >;
template class MPM_PROJECTION<VECTOR<float,3> >;
template class MPM_PROJECTION<VECTOR<double,2> >;
template class MPM_PROJECTION<VECTOR<double,3> >;
}
