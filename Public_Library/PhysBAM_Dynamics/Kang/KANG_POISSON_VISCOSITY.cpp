//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KANG_POISSON_VISCOSITY
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
#include <PhysBAM_Dynamics/Kang/KANG_POISSON_VISCOSITY.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> KANG_POISSON_VISCOSITY<TV>::
KANG_POISSON_VISCOSITY(FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters_input,const ARRAY<T,TV_INT>& old_phi_input)
    :fluids_parameters(fluids_parameters_input),old_phi(old_phi_input),print_matrix(false),test_system(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> KANG_POISSON_VISCOSITY<TV>::
~KANG_POISSON_VISCOSITY()
{
}
//#####################################################################
// Function Cell_Index
//#####################################################################
template<class TV> int KANG_POISSON_VISCOSITY<TV>::
Cell_Index(const TV_INT& cell) const
{
    const TV_INT& counts=fluids_parameters.grid->counts;
    int r=cell(1)-1;
    for(int i=2; i<=TV::m; i++)
        r=r*counts(i)+(cell(i)-1);
    return r+1;
}
//#####################################################################
// Function Pressure_Jump
//#####################################################################
template<class TV> typename TV::SCALAR KANG_POISSON_VISCOSITY<TV>::
Pressure_Jump(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const TV_INT& cell,T dt) const
{
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET LEVELSET;
    const LEVELSET& phi=fluids_parameters.particle_levelset_evolution->Levelset(1);
    T kappa=phi.Compute_Curvature(cell);
    T pj_st=dt*fluids_parameters.surface_tension*kappa;
    TV du_n;
    TV N=LEVELSET::Normal_At_Node(*fluids_parameters.grid,old_phi,cell);

    for(int d=1;d<=TV::m;d++){
        TV_INT cellp=cell,celln=cell;
        cellp(d)++;
        celln(d)--;
        TV v1=face_velocities.Cell_Centered_Average(celln);
        TV v2=face_velocities.Cell_Centered_Average(cellp);
        du_n(d)=(T).5*fluids_parameters.grid->one_over_dX(d)*TV::Dot_Product(v2-v1,N);}

    T n_du_n=TV::Dot_Product(du_n,N);
    T pj_mu=2*dt*(fluids_parameters.outside_viscosity-fluids_parameters.viscosity)*n_du_n;
    return pj_st+pj_mu;
}
//#####################################################################
// Function Project_Fluid
//#####################################################################
template<class TV> void KANG_POISSON_VISCOSITY<TV>::
Project_Fluid(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt,T time) const
{
    // TODO: phi for normals must be phi^n
    const ARRAY<T,TV_INT>& phi=fluids_parameters.particle_levelset_evolution->Levelset(1).phi;
    const GRID<TV>& grid=*fluids_parameters.grid;

    SYSTEM_MATRIX_HELPER<T> helper;
    ARRAY<TRIPLE<int,int,T> > grad_helper;
    int num_cells=grid.counts.Product();
    KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > p,rhs,q,s,r,k,z;
    p.v.Resize(num_cells);
    rhs.v.Resize(num_cells);
    q.v.Resize(num_cells);
    s.v.Resize(num_cells);
    r.v.Resize(num_cells);
    k.v.Resize(num_cells);
    z.v.Resize(num_cells);

    T beta_n=1/fluids_parameters.density,beta_p=1/fluids_parameters.outside_density;
    TV one_over_dx_2=sqr(grid.one_over_dX);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        TV_INT cell1=it.First_Cell_Index(),cell2=it.Second_Cell_Index();
        T phi1=phi(cell1),phi2=phi(cell2);
        T beta_hat,rhs1=face_velocities(it.Full_Index())*grid.one_over_dX(it.Axis());
        T beta1=(phi1<0)?beta_n:beta_p;
        if((phi1<0)==(phi2<0)){
            beta_hat=beta1;}
        else{
            T beta2=(phi2<0)?beta_n:beta_p;
            T theta=phi1/(phi1-phi2);
            beta_hat=beta1*beta2/(theta*beta1+(1-theta)*beta2);
            T pj=Pressure_Jump(face_velocities,cell1,dt)*(1-theta)+Pressure_Jump(face_velocities,cell2,dt)*theta;
            rhs1-=pj*beta_hat*grid.one_over_dX(it.Axis());}

        T A=beta_hat*one_over_dx_2(it.Axis()); // positive for diagonal, neg for off
        int index1=Cell_Index(cell1),index2=Cell_Index(cell2);
        helper.data.Append(TRIPLE<int,int,T>(index1,index1,A));
        helper.data.Append(TRIPLE<int,int,T>(index1,index2,-A));
        helper.data.Append(TRIPLE<int,int,T>(index2,index1,-A));
        helper.data.Append(TRIPLE<int,int,T>(index2,index2,A));

        T beta_G=beta_hat*grid.one_over_dX(it.Axis()); // positive for 2, neg for 1
        grad_helper.Append(TRIPLE<int,int,T>(index1,index2,beta_G));

        rhs.v(index1)-=rhs1; // Setting up negative of standard system
        rhs.v(index2)+=rhs1;}

    SPARSE_MATRIX_FLAT_NXN<T> matrix;
    helper.Compact();
    helper.Set_Matrix(num_cells,matrix);
    typedef KRYLOV::MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_NXN<T>,T,VECTOR_ND<T> > SYSTEM;
    matrix.Construct_Incomplete_Cholesky_Factorization();
    SYSTEM system(matrix);
    system.P=matrix.C;

    if(test_system) system.Test_System(rhs,k,z);

    static int solve_id=1;solve_id++;
    if(print_matrix){
        LOG::cout<<"pressure solve id "<<solve_id<<std::endl;
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",matrix);
        if(system.P) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("Z-%i.txt",solve_id).c_str()).Write("Z",*system.P);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",rhs.v);}

    CONJUGATE_RESIDUAL<T> cr;
    CONJUGATE_GRADIENT<T> cg;
    SYMMQMR<T> qm;
    KRYLOV_SOLVER<T>* solver=&cg;
    solver->restart_iterations=fluids_parameters.cg_restart_iterations;
    solver->Solve(system,p,rhs,q,s,r,k,z,fluids_parameters.incompressible_tolerance,0,fluids_parameters.incompressible_iterations);
    if(print_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",p.v);

    int i=1;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next(),i++){
        T dp=p.v(grad_helper(i).y)-p.v(grad_helper(i).x);
        face_velocities(it.Full_Index())-=grad_helper(i).z*dp;}
}
//#####################################################################
// Function Apply_Viscosity
//#####################################################################
template<class TV> void KANG_POISSON_VISCOSITY<TV>::
Apply_Viscosity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt,T time) const
{
}
template class KANG_POISSON_VISCOSITY<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KANG_POISSON_VISCOSITY<VECTOR<double,2> >;
#endif
