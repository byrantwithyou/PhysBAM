//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KANG_POISSON_VISCOSITY
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
#include <PhysBAM_Dynamics/Kang/KANG_POISSON_VISCOSITY.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}

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
// Function Pressure_Jump
//#####################################################################
template<class TV> typename TV::SCALAR KANG_POISSON_VISCOSITY<TV>::
Pressure_Jump(const TV_INT& cell,T dt) const
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
        TV v1=face_velocities_ghost.Cell_Centered_Average(celln);
        TV v2=face_velocities_ghost.Cell_Centered_Average(cellp);
        du_n(d)=(T).5*fluids_parameters.grid->one_over_dX(d)*TV::Dot_Product(v2-v1,N);}

    T n_du_n=TV::Dot_Product(du_n,N);
    T pj_mu=2*dt*(fluids_parameters.outside_viscosity-fluids_parameters.viscosity)*n_du_n;
    return pj_st+pj_mu;
}
namespace{
template<class T>
struct GRAD_HELPER
{
    int i,j;
    T x;
};
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

    int num_cells=0;
    ARRAY<int,TV_INT> cell_index(grid.Domain_Indices(3));
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        if(!psi_D(it.index))
            cell_index(it.index)=++num_cells;

    SYSTEM_MATRIX_HELPER<T> helper;
    ARRAY<GRAD_HELPER<T> > grad_helper;
    KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > p,rhs,q,s,r,k,z;
    p.v.Resize(num_cells);
    rhs.v.Resize(num_cells);
    q.v.Resize(num_cells);
    s.v.Resize(num_cells);
    r.v.Resize(num_cells);
    k.v.Resize(num_cells);
    z.v.Resize(num_cells);

    face_velocities_ghost.Resize(*fluids_parameters.grid,3,false);
    fluids_parameters.incompressible->boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time+dt,3);

    T beta_n=1/fluids_parameters.density,beta_p=1/fluids_parameters.outside_density;
    TV one_over_dX_2=sqr(grid.one_over_dX);
    T mn=FLT_MAX,mx=-mn;
    GRAD_HELPER<T> null_helper = {0,0,0};

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        TV_INT cell1=it.First_Cell_Index(),cell2=it.Second_Cell_Index();
        bool psi_n=psi_N(it.Full_Index());
        int index1=cell_index(cell1),index2=cell_index(cell2);

        if(!index1 && !index2){
            grad_helper.Append(null_helper);
            continue;}
        if(psi_n){
            T rhs1=psi_N_value(it.Full_Index())*grid.one_over_dX(it.Axis());
            if(index1) rhs.v(index1)+=rhs1; // TODO: Check signs here
            if(index2) rhs.v(index2)-=rhs1;
            grad_helper.Append(null_helper);
            continue;}

        T phi1=phi(cell1),phi2=phi(cell2);
        T beta_hat,rhs1=face_velocities(it.Full_Index())*grid.one_over_dX(it.Axis());
        T beta1=(phi1<0)?beta_n:beta_p;
        T pj=0;

        if((phi1<0)==(phi2<0)){
            beta_hat=beta1;}
        else{
            T beta2=(phi2<0)?beta_n:beta_p;
            T theta=phi1/(phi1-phi2);
            beta_hat=beta1*beta2/(theta*beta1+(1-theta)*beta2);
            pj=Pressure_Jump(cell1,dt)*(1-theta)+Pressure_Jump(cell2,dt)*theta;
            if(pj>mx) mx=pj;
            if(pj<mn) mn=pj;
            if(phi1<0) pj=-pj;
            rhs1-=pj*beta_hat*one_over_dX_2(it.Axis());}

        T A=beta_hat*one_over_dX_2(it.Axis()); // positive for diagonal, neg for off

        if(!index1) rhs1+=A*psi_D_value(cell1);
        if(!index2) rhs1-=A*psi_D_value(cell2);

        if(index1) helper.data.Append(TRIPLE<int,int,T>(index1,index1,A));
        if(index1 && index2) helper.data.Append(TRIPLE<int,int,T>(index1,index2,-A));
        if(index1 && index2) helper.data.Append(TRIPLE<int,int,T>(index2,index1,-A));
        if(index2) helper.data.Append(TRIPLE<int,int,T>(index2,index2,A));

        T beta_G=beta_hat*grid.one_over_dX(it.Axis()); // positive for 2, neg for 1
        GRAD_HELPER<T> gh = {index1,index2,beta_G};

        grad_helper.Append(gh);
        face_velocities(it.Full_Index())-=beta_G*pj;

        if(index1) rhs.v(index1)-=rhs1; // Setting up negative of standard system
        if(index2) rhs.v(index2)+=rhs1;}

    LOG::cout<<"pressure range: "<<mn<<"  "<<mx<<std::endl;
    SPARSE_MATRIX_FLAT_NXN<T> matrix;
    helper.Compact();
    helper.Set_Matrix(num_cells,matrix);
    typedef KRYLOV::MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_NXN<T>,T,KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > > SYSTEM;
    matrix.Construct_Incomplete_Cholesky_Factorization();
    SYSTEM system(matrix);
    system.P=matrix.C;

    if(test_system) system.Test_System(rhs,k,z);

    static int solve_id=0;solve_id++;
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

    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(-.3,.3,false,false,false);

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int ci=cell_index(it.index);
        T val=p.v(ci);
        Add_Debug_Particle(it.Location(), color_map(val));}

    int i=1;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next(),i++){
        T dp=0;
        if(grad_helper(i).i) dp-=p.v(grad_helper(i).i);
        if(grad_helper(i).j) dp+=p.v(grad_helper(i).j);
        face_velocities(it.Full_Index())-=grad_helper(i).x*dp;}
}
//#####################################################################
// Function Apply_Viscosity
//#####################################################################
template<class TV> void KANG_POISSON_VISCOSITY<TV>::
Apply_Viscosity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt,T time) const
{
}
template class KANG_POISSON_VISCOSITY<VECTOR<float,1> >;
template class KANG_POISSON_VISCOSITY<VECTOR<float,2> >;
template class KANG_POISSON_VISCOSITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KANG_POISSON_VISCOSITY<VECTOR<double,1> >;
template class KANG_POISSON_VISCOSITY<VECTOR<double,2> >;
template class KANG_POISSON_VISCOSITY<VECTOR<double,3> >;
#endif
