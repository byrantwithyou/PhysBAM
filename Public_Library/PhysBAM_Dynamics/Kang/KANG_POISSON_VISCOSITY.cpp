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
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
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
//#####################################################################
// Function Viscosity_Jump
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m> KANG_POISSON_VISCOSITY<TV>::
Viscosity_Jump(const FACE_INDEX<TV::m>& face) const
{
    return (T).5*(Viscosity_Jump(face.First_Cell_Index())+Viscosity_Jump(face.Second_Cell_Index()));
}
//#####################################################################
// Function Viscosity_Jump
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m> KANG_POISSON_VISCOSITY<TV>::
Viscosity_Jump(const TV_INT& cell) const
{
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET LEVELSET;
    MATRIX<T,TV::m> du;
    TV N=LEVELSET::Normal_At_Node(*fluids_parameters.grid,old_phi,cell);

    for(int d=1;d<=TV::m;d++){
        TV_INT cellp=cell,celln=cell;
        cellp(d)++;
        celln(d)--;
        TV v1=face_velocities_ghost.Cell_Centered_Average(celln);
        TV v2=face_velocities_ghost.Cell_Centered_Average(cellp);
        du.Column(d)=(T).5*fluids_parameters.grid->one_over_dX(d)*(v2-v1);}

    T dmu=fluids_parameters.outside_viscosity-fluids_parameters.viscosity;
    MATRIX<T,TV::m> NN=MATRIX<T,TV::m>::Outer_Product(N,N),TT=(T)1-NN;
    MATRIX<T,TV::m> J=dmu*(du*TT+NN*du*NN-NN*du*TT);
    LOG::cout<<cell<<"     "<<J<<std::endl;
    return J;
}
//#####################################################################
// Function Viscosity_Jump
//#####################################################################
template<class TV> typename TV::SCALAR KANG_POISSON_VISCOSITY<TV>::
Face_Phi(const FACE_INDEX<TV::m>& face) const
{
    return (T).5*(old_phi(face.First_Cell_Index())+old_phi(face.Second_Cell_Index()));
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
Project_Fluid(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt) const
{
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
    fluids_parameters.incompressible->boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,0,3);

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
            face_velocities(it.Full_Index())=psi_N_value(it.Full_Index());
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

    if(test_system) system.Test_System(r,k,z);

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
    color_map.Initialize_Colors(-.3,.3,false,true,false);

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
Apply_Viscosity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt,bool implicit) const
{
    face_velocities_ghost.Resize(*fluids_parameters.grid,3,false);
    fluids_parameters.incompressible->boundary->Fill_Ghost_Cells_Face(*fluids_parameters.grid,face_velocities,face_velocities_ghost,0,3);
    for(int d=1;d<=TV::m;d++){
        Apply_Viscosity(face_velocities,d,dt,implicit);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after viscosity direction",0,1);}
}
//#####################################################################
// Function Apply_Viscosity
//#####################################################################
template<class TV> void KANG_POISSON_VISCOSITY<TV>::
Apply_Viscosity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,int axis,T dt,bool implicit) const
{
    const GRID<TV>& grid=*fluids_parameters.grid;
    const GRID<TV> dual_grid=fluids_parameters.grid->Get_Face_Grid(axis).Get_MAC_Grid_At_Regular_Positions();

    int num_dual_cells=0;
    ARRAY<int,TV_INT> dual_cell_index(dual_grid.Domain_Indices(3));
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,0,GRID<TV>::WHOLE_REGION,0,axis);it.Valid();it.Next()){
        if(psi_N(it.Full_Index())) face_velocities(it.Full_Index())=psi_N_value(it.Full_Index());
        else if(!(psi_D(it.First_Cell_Index()) && psi_D(it.Second_Cell_Index())))
            dual_cell_index(it.index)=++num_dual_cells;}

    SYSTEM_MATRIX_HELPER<T> helper;
    KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > u,b,q,s,r,k,z;
    u.v.Resize(num_dual_cells);
    b.v.Resize(num_dual_cells);
    r.v.Resize(num_dual_cells); // used for density before solve
    q.v.Resize(num_dual_cells);
    if(implicit){
        s.v.Resize(num_dual_cells);
        k.v.Resize(num_dual_cells);
        z.v.Resize(num_dual_cells);}

    T mu_n=fluids_parameters.viscosity,mu_p=fluids_parameters.outside_viscosity;
    T rho_n=fluids_parameters.density,rho_p=fluids_parameters.outside_density;
    TV one_over_dX_2=sqr(dual_grid.one_over_dX);

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,0,GRID<TV>::WHOLE_REGION,0,axis);it.Valid();it.Next()){
        if(int index=dual_cell_index(it.index)){
            b.v(index)=face_velocities(it.Full_Index());
            r.v(index)=Face_Phi(it.Full_Index())<0?rho_n:rho_p;}}
    if(!implicit) u.v=b.v;

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(dual_grid);it.Valid();it.Next()){
        TV_INT dual_cell1=it.First_Cell_Index(),dual_cell2=it.Second_Cell_Index();
        FACE_INDEX<TV::m> face1(axis,dual_cell1),face2(axis,dual_cell2);
        int index1=dual_cell_index(dual_cell1),index2=dual_cell_index(dual_cell2);
        if(!index1 && !index2){Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));continue;}
        if(it.Axis()==axis){if(psi_D(dual_cell1)){Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));continue;}}
        else{
            if(!index1 && !psi_N(face1)){Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,1));Add_Debug_Particle(grid.Axis_X_Face(face1),VECTOR<T,3>(1,0,1));continue;}
            if(!index2 && !psi_N(face2)){Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,1));Add_Debug_Particle(grid.Axis_X_Face(face2),VECTOR<T,3>(1,0,1));continue;}}
        if(!index1 || !index2){Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,0,1));}

        T phi1=Face_Phi(face1),phi2=Face_Phi(face2);
        T mu_hat,mu1=(phi1<0)?mu_n:mu_p,muj=0;

        if((phi1<0)==(phi2<0)){
            mu_hat=mu1;}
        else{
            T mu2=(phi2<0)?mu_n:mu_p;
            T theta=phi1/(phi1-phi2);
            mu_hat=mu1*mu2/(theta*mu1+(1-theta)*mu2);
            muj=Viscosity_Jump(face1)(axis,it.Axis())*(1-theta)+Viscosity_Jump(face2)(axis,it.Axis())*theta;
            if(phi1<0) muj=-muj;
            muj*=mu_hat*dual_grid.one_over_dX(it.Axis());
            if(index1) b.v(index1)+=muj*(1-theta)/(mu2*r.v(index1));
            if(index2) b.v(index2)+=muj*theta/(mu1*r.v(index2));}

        T A=dt*mu_hat*one_over_dX_2(it.Axis()); // positive for diagonal, neg for off

        if(!index1) b.v(index2)+=A*psi_N_value(face1)/r.v(index2);
        if(!index2) b.v(index1)+=A*psi_N_value(face2)/r.v(index1);

        if(index1) helper.data.Append(TRIPLE<int,int,T>(index1,index1,A));
        if(index1 && index2) helper.data.Append(TRIPLE<int,int,T>(index1,index2,-A));
        if(index1 && index2) helper.data.Append(TRIPLE<int,int,T>(index2,index1,-A));
        if(index2) helper.data.Append(TRIPLE<int,int,T>(index2,index2,A));}

    if(implicit) for(int i=1;i<=r.v.n;i++) helper.data.Append(TRIPLE<int,int,T>(i,i,r.v(i)));

    SPARSE_MATRIX_FLAT_NXN<T> matrix;
    helper.Compact();
    helper.Set_Matrix(num_dual_cells,matrix);

    static int solve_id=0;solve_id++;
    if(print_matrix){
        LOG::cout<<"viscosity id "<<solve_id<<std::endl;
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-rho-%i.txt",solve_id).c_str()).Write("rho",r.v);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-M-%i.txt",solve_id).c_str()).Write("M",matrix);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-b-%i.txt",solve_id).c_str()).Write("b",b.v);}

    if(implicit){
        typedef KRYLOV::MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_NXN<T>,T,KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > > SYSTEM;
        matrix.Construct_Incomplete_Cholesky_Factorization();
        SYSTEM system(matrix);
        system.P=matrix.C;
        b.v*=r.v;
        if(print_matrix && system.P) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-Z-%i.txt",solve_id).c_str()).Write("Z",*system.P);
        CONJUGATE_RESIDUAL<T> cr;
        CONJUGATE_GRADIENT<T> cg;
        SYMMQMR<T> qm;
        KRYLOV_SOLVER<T>* solver=&cg;
        solver->restart_iterations=fluids_parameters.cg_restart_iterations;
        if(test_system) system.Test_System(r,k,z);
        solver->Solve(system,u,b,q,s,r,k,z,fluids_parameters.incompressible_tolerance,0,fluids_parameters.incompressible_iterations);
    }
    else{
        matrix.Times(u.v,q.v);
        q.v/=r.v;
        u.v=b.v-q.v;}
    if(print_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-u-%i.txt",solve_id).c_str()).Write("u",u.v);

    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(u.v.Min(),u.v.Max(),false,true,false);

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,0,GRID<TV>::WHOLE_REGION,0,axis);it.Valid();it.Next()){
        if(int index=dual_cell_index(it.index)){
            face_velocities(it.Full_Index())=u.v(index);
            Add_Debug_Particle(it.Location(),color_map(u.v(index)));}}
}
template class KANG_POISSON_VISCOSITY<VECTOR<float,1> >;
template class KANG_POISSON_VISCOSITY<VECTOR<float,2> >;
template class KANG_POISSON_VISCOSITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KANG_POISSON_VISCOSITY<VECTOR<double,1> >;
template class KANG_POISSON_VISCOSITY<VECTOR<double,2> >;
template class KANG_POISSON_VISCOSITY<VECTOR<double,3> >;
#endif
