#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include "ACCURACY_INFO.h"
#include "HEADER.h"
#include "POISSON_PROJECTION_SYSTEM.h"
using namespace PhysBAM;

template<class T,class TV,int d>
void Project_Incompressibility(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& callback,const ACCURACY_INFO<d>& ai,T time,T density,
    T theta_threshold,T cg_tolerance,bool verbose)
{
    typedef VECTOR<int,d> TV_INT;
    typedef KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > T_VECTOR;
    typedef KRYLOV::MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<T>,T,T_VECTOR > T_SYSTEM;

    Fill_Ghost_Cells(grid,3,3,u,callback);

    ARRAY<int,TV_INT> cell_to_index(grid.Domain_Indices());
    ARRAY<TV_INT> index_to_cell;
    ARRAY<FACE_INDEX<d> > index_to_face;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        if(callback.Inside(grid.X(it.index)))
            cell_to_index(it.index)=index_to_cell.Append(it.index);

    POISSON_PROJECTION_SYSTEM<TV> system;
    ARRAY<T> beta_inverse;
    system.gradient.Reset(index_to_cell.m);
    VECTOR<T,2> sign(-1,1);
    VECTOR_ND<T> rhs(index_to_cell.m);

    bool neumann_pocket=true;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<d> face=it.Full_Index();
        VECTOR<TV_INT,2> cell;
        VECTOR<TV,2> X;
        VECTOR<bool,2> inside;
        T dxi=grid.one_over_dX(it.Axis());
        for(int i=0;i<2;i++){
            cell(i)=face.Cell_Index(i);
            X(i)=grid.X(cell(i));
            inside(i)=callback.Inside(X(i));
            if(inside(i)) system.gradient.Append_Entry_To_Current_Row(cell_to_index(cell(i)),sign(i)*dxi);}

        if(!inside.Contains(true)) continue;
        system.gradient.Finish_Row();
        index_to_face.Append(face);

        if(!inside.Contains(false)){beta_inverse.Append(1);continue;}

        int k=inside(1)?1:2;
        T theta=1;
        TV value;
        int type=callback.Boundary_Condition(X(k),X(3-k),theta,value,time);
        if(type==BOUNDARY_CONDITIONS<TV>::noslip){
            beta_inverse.Append(0);
            rhs(cell_to_index(cell(k)))-=sign(k)*dxi*value(it.Axis());
            continue;}
        if(type==BOUNDARY_CONDITIONS<TV>::free){
            neumann_pocket=false;
            beta_inverse.Append(1/max(theta,theta_threshold));
            rhs(cell_to_index(cell(k)))-=dxi*dxi*value(it.Axis());
            continue;}
        PHYSBAM_FATAL_ERROR();}
    system.gradient.Sort_Entries();
    system.beta_inverse.Resize(beta_inverse.m);
    for(int i=0;i<beta_inverse.m;i++) system.beta_inverse(i)=beta_inverse(i)/density;
    system.Initialize();
    if(neumann_pocket){
        system.projections.Append(VECTOR_ND<T>());
        system.projections.Last().Resize(system.poisson.n);
        system.projections.Last().Fill(1);}

    T_VECTOR x,b,q,s,r,k,z;
    x.v.Resize(index_to_cell.m);
    b.v.Resize(index_to_cell.m);
    q.v.Resize(index_to_cell.m);
    s.v.Resize(index_to_cell.m);
    r.v.Resize(index_to_cell.m);
    k.v.Resize(index_to_cell.m);
    z.v.Resize(index_to_cell.m);

    VECTOR_ND<T> temp(index_to_face.m);
    for(int i=0;i<index_to_face.m;i++) temp(i)=u(index_to_face(i));
    system.gradient.Transpose_Times(temp,b.v);
    b.v-=rhs; // rhs set up based on a negative definite Poisson system.

    ARRAY<T,TV_INT> p(grid.Domain_Indices());
    for(int i=0;i<index_to_cell.m;i++) p(index_to_cell(i))=b.v(i);
    ai.Print("DIVERGENCE",p);

    static int solve_id=0;solve_id++;
    if(verbose){
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-temp-%i.txt",solve_id).c_str()).Write("temp",temp);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-M-%i.txt",solve_id).c_str()).Write("M",system,q,s);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-G-%i.txt",solve_id).c_str()).Write("G",system.gradient);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-ND-%i.txt",solve_id).c_str()).Write("ND",system.neg_divergence);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-bi-%i.txt",solve_id).c_str()).Write("bi",system.beta_inverse);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-b-%i.txt",solve_id).c_str()).Write("b",b);}

    CONJUGATE_RESIDUAL<T> solver;
    solver.Solve(system,x,b,q,s,r,k,z,cg_tolerance,1,1000000);

    if(verbose){OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-x-%i.txt",solve_id).c_str()).Write("x",x);}

    for(int i=0;i<index_to_cell.m;i++) p(index_to_cell(i))=x.v(i);
    if(neumann_pocket) p.Subset(ai.cell_samples)-=p.Subset(ai.cell_samples).Average();
    ai.Print("PRESSURE",p);

    system.gradient.Times(x.v,temp);
    temp*=system.beta_inverse;

    for(int i=0;i<index_to_face.m;i++) u(index_to_face(i))-=temp(i);
}

template void Project_Incompressibility<double,VECTOR<double,1>,1>(GRID<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> >&,const BOUNDARY_CONDITIONS<VECTOR<double,1> >&,
    const ACCURACY_INFO<1>&,double,double,double,double,bool);
template void Project_Incompressibility<double,VECTOR<double,2>,2>(GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> >&,const BOUNDARY_CONDITIONS<VECTOR<double,2> >&,
    const ACCURACY_INFO<2>&,double,double,double,double,bool);
