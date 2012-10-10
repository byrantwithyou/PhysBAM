#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include "ACCURACY_INFO.h"
#include "HEADER.h"
#include "POISSON_PROJECTION_SYSTEM.h"
using namespace PhysBAM;

// TODO: It is possible for both cell centers to be inside but that I still need to do something for the Neumann face.

template<class T,class TV> T Face_Fraction(const GRID<TV>& grid,const VECTOR<T,3>& X,const TV& Y,const BOUNDARY_CONDITIONS<TV>& callback){return 0;}

template<class T,class TV>
T Face_Fraction(const GRID<TV>& grid,const VECTOR<T,2>& X,const TV& Y,const BOUNDARY_CONDITIONS<TV>& callback)
{
    TV M=(X+Y)/2,D=(Y-X)/2;
    TV A=M+D.Orthogonal_Vector(),B=M-D.Orthogonal_Vector();
    T pa=callback.Theta(A),pb=callback.Theta(B);
    if(pa>0){
        if(pb>0) return 0;
        return pb/(pb-pa);}
    if(pb<=0) return 1;
    return pa/(pa-pb);
}

template<class T,class TV>
T Face_Fraction(const GRID<TV>& grid,const VECTOR<T,1>& X,const TV& Y,const BOUNDARY_CONDITIONS<TV>& callback)
{
    TV M=(X+Y)/2;
    return callback.Inside(M)?1:0;
}

template<class T>
T L_From_Phi(T a,T b)
{
    if(a<0 && b>0) return a/(a-b);
    if(a>0 && b<0) return b/(b-a);
    if(a<0 && b<0) return 1;
    return 0;
}

template<class T,class TV,class TV_INT,int d>
void Test_Gibou(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<d> >& u,const ARRAY<T,TV_INT>& p,ARRAY<T,FACE_INDEX<d> >& du,const BOUNDARY_CONDITIONS<TV>& callback,T density)
{
    PHYSBAM_FATAL_ERROR();
}

template<class T,class TV,class TV_INT>
void Test_Gibou(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<2> >& u,const ARRAY<T,TV_INT>& p,ARRAY<T,FACE_INDEX<2> >& du,const BOUNDARY_CONDITIONS<TV>& callback,T density)
{
    ERROR_COLOR_MAP<T> color(1e-12,1,true,true,true);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
    {
        TV X=grid.X(it.index);
        T dx=grid.dX.x;
        T dy=grid.dX.y;

        T phi_pp=callback.Theta(X+TV(dx/2,dy/2));
        T phi_pn=callback.Theta(X+TV(dx/2,-dy/2));
        T phi_np=callback.Theta(X+TV(-dx/2,dy/2));
        T phi_nn=callback.Theta(X+TV(-dx/2,-dy/2));

        T L_ne=L_From_Phi(phi_nn,phi_np)*dy;
        T L_pe=L_From_Phi(phi_pn,phi_pp)*dy;
        T L_en=L_From_Phi(phi_nn,phi_pn)*dx;
        T L_ep=L_From_Phi(phi_np,phi_pp)*dx;

        T rho_ne=density;
        T rho_pe=density;
        T rho_en=density;
        T rho_ep=density;

        T p_ee=p(it.index);
        T p_pe=p(it.index+TV_INT(1,0));
        T p_ne=p(it.index+TV_INT(-1,0));
        T p_ep=p(it.index+TV_INT(0,1));
        T p_en=p(it.index+TV_INT(0,-1));

        T u_ne=u(FACE_INDEX<TV::m>(1,it.index));
        T u_pe=u(FACE_INDEX<TV::m>(1,it.index+TV_INT(1,0)));
        T u_en=u(FACE_INDEX<TV::m>(2,it.index));
        T u_ep=u(FACE_INDEX<TV::m>(2,it.index+TV_INT(0,1)));

        T du_ne=1/rho_ne*(p_ee-p_ne)/dx;
        T du_pe=1/rho_pe*(p_pe-p_ee)/dx;
        T du_en=1/rho_en*(p_ee-p_en)/dy;
        T du_ep=1/rho_ep*(p_ep-p_ee)/dy;

        T lhs_ne=L_ne*du_ne;
        T lhs_pe=-L_pe*du_pe;
        T lhs_en=L_en*du_en;
        T lhs_ep=-L_ep*du_ep;

        T rhs_ne=L_ne*u_ne;
        T rhs_pe=-L_pe*u_pe;
        T rhs_en=L_en*u_en;
        T rhs_ep=-L_ep*u_ep;

        T lhs=lhs_ne+lhs_pe+lhs_en+lhs_ep;
        T rhs=rhs_ne+rhs_pe+rhs_en+rhs_ep;

        LOG::cout<<"Cell "<<it.index<<"   error "<<(rhs-lhs)<<"   relative "<<(rhs-lhs)/std::max((T)1e-30,std::max(std::abs(rhs),std::abs(lhs)))<<std::endl;

        if(!du(FACE_INDEX<TV::m>(1,it.index))) du(FACE_INDEX<TV::m>(1,it.index))=du_ne;
        else LOG::cout<<"diff-du A "<<(du(FACE_INDEX<TV::m>(1,it.index))-du_ne)<<std::endl;
        if(!du(FACE_INDEX<TV::m>(1,it.index+TV_INT(1,0)))) du(FACE_INDEX<TV::m>(1,it.index+TV_INT(1,0)))=du_pe;
        else LOG::cout<<"diff-du B "<<(du(FACE_INDEX<TV::m>(1,it.index+TV_INT(1,0)))-du_pe)<<std::endl;
        if(!du(FACE_INDEX<TV::m>(2,it.index))) du(FACE_INDEX<TV::m>(2,it.index))=du_en;
        else LOG::cout<<"diff-du C "<<(du(FACE_INDEX<TV::m>(2,it.index))-du_en)<<std::endl;
        if(!du(FACE_INDEX<TV::m>(2,it.index+TV_INT(0,1)))) du(FACE_INDEX<TV::m>(2,it.index+TV_INT(0,1)))=du_ep;
        else LOG::cout<<"diff-du D "<<(du(FACE_INDEX<TV::m>(2,it.index+TV_INT(0,1)))-du_ep)<<std::endl;

        T u2_ne=u_ne-du_ne;
        T u2_pe=u_pe-du_pe;
        T u2_en=u_en-du_en;
        T u2_ep=u_ep-du_ep;

        T rhs2_ne=L_ne*u2_ne;
        T rhs2_pe=-L_pe*u2_pe;
        T rhs2_en=L_en*u2_en;
        T rhs2_ep=-L_ep*u2_ep;

        T rhs2=rhs2_ne+rhs2_pe+rhs2_en+rhs2_ep;

        Add_Debug_Particle(X,color(abs(rhs2)));

        LOG::cout<<"resid "<<rhs2<<std::endl;
    }
    Dump_Frame<RW>(u,"test");
}

template<class T,class TV,int d>
void Project_Incompressibility_Gibou(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& callback,const ACCURACY_INFO<d>& ai,T time,T density,
    T theta_threshold,T cg_tolerance,bool verbose)
{
    typedef VECTOR<int,d> TV_INT;
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > T_VECTOR;
    typedef KRYLOV::MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<T>,T,T_VECTOR > T_SYSTEM;

    ARRAY<int,TV_INT> cell_to_index(grid.Domain_Indices(1));
    ARRAY<TV_INT> index_to_cell;
    ARRAY<FACE_INDEX<d> > index_to_face;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        if(callback.Inside(grid.X(it.index)))
            cell_to_index(it.index)=index_to_cell.Append(it.index);

    POISSON_PROJECTION_SYSTEM<TV> system;
    ARRAY<T> beta_inverse;
    VECTOR<T,2> sign(-1,1);
    SPARSE_MATRIX_FLAT_MXN<T> neg_div_t;
    neg_div_t.Reset(0);
    system.gradient.Reset(0);

    bool neumann_pocket=true;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<d> face=it.Full_Index();
        VECTOR<TV,2> X;
        T dxi=grid.one_over_dX(it.Axis());
        for(int i=0;i<2;i++) X(i)=grid.X(face.Cell_Index(i));
        T fraction=Face_Fraction(grid,X(1),X(2),callback);
        if(fraction<1e-8) continue;
        index_to_face.Append(face);
        for(int i=0;i<2;i++){
            int& index=cell_to_index(face.Cell_Index(i));
            if(index<0) index=index_to_cell.Append(face.Cell_Index(i));
            system.gradient.Append_Entry_To_Current_Row(index,sign(i)*dxi);
            neg_div_t.Append_Entry_To_Current_Row(index,sign(i)*fraction*dxi);}
        system.gradient.Finish_Row();
        neg_div_t.Finish_Row();}
    system.gradient.n=neg_div_t.n=index_to_cell.m;
    system.gradient.Sort_Entries();
    neg_div_t.Sort_Entries();
    neg_div_t.Transpose(system.neg_divergence);
    system.beta_inverse.Resize(system.gradient.m);
    system.beta_inverse.Fill(1/density);

    system.Initialize();
    system.projections.Append(ARRAY<T>());
    system.projections.Last().Resize(system.poisson.n);
    system.projections.Last().Fill(1/sqrt((T)system.poisson.n));

    T_VECTOR x,b,q,s,r,k,z;
    x.v.Resize(index_to_cell.m);
    b.v.Resize(index_to_cell.m);
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    ARRAY<T> temp(index_to_face.m);
    for(int i=0;i<index_to_face.m;i++) temp(i)=u(index_to_face(i));
    ARRAY<T> old_u=temp;
    system.neg_divergence.Times(temp,b.v);

    ARRAY<T,TV_INT> p(grid.Domain_Indices(1));
    for(int i=0;i<index_to_cell.m;i++) p(index_to_cell(i))=b.v(i);
    ERROR_COLOR_MAP<T> color(1e-12,1,true,true,true);
    for(int i=0;i<index_to_cell.m;i++) Add_Debug_Particle(grid.X(index_to_cell(i)),color(p(index_to_cell(i))));
    Dump_Frame<RW>(u,"divergence");
    ai.Print("DIVERGENCE",p);
    s=b;
    system.Project(s);
    s-=b;
    LOG::cout<<"Nullspace: "<<sqrt(system.Inner_Product(s,s)/system.Inner_Product(b,b))<<std::endl;

    static int solve_id=0;solve_id++;
    if(verbose){
        KRYLOV_SOLVER<T>::Ensure_Size(vectors,x,2);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-temp-%i.txt",solve_id).c_str()).Write("temp",temp);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-M-%i.txt",solve_id).c_str()).Write("M",system,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-G-%i.txt",solve_id).c_str()).Write("G",system.gradient);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-ND-%i.txt",solve_id).c_str()).Write("ND",system.neg_divergence);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-bi-%i.txt",solve_id).c_str()).Write("bi",system.beta_inverse);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-b-%i.txt",solve_id).c_str()).Write("b",b);}

    CONJUGATE_RESIDUAL<T> solver;
    solver.Solve(system,x,b,vectors,cg_tolerance,1,1000000);

    if(verbose){OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("proj-x-%i.txt",solve_id).c_str()).Write("x",x);}

    LOG::cout<<"P-inf "<<INTERVAL<T>::Bounding_Box(x.v).Size()<<std::endl;
    ARRAY<T> tmp(x.v);
    LOG::cout<<"P-1 "<<(tmp-=x.v.Average()).Sum_Abs()/x.v.m<<std::endl;
    for(int i=0;i<index_to_cell.m;i++) p(index_to_cell(i))=x.v(i);
    if(neumann_pocket) p.Subset(ai.cell_samples)-=p.Subset(ai.cell_samples).Average();
//    ai.Print("PRESSURE",p);

    for(int i=0;i<index_to_cell.m;i++) Add_Debug_Particle(grid.X(index_to_cell(i)),color(p(index_to_cell(i))));
    Dump_Frame<RW>(u,"pressures");

//    ARRAY<T,FACE_INDEX<d> > new_u(u.Domain_Indices(),true);
//    Test_Gibou(grid,u,p,new_u,callback,density);
//    new_u-=u;

    system.gradient.Times(x.v,temp);
    temp*=system.beta_inverse;

    for(int i=0;i<index_to_face.m;i++) Add_Debug_Particle(grid.Axis_X_Face(index_to_face(i)),color(temp(i)));
    Dump_Frame<RW>(u,"du");

    for(int i=0;i<index_to_face.m;i++) u(index_to_face(i))-=temp(i);
    old_u-=temp;

/*
    for(int i=0;i<index_to_face.m;i++){
        T t=temp(i);
        T nu=new_u(index_to_face(i));
        T r=(t-nu)/std::max((T)1e-30,std::max(std::abs(t),std::abs(nu)));
        LOG::cout<<"Face "<<index_to_face(i)<<"  t "<<t<<"   nu "<<nu<<"   r "<<r<<std::endl;
        new_u(index_to_face(i))=0;
    }

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        if(!new_u(it.Full_Index())) continue;
        LOG::cout<<"Extra "<<it.Full_Index()<<"  du "<<new_u(it.Full_Index())<<"   phi "<<callback.Theta(it.Location())/grid.dX.Max()<<std::endl;
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
     }
*/

    for(int i=0;i<index_to_face.m;i++){
        if(callback.Inside(grid.Axis_X_Face(index_to_face(i))))
            Add_Debug_Particle(grid.Axis_X_Face(index_to_face(i)),VECTOR<T,3>(0,0,1));
        else Add_Debug_Particle(grid.Axis_X_Face(index_to_face(i)),VECTOR<T,3>(1,1,0));}

    for(int i=0;i<index_to_cell.m;i++) Add_Debug_Particle(grid.X(index_to_cell(i)),VECTOR<T,3>(0,1,0));
    Dump_Frame<RW>(u,"in-out tests");

    for(int i=0;i<index_to_face.m;i++){
        TV X=grid.Axis_X_Face(index_to_face(i));
        T uu=u(index_to_face(i));
        T au=callback.Analytic_Velocity(X,time)(index_to_face(i).axis);
        Add_Debug_Particle(X,color(abs(uu-au)));}
    Dump_Frame<RW>(u,"direct error");

    for(int i=0;i<index_to_face.m;i++) old_u(i)=u(index_to_face(i));
    system.neg_divergence.Times(old_u,s.v);
    for(int i=0;i<index_to_cell.m;i++) Add_Debug_Particle(grid.X(index_to_cell(i)),color(s.v(i)));
    Dump_Frame<RW>(u,"residual divergence");
    LOG::cout<<"MAX RESIDUAL BEFORE EXTRAP "<<s.v.Maximum_Magnitude()<<std::endl;

//     Fill_Ghost_Cells(grid,3,3,u,callback);
//     for(int i=0;i<index_to_face.m;i++) old_u(i)=u(index_to_face(i));
//     system.neg_divergence.Times(old_u,s.v);
//     for(int i=0;i<index_to_cell.m;i++) Add_Debug_Particle(grid.X(index_to_cell(i)),color(s.v(i)));
//     Dump_Frame<RW>(u,"residual divergence (after extrap)");
//     LOG::cout<<"MAX RESIDUAL AFTER EXTRAP "<<s.v.Maximum_Magnitude()<<std::endl;

    LOG::cout<<"G-inf "<<INTERVAL<T>::Bounding_Box(temp).Size()<<std::endl;
    LOG::cout<<"G-1 "<<(temp-=temp.Average()).Sum_Abs()/temp.m<<std::endl;
}

template void Project_Incompressibility_Gibou<double,VECTOR<double,1>,1>(GRID<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> >&,const BOUNDARY_CONDITIONS<VECTOR<double,1> >&,
    const ACCURACY_INFO<1>&,double,double,double,double,bool);
template void Project_Incompressibility_Gibou<double,VECTOR<double,2>,2>(GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> >&,const BOUNDARY_CONDITIONS<VECTOR<double,2> >&,
    const ACCURACY_INFO<2>&,double,double,double,double,bool);
