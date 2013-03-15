#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
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

template<class T,class TV,int d>
void Project_Incompressibility_Slip(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& callback,const ACCURACY_INFO<d>& ai,T time,T density,
    T theta_threshold,T cg_tolerance,bool verbose)
{
    typedef VECTOR<int,d> TV_INT;
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > T_VECTOR;

    ARRAY<int,TV_INT> cell_to_index(grid.Domain_Indices(1));
    ARRAY<TV_INT> index_to_cell;
    ARRAY<FACE_INDEX<d> > index_to_face;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next())
        if(callback.Inside(grid.X(it.index)))
            cell_to_index(it.index)=index_to_cell.Append(it.index);

    POISSON_PROJECTION_SYSTEM<TV> system;
    VECTOR<T,2> sign(-1,1);
    system.gradient.Reset(0);
    T cell_vol=grid.dX.Product();
    ARRAY<T> fractions;

    bool neumann_pocket=true;
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<d> face=it.Full_Index();
        VECTOR<TV,2> X;
        T dxi=grid.one_over_dX(it.Axis());
        for(int i=0;i<2;i++) X(i)=grid.X(face.Cell_Index(i));
        T fraction=Face_Fraction(grid,X(0),X(1),callback);
        if(fraction<1e-8) continue;
        index_to_face.Append(face);
        fractions.Append(fraction);
        for(int i=0;i<2;i++){
            int& index=cell_to_index(face.Cell_Index(i));
            if(index<0) index=index_to_cell.Append(face.Cell_Index(i));
            system.gradient.Append_Entry_To_Current_Row(index,sign(i)*fraction*dxi*cell_vol);}
        system.gradient.Finish_Row();}
    system.gradient.n=index_to_cell.m;
    system.gradient.Sort_Entries();
    system.gradient.Transpose(system.neg_divergence);
    system.beta_inverse.Resize(fractions.m);
    for(int i=0;i<fractions.m;i++) system.beta_inverse(i)=1/(density*fractions(i)*cell_vol);

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
    system.neg_divergence.Times(temp,b.v);

    ARRAY<T,TV_INT> p(grid.Domain_Indices(1));
    for(int i=0;i<index_to_cell.m;i++) p(index_to_cell(i))=b.v(i);
    ERROR_COLOR_MAP<T> color(1e-12,1,true,true,true);
    for(int i=0;i<index_to_cell.m;i++) Add_Debug_Particle(grid.X(index_to_cell(i)),color(p(index_to_cell(i))));
    Flush_Frame(u,"divergence");
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
    Flush_Frame(u,"pressures");

    system.gradient.Times(x.v,temp);
    temp*=system.beta_inverse;

    for(int i=0;i<index_to_face.m;i++) Add_Debug_Particle(grid.Face(index_to_face(i)),color(temp(i)));
    Flush_Frame(u,"du");

    for(int i=0;i<index_to_face.m;i++) u(index_to_face(i))-=temp(i);

    LOG::cout<<"G-inf "<<INTERVAL<T>::Bounding_Box(temp).Size()<<std::endl;
    LOG::cout<<"G-1 "<<(temp-=temp.Average()).Sum_Abs()/temp.m<<std::endl;
}

template void Project_Incompressibility_Slip<double,VECTOR<double,1>,1>(GRID<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> >&,const BOUNDARY_CONDITIONS<VECTOR<double,1> >&,
    const ACCURACY_INFO<1>&,double,double,double,double,bool);
template void Project_Incompressibility_Slip<double,VECTOR<double,2>,2>(GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> >&,const BOUNDARY_CONDITIONS<VECTOR<double,2> >&,
    const ACCURACY_INFO<2>&,double,double,double,double,bool);
