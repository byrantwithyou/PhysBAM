#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include "HEADER.h"
using namespace PhysBAM;

template<class T,class TV> T Face_Fraction(const GRID<TV>& grid,const VECTOR<T,3>& X,const TV& Y,BOUNDARY_CONDITIONS<TV>& callback){return 0;}
template<class T,class TV> T Face_Fraction(const GRID<TV>& grid,const VECTOR<T,1>& X,const TV& Y,BOUNDARY_CONDITIONS<TV>& callback){return 0;}

template<class T,class TV>
T Face_Fraction(const GRID<TV>& grid,const VECTOR<T,2>& X,const TV& Y,BOUNDARY_CONDITIONS<TV>& callback)
{
    TV M=(X+Y)/2,D=(Y-X)/2;
    TV A=M+D.Orthogonal_Vector(),B=M-D.Orthogonal_Vector();
    bool ain=callback.Inside(A),bin=callback.Inside(B);
    if(ain && bin) return 1;
    if(!ain && !bin) return 0;
    T theta=1,value=0;
    PHYSBAM_ASSERT(callback.Boundary_Conditions(A,B,theta,value)==BOUNDARY_CONDITIONS<TV>::free);
    if(ain) return theta;
    return 1-theta;
}

template<class T,class TV,int d>
void Apply_Viscosity(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& callback,T dt,T time,T viscosity,T density,int axis,T theta_threshold,T cg_tolerance,bool verbose)
{
    typedef KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > T_VECTOR;
    typedef KRYLOV::MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_NXN<T>,T,T_VECTOR > T_SYSTEM;

    ARRAY<int,FACE_INDEX<d> > index(u.Domain_Indices());
    ARRAY<FACE_INDEX<d> > faces;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,0,GRID<TV>::WHOLE_REGION,-1,axis);it.Valid();it.Next())
        if(callback.Inside(grid.Axis_X_Face(it.Full_Index())))
            index(it.Full_Index())=faces.Append(it.Full_Index());

    SPARSE_MATRIX_FLAT_NXN<T> P;
    P.Reset();
    T_VECTOR x,b,q,s,r,k,z;
    x.v.Resize(faces.m);
    b.v.Resize(faces.m);
    q.v.Resize(faces.m);
    s.v.Resize(faces.m);
    r.v.Resize(faces.m);
    k.v.Resize(faces.m);
    z.v.Resize(faces.m);

    TV one_over_dX2=grid.one_over_dX*grid.one_over_dX;
    for(int i=0;i<faces.m;i++){
        FACE_INDEX<d> face=faces(i);
        TV X=grid.Axis_X_Face(face);
        if(!callback.Inside(X)) continue;
        T middle_value=-2*one_over_dX2.Sum();
        for(int a=0;a<d;a++){
            for(int s=-1;s<=1;s+=2){
                FACE_INDEX<d> neighbor=face;
                neighbor.index(a)+=s;
                TV Y=grid.Axis_X_Face(neighbor);
                if(callback.Inside(Y)){P.Append_Entry_To_Current_Row(index(neighbor),one_over_dX2(a));continue;}

                T theta=1;
                TV value;
                int type=callback.Boundary_Condition(X,Y,theta,value,time+dt);
                if(type==BOUNDARY_CONDITIONS<TV>::free){middle_value+=one_over_dX2(a)*(1+value(axis)*grid.dX(a));}
                else if(type==BOUNDARY_CONDITIONS<TV>::noslip){
                    if(theta<theta_threshold) theta=theta_threshold;
                    middle_value+=one_over_dX2(a)*(1-1/theta);
                    b.v(i)-=one_over_dX2(a)*value(axis)/theta;}
                else PHYSBAM_FATAL_ERROR();}}
        P.Append_Entry_To_Current_Row(i,middle_value);
        P.Finish_Row();}
    P.Sort_Entries();

    b*=-dt*viscosity/density;
    P*=-dt*viscosity/density;
    P+=(T)1;

    for(int i=0;i<faces.m;i++) b.v(i)+=u(faces(i));

    T_SYSTEM system(P);
    T_VECTOR system_vector;
    system_vector.v.Resize(P.n);
    P.Construct_Incomplete_Cholesky_Factorization();
    system.Set_Preconditioner(*P.C,system_vector);

    static int solve_id=0;solve_id++;
    if(verbose){
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-M-%i.txt",solve_id).c_str()).Write("M",system,q,s);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-b-%i.txt",solve_id).c_str()).Write("b",b);}

    CONJUGATE_GRADIENT<T> cg;
    x=b;
    cg.Solve(system,x,b,q,s,r,k,z,cg_tolerance,1,1000000);

    if(verbose){OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-x-%i.txt",solve_id).c_str()).Write("x",x);}

    for(int i=0;i<faces.m;i++) u(faces(i))=x.v(i);
}

template void Apply_Viscosity<double,VECTOR<double,1>,1>(GRID<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> >&,const BOUNDARY_CONDITIONS<VECTOR<double,1> >&,double,double,double,double,int,
    double,double,bool);
template void Apply_Viscosity<double,VECTOR<double,2>,2>(GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> >&,const BOUNDARY_CONDITIONS<VECTOR<double,2> >&,double,double,double,double,int,
    double,double,bool);

