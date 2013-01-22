#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include "POISSON_PROJECTION_SYSTEM.h"
#include <boost/function.hpp>

using namespace PhysBAM;

template<class T,class TV>
void Project(const GRID<TV>& grid,int ghost,boost::function<T(TV X)> phi,boost::function<TV(TV X)> u_star,
    boost::function<TV(TV X)> u_projected,boost::function<T(TV X)> p,T density,T theta_threshold,T cg_tolerance,bool use_p_null_mode)
{
    typedef VECTOR<int,TV::m> TV_INT;
    ARRAY<TV> u_loc,p_loc;
    ARRAY<T> us,u_proj;
    ARRAY<FACE_INDEX<TV::m> > u_face;
    ARRAY<T> S,R;
    SPARSE_MATRIX_FLAT_MXN<T> G_hat;
    G_hat.Reset(0);
    HASHTABLE<TV_INT,int> cell_index;
    int next_cell=0;
    INTERPOLATED_COLOR_MAP<T> color;
    color.Initialize_Colors(1e-12,1,true,true,true);

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face=it.Full_Index();
        TV X=grid.X(face.Cell_Index(0)),Y=grid.X(face.Cell_Index(1));
        TV M=(X+Y)/2,D=(Y-X)/2,A=M+D.Orthogonal_Vector(),B=M-D.Orthogonal_Vector();
        T pa=phi(A),pb=phi(B);
        if(abs(pa)<theta_threshold) pa=sign_nonzero(pa)*theta_threshold;
        if(abs(pb)<theta_threshold) pb=sign_nonzero(pb)*theta_threshold;
        if(pa>0){
            if(pb>0) continue;
            T th=pb/(pb-pa);
            S.Append(th);
            u_loc.Append(B+th/2*(A-B));}
        else if(pb<=0){
            S.Append(1);
            u_loc.Append(M);}
        else{
            T th=pa/(pa-pb);
            S.Append(th);
            u_loc.Append(A+th/2*(B-A));}
        u_face.Append(face);
        R.Append(density);
        us.Append(u_star(u_loc.Last())(it.Axis()));

        for(int s=0;s<2;s++){
            int index=0;
            if(!cell_index.Get(face.Cell_Index(s),index)){
                p_loc.Append(grid.Center(face.Cell_Index(s)));
                cell_index.Set(face.Cell_Index(s),index=next_cell++);}
            G_hat.Append_Entry_To_Current_Row(index,(s*2-1)*grid.one_over_dX(it.Axis()));}
        G_hat.Finish_Row();}
    G_hat.n=next_cell;
    G_hat.Sort_Entries();
    S*=grid.dX.Product();

    LOG::cout<<S<<std::endl;
    LOG::cout<<G_hat<<std::endl;

    POISSON_PROJECTION_SYSTEM<TV> system;
    system.gradient=G_hat;
    system.gradient.Set_Times_Diagonal(S);
    system.beta_inverse=S/R;
    system.Initialize();
    if(use_p_null_mode){
        system.projections.Append(ARRAY<T>());
        system.projections.Last().Resize(system.poisson.n);
        system.projections.Last().Fill(1/sqrt((T)system.poisson.n));}

    KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > x,b;
    x.v.Resize(system.poisson.n);
    b.v.Resize(system.poisson.n);
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    system.gradient.Transpose_Times(us,b.v);
    b.v*=-(T)1;

    CONJUGATE_GRADIENT<T> solver;
    solver.Solve(system,x,b,vectors,cg_tolerance,1,1000000);

    T li=0,l1=0;
    int cnt=0;
    for(int i=0;i<x.v.m;i++){
        T z=abs(x.v(i)-p(p_loc(i)));
        li=max(li,z);
        l1+=z;
        cnt++;
        Add_Debug_Particle(p_loc(i),color(z));}
    printf("p %g %g\n", li, l1/cnt);

    li=0;
    l1=0;
    cnt=0;
    u_proj.Resize(system.gradient.m);
    system.gradient.Times(x.v,u_proj);
    u_proj=us-system.beta_inverse*u_proj;
    ARRAY<T,FACE_INDEX<TV::m> > u_error(grid,3);
    for(int i=0;i<u_proj.m;i++){
        T z=abs(u_proj(i)-u_projected(u_loc(i))(u_face(i).axis));
        li=max(li,z);
        l1+=z;
        cnt++;
        u_error(u_face(i))=z;}

    printf("u %g %g\n", li, l1/cnt);

    Flush_Frame(u_error,"errors");
}

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    PARSE_ARGS parse_args(argc,argv);
    int refine=1,resolution=32,interface=0,velocity_field=0;
    T rho=1,kg=1,m=1,s=1;
    bool test_analytic_diff=false,use_p_null_mode=false,dump_matrix=false;
    parse_args.Extra(&interface,"number","interface to use");
    parse_args.Extra(&velocity_field,"number","velocity to use");
    parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
    parse_args.Add("-dump_matrix",&dump_matrix,"dump out system and rhs");
    parse_args.Add("-rho",&rho,"density","density for first fluid region");
    parse_args.Add("-m",&m,"scale","meter scale");
    parse_args.Add("-s",&s,"scale","second scale");
    parse_args.Add("-kg",&kg,"scale","kilogram scale");
    parse_args.Add("-test_diff",&test_analytic_diff,"test analytic derivatives");
    parse_args.Add("-refine",&refine,"num","Refine space/time by this factor");
    parse_args.Add("-null_p",&use_p_null_mode,"Assume pressure null mode and project it out");
    parse_args.Parse();

    GRID<TV> grid(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),grid,"output");

    boost::function<T(TV X)> phi,p;
    boost::function<TV(TV X)> u_star,u_projected;

    switch(interface){
        case 0:phi=[](TV X){return (X-.5).Magnitude()-.3;};break;
        default: PHYSBAM_FATAL_ERROR("Unrecognized interface");}

    switch(velocity_field){
        case 0:
            u_star=[](TV X){return TV(0,0);};
            u_projected=[](TV X){return TV(0,0);};
            p=[](TV X){return 0;};
            break;
        case 1:
            u_star=[](TV X){return TV(2,1);};
            u_projected=[](TV X){return TV(2,1);};
            p=[](TV X){return 0;};
            break;
        default: PHYSBAM_FATAL_ERROR("Unrecognized velocity");}

    Project<T,TV>(grid,3,phi,u_star,u_projected,p,rho,(T)1e-8,(T)1e-12,use_p_null_mode);
    Flush_Frame<TV>("flush");
    return 0;
}


