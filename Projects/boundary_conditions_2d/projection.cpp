#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <boost/function.hpp>

using namespace PhysBAM;

template<class T,class TV>
void Project(const GRID<TV>& grid,int ghost,boost::function<T(TV X)> phi,boost::function<TV(TV X)> u_star,boost::function<TV(TV X)> u_projected,boost::function<T(TV X)> p,T density,T theta_threshold,T cg_tolerance,bool use_p_null_mode)
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
    INTERPOLATED_COLOR_MAP<T> color(1e-12,1,true,true,true);

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

    for(int i=0;i<x.v.m;i++)
        Add_Debug_Particle(p_loc(i),color(abs(x.v-p(p_loc(i)))));

    system.gradient.Times(u_proj,x.v);
    u_proj=us-system.beta_inverse*u_proj;
    ARRAY<T,FACE_INDEX<TV::m> > u_error(grid,3);
    for(int i=0;i<u_proj.m;i++)
        u_error(u_face(i))=abs(u_proj(i)-u_projected(u_loc(i)));
    Flush_Frame(u_error,"");
}

int main(int argc,char* argv[])
{
    typedef double T;
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    PARSE_ARGS parse_args(argc,argv);
    int test_number=0,refine=1,resolution=32;
    T rho=1,kg=1,m=1,s=1;
    bool test_analytic_diff=false,use_p_null_mode=false,dump_matrix=false;
    parse_args.Extra(&test_number,"example number","example number to run");
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
    Project<T,TV>(grid,3,[](TV X){return (X-.5).Magnitude()-.3;},[](TV X){return TV(2,1);},[](TV X){return TV(2,1);},[](TV X){return 0;},rho,(T)1e-8,(T)1e-12,use_p_null_mode);


}













