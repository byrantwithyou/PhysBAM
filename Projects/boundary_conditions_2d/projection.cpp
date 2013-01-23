#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include "POISSON_PROJECTION_SYSTEM.h"
#include <boost/function.hpp>

using namespace PhysBAM;

template<class TV,class T>
TV Centroid(const ARRAY<VECTOR<TV,TV::m> >& ar,T* pA,TV* N=0)
{
    TV X;
    T A=0;
    if(N) *N=TV();
    for(int i=0;i<ar.m;i++){
        typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE face(ar(i));
        T a=face.Size();
        X+=a*face.Center();
        if(N) *N+=a*face.Normal();
        A+=a;}
    if(pA) *pA=A;
    if(N) *N/=A;
    return X/A;
}

// phi at nodes
template<class T,class TV,class TV_INT>
void Project(const GRID<TV>& grid,int ghost,const ARRAY<T,TV_INT>& phi,boost::function<TV(TV X)> u_star,
    boost::function<TV(TV X)> u_projected,boost::function<T(TV X)> p,T density,T theta_threshold,T cg_tolerance,bool use_p_null_mode)
{
    ARRAY<TV> u_loc,p_loc;
    ARRAY<T> us,u_proj,S,R,u_bc;
    ARRAY<FACE_INDEX<TV::m> > u_face;
    ARRAY<T,FACE_INDEX<TV::m> > us_grid(grid);
    SPARSE_MATRIX_FLAT_MXN<T> neg_div,G_hat;
    neg_div.Reset(0);
    HASHTABLE<TV_INT,int> cell_index;
    HASHTABLE<FACE_INDEX<TV::m>,int> face_index;
    INTERPOLATED_COLOR_MAP<T> color;
    color.Initialize_Colors(1e-12,10,true,true,true);

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        ARRAY<VECTOR<TV,TV::m> > surface;
        VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >,2>,2*TV::m> boundary;
        VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >*,2>,2*TV::m> pboundary;
        for(int f=0;f<2*TV::m;f++)
            for(int s=1;s<2;s++) // Inside only
                pboundary(f)(s)=&boundary(f)(s);
        MARCHING_CUBES<TV>::Get_Elements_For_Cell(surface,pboundary,phi,it.index);
        bool used_cell=false;
        for(int a=0;a<TV::m;a++)
            for(int s=0;s<2;s++){
                FACE_INDEX<TV::m> face(a,it.index);
                face.index(a)+=s;
                if(boundary(2*a+s)(1).m && !face_index.Contains(face)){
                    for(int t=0;t<boundary(2*a+s)(1).m;t++)
                        boundary(2*a+s)(1)(t)*=grid.dX;
                    T area=0;
                    TV centroid=Centroid(boundary(2*a+s)(1),&area);
                    if(!area) continue;
                    face_index.Insert(face,u_face.Append(face));
                    S.Append(area);
                    u_loc.Append(centroid+it.Location()-grid.dX/2);
                    us.Append(u_star(u_loc.Last())(a));
                    Add_Debug_Particle(u_loc.Last(),VECTOR<T,3>(0,1,1));
                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,us.Last()*TV::Axis_Vector(a));
                    us_grid(face)=us.Last();
                    R.Append(density);}
                int fi=-1;
                if(face_index.Get(face,fi)){
                    used_cell=true;
                    neg_div.Append_Entry_To_Current_Row(fi,-(s*2-1)*grid.one_over_dX(a));}}
        if(used_cell){
            p_loc.Append(it.Location());
            neg_div.Finish_Row();

            if(surface.m){
                for(int t=0;t<surface.m;t++)
                    surface(t)*=grid.dX;
                T area=0;
                TV N,X=Centroid(surface,&area,&N)+it.Location()-grid.dX/2;
                T un=u_projected(X).Dot(N);
                u_bc.Append(un*area);
                Add_Debug_Particle(X,VECTOR<T,3>(1,0,1));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,un*N);}
            else u_bc.Append(0);}}
    neg_div.n=u_loc.m;
    neg_div.Sort_Entries();
    neg_div.Transpose(G_hat);

    Flush_Frame(us_grid,"disc");

    POISSON_PROJECTION_SYSTEM<TV> system;
    system.gradient=G_hat;
    system.gradient.Set_Diagonal_Times(S);
    system.beta_inverse.Resize(S.m);
    for(int i=0;i<system.beta_inverse.m;i++) system.beta_inverse(i)=1/(S(i)*R(i));
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

    for(int i=0;i<b.v.m;i++){
        Add_Debug_Particle(p_loc(i),color(abs(b.v(i))));
        Add_Debug_Particle(p_loc(i)+grid.dX/4,color(abs(u_bc(i))));
    }
    Flush_Frame<TV>("divergence");

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

    Dump_Levelset(grid,phi,true,VECTOR<T,3>(1,1,0));
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

    ARRAY<T,TV_INT> node_phi(grid.Node_Indices());
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next())
        node_phi(it.index)=phi(it.Location());

    Project<T,TV>(grid,3,node_phi,u_star,u_projected,p,rho,(T)1e-8,(T)1e-12,use_p_null_mode);
    Flush_Frame<TV>("flush");
    return 0;
}


