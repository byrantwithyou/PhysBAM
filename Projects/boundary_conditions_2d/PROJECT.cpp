#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include "POISSON_PROJECTION_SYSTEM.h"
#include "PROJECT.h"
using namespace PhysBAM;


template<class TV,class T>
TV Centroid(const ARRAY<VECTOR<TV,TV::m> >& ar,T& A,TV& N)
{
    TV X;
    A=0;
    N=TV();
    for(int i=0;i<ar.m;i++){
        typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE face(ar(i));
        T a=face.Size();
        X+=a*face.Center();
        N+=a*face.Normal();
        A+=a;}
    N/=A;
    return X/A;
}

// phi at nodes
template<class T,class TV,class TV_INT>
void Project(const GRID<TV>& grid,int ghost,const ARRAY<T,TV_INT>& phi,boost::function<TV(TV X)> u_star,
    boost::function<TV(TV X)> u_projected,boost::function<T(TV X)> p,T density,T theta_threshold,T cg_tolerance,
    bool use_p_null_mode,bool use_bc,bool print_matrix)
{
    ARRAY<TV> u_loc,p_loc;
    ARRAY<T> us,u_proj,S,u_bc;
    ARRAY<FACE_INDEX<TV::m> > u_face;
    ARRAY<T,FACE_INDEX<TV::m> > us_grid(grid);
    HASHTABLE<TV_INT,int> cell_index;
    HASHTABLE<FACE_INDEX<TV::m>,int> face_index;
    INTERPOLATED_COLOR_MAP<T> color;
    color.Initialize_Colors(1e-12,10,true,true,true);
    POISSON_PROJECTION_SYSTEM<TV> system;
    system.neg_divergence.Reset(0);

    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
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
                    TV N,centroid=Centroid(boundary(2*a+s)(1),area,N);
                    if(!area) continue;
                    face_index.Insert(face,u_face.Append(face));
                    S.Append(area);
                    u_loc.Append(centroid+it.Location()-grid.dX/2);
//                    u_loc.Last()=grid.Face(face);
                    us.Append(u_star(u_loc.Last())(a));
                    Add_Debug_Particle(u_loc.Last(),VECTOR<T,3>(0,1,1));
                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,us.Last()*TV::Axis_Vector(a));
                    us_grid(face)=us.Last();
                    system.beta_inverse.Append(1/(density*grid.dX(a)*area));}
                int fi=-1;
                if(face_index.Get(face,fi)){
                    used_cell=true;
                    system.neg_divergence.Append_Entry_To_Current_Row(fi,-(s*2-1)*S(fi));}}
        if(used_cell){
            p_loc.Append(it.Location());
            system.neg_divergence.Finish_Row();

            T bc=0;
            if(use_bc)
                for(int i=0;i<surface.m;i++){
                    surface(i)*=grid.dX;
                    typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE face(surface(i));
                    TV C=face.Center()+it.Location()-grid.dX/2,u=u_projected(C),n=face.Normal();
                    T un=u.Dot(n);
                    bc+=un*face.Size();
                    Add_Debug_Particle(C,VECTOR<T,3>(1,0,1));
                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,un*n);}
            u_bc.Append(bc);}}
    system.neg_divergence.n=u_loc.m;
    system.neg_divergence.Sort_Entries();
    system.Initialize();

    Dump_Levelset(grid,phi,VECTOR<T,3>(1,1,0));
    Flush_Frame(us_grid,"disc");

    if(use_p_null_mode){
        system.projections.Append(ARRAY<T>());
        system.projections.Last().Resize(system.poisson.n);
        system.projections.Last().Fill(1/sqrt((T)system.poisson.n));}

    KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > x,b;
    x.v.Resize(system.poisson.n);
    b.v.Resize(system.poisson.n);
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    system.gradient.Transpose_Times(us,b.v);
    b.v-=u_bc;

    for(int i=0;i<b.v.m;i++) Add_Debug_Particle(p_loc(i),color(abs(b.v(i))));
    Flush_Frame<TV>("divergence");

    KRYLOV_SOLVER<T>::Ensure_Size(vectors,b,2);
    if(print_matrix){
        OCTAVE_OUTPUT<T> oo("matrices.txt");
        oo.Write("A",system.poisson);
        oo.Write("Mi",system.beta_inverse);
        oo.Write("G",system.gradient);
        oo.Write("ND",system.neg_divergence);
        oo.Write("b",b);}

    CONJUGATE_RESIDUAL<T> solver;
    solver.Solve(system,x,b,vectors,cg_tolerance,1,1000000);
    if(print_matrix) OCTAVE_OUTPUT<T>("matrices.txt",true).Write("x",x);

    T li=0,l1=0,ave=0;
    int cnt=0;
    for(int i=0;i<x.v.m;i++){
        cnt++;
        ave+=x.v(i)-p(p_loc(i));}
    ave/=cnt;

    for(int i=0;i<x.v.m;i++){
        T z=abs(x.v(i)-p(p_loc(i))-ave);
        li=max(li,z);
        l1+=z;
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
//        if((u_loc(i)-grid.Face(u_face(i))).Magnitude()>1e-8) continue;
        T z=abs(u_proj(i)-u_projected(u_loc(i))(u_face(i).axis));
        li=max(li,z);
        l1+=z;
        cnt++;
        u_error(u_face(i))=z;}

    printf("u %g %g\n", li, l1/cnt);

    Dump_Levelset(grid,phi,VECTOR<T,3>(1,1,0));
    Flush_Frame(u_error,"errors");
}
template void Project<float,VECTOR<float,2>,VECTOR<int,2> >(GRID<VECTOR<float,2> > const&,int,ARRAY<float,VECTOR<int,2> > const&,
    boost::function<VECTOR<float,2> (VECTOR<float,2>)>,boost::function<VECTOR<float,2> (VECTOR<float,2>)>,
    boost::function<float (VECTOR<float,2>)>,float,float,float,bool,bool,bool);
template void Project<float,VECTOR<float,3>,VECTOR<int,3> >(GRID<VECTOR<float,3> > const&,int,ARRAY<float,VECTOR<int,3> > const&,
    boost::function<VECTOR<float,3> (VECTOR<float,3>)>,boost::function<VECTOR<float,3> (VECTOR<float,3>)>,
    boost::function<float (VECTOR<float,3>)>,float,float,float,bool,bool,bool);
template void Project<double,VECTOR<double,2>,VECTOR<int,2> >(GRID<VECTOR<double,2> > const&,int,ARRAY<double,VECTOR<int,2> > const&,
    boost::function<VECTOR<double,2> (VECTOR<double,2>)>,boost::function<VECTOR<double,2> (VECTOR<double,2>)>,
    boost::function<double (VECTOR<double,2>)>,double,double,double,bool,bool,bool);
template void Project<double,VECTOR<double,3>,VECTOR<int,3> >(GRID<VECTOR<double,3> > const&,int,ARRAY<double,VECTOR<int,3> > const&,
    boost::function<VECTOR<double,3> (VECTOR<double,3>)>,boost::function<VECTOR<double,3> (VECTOR<double,3>)>,
    boost::function<double (VECTOR<double,3>)>,double,double,double,bool,bool,bool);
