//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Collisions_And_Interactions/LEVELSET_VOLUME_COLLISIONS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
using namespace HETERO_DIFF;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_VOLUME_COLLISIONS<TV>::
LEVELSET_VOLUME_COLLISIONS(DEFORMABLE_PARTICLES<TV>& particles,T stiffness)
    :COLLISION_FORCE<TV>(particles),stiffness(stiffness),pe(0)
{
    undeformed_phi.Resize(particles.X.m);
}
//#####################################################################
// Function Add_Mesh
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Add_Mesh(OBJECT& object,const IMPLICIT_OBJECT<TV>& implicit_surface)
{
    object.Initialize_Hierarchy();
    collision_bodies.Append(&object);
    ARRAY<int> unique_particles(object.mesh.elements.Flattened());
    unique_particles.Prune_Duplicates();
    for(int i=0;i<unique_particles.m;i++)
        undeformed_phi(unique_particles(i))=implicit_surface(particles.X(unique_particles(i)));
    undeformed_phi.Subset(object.Get_Boundary_Object().mesh.elements.Flattened()).Fill(0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_VOLUME_COLLISIONS<TV>::
~LEVELSET_VOLUME_COLLISIONS()
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Simplex_Intersection
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Simplex_Intersection(const VECTOR<TV,TV::m+1>& s,const ARRAY<HYPER_PLANE>& f,POLYTOPE& polytope)
{
    VECTOR<TV,max_pts> polytope_vertex(s);
    for(int j=0;j<TV::m+1;j++)
        polytope.poly(j)=(1<<(TV::m+1))-1-(1<<j);
    polytope.size=TV::m+1;
    for(int i=0;i<f.m;i++){
        const HYPER_PLANE& p=f(i);
        POLYTOPE in,out;
        VECTOR<TV,max_pts> in_vert,out_vert;
        VECTOR<T,max_pts> in_phi,out_phi;
        for(int j=0;j<polytope.size;j++){
            T phi=p.Signed_Distance(polytope_vertex(j));
            if(phi<=0){
                in_phi(in.size)=phi;
                in_vert(in.size)=polytope_vertex(j);
                in.poly(in.size++)=polytope.poly(j);}
            else{
                out_phi(out.size)=phi;
                out_vert(out.size)=polytope_vertex(j);
                out.poly(out.size++)=polytope.poly(j);}}
        polytope=in;
        polytope_vertex=in_vert;
        if(!in.size) return;
        for(int j=0;j<in.size;j++)
            for(int k=0;k<out.size;k++){
                int b=in.poly(j)&out.poly(k);
                if(count_bits(b)!=2) continue;
                T lambda=in_phi(j)/(in_phi(j)-out_phi(k));
                polytope_vertex(polytope.size)=in_vert(j)+lambda*(out_vert(k)-in_vert(j));
                polytope.poly(polytope.size++)=b|p.plane;}}
}
//#####################################################################
// Function Triangulate
//#####################################################################
template<int n> static void
Triangulate(int number_of_planes,const LEVELSET_VOLUME_COLLISIONS_POLYTOPE& polytope,ARRAY<VECTOR<int,n> >& triangulation)
{
    if(!polytope.size) return;
    int x=polytope.poly(0);
    for(int i=0;i<number_of_planes;i++){
        int plane=1<<i;
        if(x&plane) continue;
        LEVELSET_VOLUME_COLLISIONS_POLYTOPE p;
        ARRAY<VECTOR<int,n-1> > t;
        for(int j=0;j<polytope.size;j++)
            if(polytope.poly(j)&plane)
                p.poly(p.size++)=polytope.poly(j)&~plane;
        if(p.size==0) continue;
        Triangulate(number_of_planes,p,t);
        for(int j=0;j<t.m;j++){
            VECTOR<int,n-1>& s=t(j);
            for(int k=0;k<s.m;k++)
                s(k)|=plane;
            triangulation.Append(s.Append(x));}}
}
//#####################################################################
// Function Triangulate
//#####################################################################
template<> void
Triangulate<1>(int number_of_planes,const LEVELSET_VOLUME_COLLISIONS_POLYTOPE& polytope,ARRAY<VECTOR<int,1> >& triangulation)
{
    triangulation.Append(VECTOR<int,1>(0));
}
//#####################################################################
// Function Add_Plane_Edge_Intersection
//#####################################################################
template<class T,class TV_ARRAY> static void
Add_Plane_Edge_Intersection(const VECTOR<int,5>& nodes, const TV_ARRAY& X, VECTOR<T,3>& v, VECTOR<MATRIX<T,3>,5>& dv, SYMMETRIC_MATRIX<TENSOR<VECTOR<T,3> >,5>& ddv)
{
    auto p0=From_Var<5,0>(X(nodes(0)));
    auto p1=From_Var<5,1>(X(nodes(1)));
    auto p2=From_Var<5,2>(X(nodes(2)));
    auto e0=From_Var<5,3>(X(nodes(3)));
    auto e1=From_Var<5,4>(X(nodes(4)));
    auto normal=(p1-p0).Cross(p2-p0);
    auto d0=normal.Dot(e0-p0);
    auto d1=normal.Dot(e1-p0);
    auto lambda=d0/(d0-d1);
    auto z=e0+lambda*(e1-e0);
    v=z.x;
    for(int i=0;i<5;i++)
        Get(dv(i),z.dx.x,i);
    for(int i=0;i<5;i++)
        for(int j=i;j<5;j++)
            Get(ddv(i,j),z.ddx.x,i,j);
}
//#####################################################################
// Function Add_Plane_Edge_Intersection
//#####################################################################
template<class T,class TV_ARRAY> static void
Add_Plane_Edge_Intersection(const VECTOR<int,4>& nodes, const TV_ARRAY& X, VECTOR<T,2>& v, VECTOR<MATRIX<T,2>,4>& dv, SYMMETRIC_MATRIX<VECTOR<MATRIX<T,2>,2>,4>& ddv)
{
    auto p0=From_Var<4,0>(X(nodes(0)));
    auto p1=From_Var<4,1>(X(nodes(1)));
    auto e0=From_Var<4,2>(X(nodes(2)));
    auto e1=From_Var<4,3>(X(nodes(3)));
    auto normal=(p1-p0);
    auto d0=normal.Dot(e0-p0);
    auto d1=normal.Dot(e1-p0);
    auto z=e0+d0/(d0-d1+1e16)*(e1-e0);
    v=z.x;
    for(int i=0;i<4;i++)
        dv(i)=z.dx(i);
    for(int i=0;i<4;i++)
        for(int j=i;j<4;j++)
            ddv(i,j)=z.ddx(i,j);
}
//#####################################################################
// Function Integrate_Levelset
//#####################################################################
template<class T,class T_ARRAY> static void
Integrate_Levelset(const VECTOR<VECTOR<T,2>,6>& v,const T_ARRAY& undeformed_phi,T stiffness,T& pe,VECTOR<VECTOR<T,2>,6>& dpe, SYMMETRIC_MATRIX<MATRIX<T,2>,6>& ddpe)
{
    auto a=From_Var<6,0>(v(0));
    auto b=From_Var<6,1>(v(1));
    auto c=From_Var<6,2>(v(2));
    auto d=From_Var<6,3>(v(3));
    auto e=From_Var<6,4>(v(4));
    auto f=From_Var<6,5>(v(5));
    auto X=(d+e+f)/3-c;
    auto A=a-c;
    auto B=b-c;
    //TODO Finish 2D
    /* auto Det=A.Dot(B_cross_C); */
    /* auto Y=X/Det; */
    /* auto wA=Y.Dot(B_cross_C); */
    /* auto C_cross_A=C.Cross(A); */
    /* auto wB=Y.Dot(C_cross_A); */
    /* auto A_cross_B=A.Cross(B); */
    /* auto wC=Y.Dot(A_cross_B); */
    /* auto wD=1-(wA+wB+wC); */
    /* auto phi=undeformed_phi(0)*wA+undeformed_phi(1)*wB+undeformed_phi(2)*wC+undeformed_phi(3)*wD; */
    /* auto integral=abs((e-h).Dot((f-h).Cross(g-h))*phi); */
    /* pe+=integral.x; */
    /* for(int i=0;i<8;i++) */
    /*     dpe(i)=integral.dx(i); */
    /* for(int i=0;i<8;i++) */
    /*     for(int j=0;j<8;j++) */
    /*         ddpe(i,j)=integral.ddx(i,j); */
}
//#####################################################################
// Function Integrate_Levelset
//#####################################################################
template<class T,class T_ARRAY> static void
Integrate_Levelset(const VECTOR<VECTOR<T,3>,8>& v,const T_ARRAY& undeformed_phi,T stiffness,T& pe,VECTOR<VECTOR<T,3>,8>& dpe, SYMMETRIC_MATRIX<MATRIX<T,3>,8>& ddpe)
{
    auto a=From_Var<8,0>(v(0));
    auto b=From_Var<8,1>(v(1));
    auto c=From_Var<8,2>(v(2));
    auto d=From_Var<8,3>(v(3));
    auto e=From_Var<8,4>(v(4));
    auto f=From_Var<8,5>(v(5));
    auto g=From_Var<8,6>(v(6));
    auto h=From_Var<8,7>(v(7));
    auto X=(e+f+g+h)/4-d;
    auto A=a-d;
    auto B=b-d;
    auto C=c-d;
    auto B_cross_C=B.Cross(C);
    auto Det=A.Dot(B_cross_C);
    auto Y=X/Det;
    auto wA=Y.Dot(B_cross_C);
    auto C_cross_A=C.Cross(A);
    auto wB=Y.Dot(C_cross_A);
    auto A_cross_B=A.Cross(B);
    auto wC=Y.Dot(A_cross_B);
    auto wD=1-(wA+wB+wC);
    auto phi=undeformed_phi(0)*wA+undeformed_phi(1)*wB+undeformed_phi(2)*wC+undeformed_phi(3)*wD;
    auto integral=stiffness*abs((e-h).Dot((f-h).Cross(g-h))*phi);
    pe+=integral.x;
    if(integral.x<=0) return;
    for(int i=0;i<8;i++)
        dpe(i)=integral.dx(i);
    for(int i=0;i<8;i++)
        for(int j=0;j<8;j++)
            ddpe(i,j)=integral.ddx(i,j);
}
//#####################################################################
// Function Calculate_Vertex_Dependencies
//#####################################################################
static void
Calculate_Vertex_Dependencies(const VECTOR<int,4>& simplex,const VECTOR<int,8>& particle_indices,ARRAY<int>& vertex_dependencies,ARRAY<VECTOR<int,5> >& non_particle_vertex_nodes)
{
    const int n=4;
    int mask=(1<<n)-1;
    for(int j=0;j<n;j++)
        vertex_dependencies.Append(j);
    for(int j=0;j<n;j++){
        int vertex=simplex(j);
        int p1=vertex&mask;
        int p2=vertex>>n;
        switch(count_bits(vertex&mask)){
            case 0:{
                int node=n+integer_log_exact(p2^mask);
                vertex_dependencies.Append(node);
                break;}
            case 1:{
                VECTOR<int,3> p=VECTOR<int,n>(0,1,2,3).Remove_Index(integer_log_exact(p1));
                p2^=mask;
                int v0=p2&-p2;
                int v1=p2-v0;
                VECTOR<int,2> e(integer_log_exact(v0)+n,integer_log_exact(v1)+n);
                int match=-1;
                for(int i=0;i<3;i++)
                    for(int k=0;k<2;k++)
                        if(particle_indices(p(i))==particle_indices(e(k)))
                            match=p(i);
                if(match>=0) vertex_dependencies.Append(match);
                else{
                    non_particle_vertex_nodes.Append(p.Append_Elements(e));
                    vertex_dependencies.Append(-non_particle_vertex_nodes.m);}
                break;}
            case 2:{
                VECTOR<int,3> p=VECTOR<int,n>(4,5,6,7).Remove_Index(integer_log_exact(p2));
                p1^=mask;
                int v0=p1&-p1;
                int v1=p1-v0;
                VECTOR<int,2> e(integer_log_exact(v0),integer_log_exact(v1));
                int match=-1;
                for(int i=0;i<3;i++)
                    for(int k=0;k<2;k++)
                        if(particle_indices(p(i))==particle_indices(e(k)))
                            match=p(i);
                if(match>=0) vertex_dependencies.Append(match);
                else{
                    non_particle_vertex_nodes.Append(p.Append_Elements(e));
                    vertex_dependencies.Append(-non_particle_vertex_nodes.m);}
                break;}
            case 3:{
                int node=integer_log_exact(p1^mask);
                vertex_dependencies.Append(node);
                break;}
            default:
                PHYSBAM_ASSERT(count_bits(vertex)==n);
                break;
        }
    }
}
//#####################################################################
// Function Calculate_Vertex_Dependencies
//#####################################################################
#if 0
static void
Calculate_Vertex_Dependencies(const VECTOR<int,3>& simplex,ARRAY<int>& vertex_dependencies,ARRAY<VECTOR<int,4> >& non_particle_vertex_nodes)
{
    const int n=3;
    int mask=(1<<n)-1;
    for(int j=0;j<n;j++)
        vertex_dependencies.Append(j);
    for(int j=0;j<n;j++){
        int vertex=simplex(j);
        int p1=vertex&mask;
        int p2=vertex>>n;
        switch(count_bits(vertex&mask)){
            case 0:{
                int node=2*n-integer_log_exact(p2^mask);
                vertex_dependencies.Append(node);
                break;}
            case 1:{
                VECTOR<int,2> p=VECTOR<int,n>(0,1,2).Remove_Index(integer_log_exact(p1));
                VECTOR<int,2> e=VECTOR<int,n>(3,4,5).Remove_Index(integer_log_exact(p2));
                non_particle_vertex_nodes.Append(p.Append_Elements(e));
                vertex_dependencies.Append(-non_particle_vertex_nodes.m);
                break;}
            case 2:{
                int node=n-integer_log_exact(p1^mask);
                vertex_dependencies.Append(node);
                break;}
            default:
                PHYSBAM_ASSERT(count_bits(vertex)==n);
                break;
        }
    }
}
#endif
//#####################################################################
// Function Integrate_Simplex
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Integrate_Simplex(const VECTOR<int,TV::m+1>& simplex,const X_ARRAY& X,const PHI_ARRAY& nodewise_undeformed_phi,VECTOR<TV,2*TV::m+2>& df,MATRIX<MATRIX<T,TV::m>,2*TV::m+2>& ddf)
{
    ARRAY<int> vertex_dependencies;
    ARRAY<VECTOR<int,TV::m+2> > non_particle_vertex_nodes;
    Calculate_Vertex_Dependencies(simplex,X.indices,vertex_dependencies,non_particle_vertex_nodes);
    VECTOR<TV,2*TV::m+2> vertex;
    ARRAY<VECTOR<MATRIX<T,TV::m>,TV::m+2> > non_particle_vertex_jacobian;
    ARRAY<SYMMETRIC_MATRIX<TENSOR<VECTOR<T,TV::m> >,TV::m+2> > non_particle_vertex_tensor;
    for(int i=0;i<vertex_dependencies.m;i++){
        int k=vertex_dependencies(i);
        if(k<0){
            TV v;
            VECTOR<MATRIX<T,TV::m>,TV::m+2> dv;
            SYMMETRIC_MATRIX<TENSOR<VECTOR<T,TV::m> >,TV::m+2> ddv;
            Add_Plane_Edge_Intersection(non_particle_vertex_nodes(~k),X,v,dv,ddv);
            vertex(i)=v;
            non_particle_vertex_jacobian.Append(dv);
            non_particle_vertex_tensor.Append(ddv);}
        else
            vertex(i)=X(k);}
    VECTOR<TV,2*TV::m+2> dintegral;
    SYMMETRIC_MATRIX<MATRIX<T,TV::m>,2*TV::m+2> ddintegral;
    Integrate_Levelset(vertex,nodewise_undeformed_phi,stiffness,pe,dintegral,ddintegral);
    for(int i=0;i<vertex_dependencies.m;i++){
        int vdi=vertex_dependencies(i);
        if(vdi>=0){
            df(vdi)+=dintegral(i);
            for(int j=i;j<vertex_dependencies.m;j++){
                int vdj=vertex_dependencies(j);
                if(vdj>=0)
                    ddf(vdi,vdj)+=ddintegral(i,j);
                else{
                    const VECTOR<int,TV::m+2>& npj=non_particle_vertex_nodes(~vdj);
                    for(int m=0;m<TV::m+2;m++)
                        ddf(vdi,npj(m))+=ddintegral(i,j)*non_particle_vertex_jacobian(~vdj)(m);}}}
        else{
            const VECTOR<int,TV::m+2>& npi=non_particle_vertex_nodes(~vdi);
            for(int l=0;l<TV::m+2;l++){
                df(npi(l))+=non_particle_vertex_jacobian(~vdi)(l).Transpose_Times(dintegral(i));
                for(int j=0;j<vertex_dependencies.m;j++){
                    int vdj=vertex_dependencies(j);
                    if(vdj>=0)
                        ddf(npi(l),vdj)+=ddintegral(i,j)*non_particle_vertex_jacobian(~vdi)(l);
                    else{
                        const VECTOR<int,TV::m+2>& npj=non_particle_vertex_nodes(~vdj);
                        for(int m=0;m<TV::m+2;m++){
                            ddf(npi(l),npj(m))+=ddintegral(i,j)*non_particle_vertex_jacobian(~vdi)(l)*non_particle_vertex_jacobian(~vdj)(m);
                            ddf(npi(l),npi(m))+=Contract_0(non_particle_vertex_tensor(~vdi)(l,m),dintegral(i));}}}}}}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    pe=0;
    grad_pe.Remove_All();
    H_pe.Remove_All();
    overlapping_particles.Remove_All();
    for(int i=0;i<collision_bodies.m;i++)
        collision_bodies(i)->hierarchy->Update_Boxes(1e-15);
    for(int i=0;i<collision_bodies.m;i++)
        for(int j=0;j<collision_bodies.m;j++)
            Update_Position_Based_State_Pair(*collision_bodies(i),*collision_bodies(j));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Update_Position_Based_State_Pair(const OBJECT& o0,const OBJECT& o1)
{
    for(int e=0;e<o0.mesh.elements.m;e++){
        const VECTOR<int,TV::m+1>& nodes0=o0.mesh.elements(e);
        ARRAY<int> intersection_list;
        VECTOR<TV,TV::m+1> s(particles.X.Subset(nodes0));
        o1.hierarchy->Intersection_List(RANGE<TV>::Bounding_Box(s),intersection_list);
        for(int n=0;n<intersection_list.m;n++){
            int t=intersection_list(n);
            if(t==e && &o0==&o1) continue;
            const VECTOR<int,TV::m+1>& nodes1=o1.mesh.elements(t);
            int num_match=0;
            for(int i=0;i<TV::m+1;i++)
                if(nodes1.Contains(nodes0(i)))
                    num_match++;
            if(num_match==TV::m) continue;
            POLYTOPE polytope;
            VECTOR<TV,TV::m+1> v(particles.X.Subset(nodes1));
            ARRAY<HYPER_PLANE> f;
            for(int i=0;i<TV::m+1;i++){
                int p=1<<(TV::m+1+i);
                VECTOR<TV,TV::m> w=v.Remove_Index(i);
                TV n=TRIANGLE_3D<T>::Normal(w);
                T s=w(0).Dot(n);
                if(v(i).Dot(n)>s){n=-n;s=-s;}
                f.Append({p,n,s});}
            Simplex_Intersection(s,f,polytope);
            if(polytope.size<=TV::m)
                continue;
            ARRAY<VECTOR<int,TV::m+1> > triangulation;
            Triangulate(2*(TV::m+1),polytope,triangulation);
            VECTOR<TV,2*TV::m+2> df;
            MATRIX<MATRIX<T,TV::m>,2*TV::m+2> ddf;
            VECTOR<int,2*TV::m+2> nodes=nodes0.Append_Elements(nodes1);
            for(int i=0;i<triangulation.m;i++)
                Integrate_Simplex(triangulation(i),particles.X.Subset(nodes),undeformed_phi.Subset(nodes0),df,ddf);
            overlapping_particles.Append(nodes);
            grad_pe.Append(df);
            H_pe.Append(ddf);}}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int pp=0;pp<overlapping_particles.m;pp++)
        F.Subset(overlapping_particles(pp))-=grad_pe(pp);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int LEVELSET_VOLUME_COLLISIONS<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const
{
    for(int pp=0;pp<overlapping_particles.m;pp++){
        const VECTOR<int,2*TV::m+2>& n=overlapping_particles(pp);
        for(int i=0;i<n.m;i++){
            int p=n(i);
            for(int j=0;j<n.m;j++)
                F(p)-=H_pe(pp)(i,j)*V(n(j))*scale;}}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_VOLUME_COLLISIONS<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_VOLUME_COLLISIONS<TV>::
Potential_Energy(const T time) const
{
    return pe;
}
//#####################################################################
// Function Apply_Friction
//#####################################################################
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Apply_Friction(ARRAY_VIEW<TV> V,const T time) const
{
}
/* template class LEVELSET_VOLUME_COLLISIONS<VECTOR<float,2> >; */
/* template class LEVELSET_VOLUME_COLLISIONS<VECTOR<double,2> >; */
template class LEVELSET_VOLUME_COLLISIONS<VECTOR<float,3> >;
template class LEVELSET_VOLUME_COLLISIONS<VECTOR<double,3> >;
}
