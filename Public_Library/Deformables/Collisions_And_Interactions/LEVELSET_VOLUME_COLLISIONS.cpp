//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Collisions_And_Interactions/LEVELSET_VOLUME_COLLISIONS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_VOLUME_COLLISIONS<TV>::
    LEVELSET_VOLUME_COLLISIONS(DEFORMABLE_PARTICLES<TV>& particles,OBJECT& collision_body,IMPLICIT_OBJECT<TV>& implicit_surface,
        T stiffness)
    :COLLISION_FORCE<TV>(particles),collision_body(collision_body),stiffness(stiffness),pe(0)
{
    collision_body.Initialize_Hierarchy();
    undeformed_phi.Resize(particles.X.m);
    for(int i=0;i<particles.X.m;i++)
        undeformed_phi(i)=implicit_surface(particles.X(i));
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
template<class TV> void LEVELSET_VOLUME_COLLISIONS<TV>::
Simplex_Intersection(const SIMPLEX& s,const ARRAY<HYPER_PLANE>& f,ARRAY<unsigned char>& polytope)
{
    ARRAY<TV> polytope_vertex(TV::m+1,false);
    for(int j=0;j<TV::m+1;j++){
        polytope(j)=(1u<<(TV::m+1))-1u-(1u<<(TV::m+1-j));
        polytope_vertex(j)=s.X(j);}
    ARRAY<T> d(TV::m+1);
    for(int i=0;i<f.m;i++){
        const HYPER_PLANE& p=f(i);
        d.Resize(polytope.m);
        for(int j=0;j<polytope_vertex.m-1;j++)
            d(j)=p.Signed_Distance(polytope_vertex(j));
        for(int j=polytope_vertex.m-1;j>=0;j--){
            T dj = d(j);
            if(dj<=0) continue;
            TV u=polytope_vertex(j);
            unsigned char planes=polytope(j);
            polytope_vertex.Remove_Index_Lazy(j);
            polytope.Remove_Index_Lazy(j);
            d.Remove_Index_Lazy(j);
            for(int k=polytope_vertex.m-1;k>=0;k--){
                if(d(k)>=0) continue;
                unsigned char edge=planes&polytope(k);
                if(count_bits(edge)!=2) continue;
                const TV& v=polytope_vertex(k);
                polytope_vertex.Append((d(k)*v-dj*u)/(d(k)-dj));
                polytope.Append((1u<<(TV::m+1+i))|edge);}}}
}

template<class T> VECTOR<T,3>
Normal_To_Face(int i,const INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<T,3>,int>,VECTOR<int,4>&>& v)
{
    VECTOR<int,3> nodes=VECTOR<int,4>(3,2,1,0).Remove_Index(i);
    VECTOR<T,3> n=(v(nodes(2))-v(nodes(0))).Cross(v(nodes(1))-v(nodes(0)));
    if(n.Dot(v(3-i)-v(nodes(0)))>0)
        n*=-1;
    return n.Normalized();
}

template<class T> VECTOR<T,2>
Normal_To_Face(int i,const INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<T,2>,int>,VECTOR<int,3>&>& v)
{
    VECTOR<int,2> nodes=VECTOR<int,3>(2,1,0).Remove_Index(i);
    VECTOR<T,2> n=(v(nodes(1))-v(nodes(0))).Perpendicular();
    if(n.Dot(v(2-i)-v(nodes(0)))>0)
        n*=-1;
    return n.Normalized();
}

template<int n> void 
Triangulate(int number_of_planes,const ARRAY<unsigned char>& polytope,ARRAY<VECTOR<unsigned char,n> >& triangulation)
{
    unsigned char x=polytope(0);
    for(int i=0;i<number_of_planes;i++){
        unsigned char plane=1u<<i;
        if(x&plane) continue;
        ARRAY<unsigned char> p;
        ARRAY<VECTOR<unsigned char,n-1> > t;
        for(int j=0;j<polytope.m;j++)
            if(polytope(j)&plane)
                p.Append(polytope(j)&~plane);
        Triangulate(number_of_planes,p,t);
        for(int j=0;j<t.m;j++)
            triangulation.Append(t(j).Append(x));}
}

template<> void 
Triangulate<1>(int number_of_planes,const ARRAY<unsigned char>& polytope,ARRAY<VECTOR<unsigned char,1> >& triangulation)
{
        for(int i=0;i<polytope.m;i++)
            triangulation.Append(VECTOR<unsigned char,1>(polytope(i)));
        return;
}

template<class T,class TV_ARRAY> void 
Add_Plane_Edge_Intersection(const VECTOR<int,5> nodes, const TV_ARRAY& X, VECTOR<T,3>& v, VECTOR<MATRIX<T,3>,5>& dv, SYMMETRIC_MATRIX<VECTOR<MATRIX<T,3>,3>,5>& ddv)
{
    auto p0=From_Var<5,0>(X(nodes(0)));
    auto p1=From_Var<5,1>(X(nodes(1)));
    auto p2=From_Var<5,2>(X(nodes(2)));
    auto e0=From_Var<5,3>(X(nodes(3)));
    auto e1=From_Var<5,4>(X(nodes(4)));
    auto normal=(p1-p0).Cross(p2-p0);
    auto d0=normal.Dot(e0-p0);
    auto d1=normal.Dot(e1-p0);
    auto z=(e0*d1-e1*d0)/(d1-d0);
    for(int i=0;i<5;i++)
        dv(i)=z.dx(i);
    for(int i=0;i<5;i++)
        for(int j=i;j<5;j++)
            ddv(i,j)=z.ddx(i,j);
}

template<class T,class TV_ARRAY> void 
Add_Plane_Edge_Intersection(const VECTOR<int,4> nodes, const TV_ARRAY& X, VECTOR<T,2>& v, VECTOR<MATRIX<T,2>,4>& dv, SYMMETRIC_MATRIX<VECTOR<MATRIX<T,2>,2>,4>& ddv)
{
    auto p0=From_Var<4,0>(X(nodes(0)));
    auto p1=From_Var<4,1>(X(nodes(1)));
    auto e0=From_Var<4,2>(X(nodes(2)));
    auto e1=From_Var<4,3>(X(nodes(3)));
    auto normal=(p1-p0);
    auto d0=normal.Dot(e0-p0);
    auto d1=normal.Dot(e1-p0);
    auto z=(e0*d1-e1*d0)/(d1-d0);
    for(int i=0;i<4;i++)
        dv(i)=z.dx(i);
    for(int i=0;i<4;i++)
        for(int j=i;j<4;j++)
            ddv(i,j)=z.ddx(i,j);
}

template<class T,class T_ARRAY>
void Integrate_Levelset(const VECTOR<VECTOR<T,2>,6>& v,const T_ARRAY& undeformed_phi,T& pe,VECTOR<VECTOR<T,2>,6>& dpe, SYMMETRIC_MATRIX<MATRIX<T,2>,6>& ddpe)
{
    auto a=From_Var<8,0>(v(0));
    auto b=From_Var<8,1>(v(1));
    auto c=From_Var<8,2>(v(2));
    auto d=From_Var<8,3>(v(3));
    auto e=From_Var<8,4>(v(4));
    auto f=From_Var<8,5>(v(5));
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

template<class T,class T_ARRAY>
void Integrate_Levelset(const VECTOR<VECTOR<T,3>,8>& v,const T_ARRAY& undeformed_phi,T& pe,VECTOR<VECTOR<T,3>,8>& dpe, SYMMETRIC_MATRIX<MATRIX<T,3>,8>& ddpe)
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
    auto integral=abs((e-h).Dot((f-h).Cross(g-h))*phi);
    pe+=integral.x;
    for(int i=0;i<8;i++)
        dpe(i)=integral.dx(i);
    for(int i=0;i<8;i++)
        for(int j=0;j<8;j++)
            ddpe(i,j)=integral.ddx(i,j);
}

void Calculate_Vertex_Dependencies(VECTOR<unsigned char,4> simplex,ARRAY<int> vertex_dependencies,ARRAY<VECTOR<int,5> > non_particle_vertex_nodes)
{
    const int n=4;
    unsigned char mask=(1u<<n)-1u;
    for(int j=0;j<n;j++)
        vertex_dependencies.Append(j);
    for(int j=0;j<n;j++){
        unsigned char vertex=simplex(j);
        unsigned char p1=vertex&mask;
        unsigned char p2=vertex>>n;
        switch(count_bits(vertex&mask)){
            case 0:{
                int node=2*n-integer_log_exact(p2^mask);
                vertex_dependencies.Append(node);
                break;}
            case 1:{
                VECTOR<int,3>  p=VECTOR<int,n>(3,2,1,0).Remove_Index(integer_log_exact(p1));
                unsigned char v0=p2&(unsigned char)-(signed char)p2;
                unsigned char v1=p2-v0;
                VECTOR<int,2> e(integer_log_exact(v0)+n,integer_log_exact(v1)+n);
                non_particle_vertex_nodes.Append(p.Append_Elements(e));
                vertex_dependencies.Append(-non_particle_vertex_nodes.m);
                break;}
            case 2:{
                VECTOR<int,3>  p=VECTOR<int,n>(7,6,5,4).Remove_Index(integer_log_exact(p2));
                unsigned char v0=p1&(unsigned char)-(signed char)p1;
                unsigned char v1=p1-v0;
                VECTOR<int,2> e(integer_log_exact(v0),integer_log_exact(v1));
                non_particle_vertex_nodes.Append(p.Append_Elements(e));
                vertex_dependencies.Append(-non_particle_vertex_nodes.m);
                break;}
            case 3:{
                int node=n-integer_log_exact(p1^mask);
                vertex_dependencies.Append(node);
                break;}
            default:
                PHYSBAM_ASSERT(count_bits(vertex)==n);
                break;
        }
    }
}

void Calculate_Vertex_Dependencies(VECTOR<unsigned char,3> simplex,ARRAY<int> vertex_dependencies,ARRAY<VECTOR<int,4> > non_particle_vertex_nodes)
{
    const int n=3;
    unsigned char mask=(1u<<n)-1u;
    for(int j=0;j<n;j++)
        vertex_dependencies.Append(j);
    for(int j=0;j<n;j++){
        unsigned char vertex=simplex(j);
        unsigned char p1=vertex&mask;
        unsigned char p2=vertex>>n;
        switch(count_bits(vertex&mask)){
            case 0:{
                int node=2*n-integer_log_exact(p2^mask);
                vertex_dependencies.Append(node);
                break;}
            case 1:{
                VECTOR<int,2>  p=VECTOR<int,n>(2,1,0).Remove_Index(integer_log_exact(p1));
                VECTOR<int,2>  e=VECTOR<int,n>(5,4,3).Remove_Index(integer_log_exact(p2));
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

template<class TV> void
LEVELSET_VOLUME_COLLISIONS<TV>::Integrate_Simplex(VECTOR<unsigned char,TV::m+1> simplex,const X_ARRAY& X,const PHI_ARRAY& nodewise_undeformed_phi,VECTOR<TV,2*TV::m+2>& df,MATRIX<MATRIX<T,TV::m>,2*TV::m+2>& ddf)
{
    ARRAY<int> vertex_dependencies;
    ARRAY<VECTOR<int,TV::m+2> > non_particle_vertex_nodes;
    VECTOR<TV,2*TV::m+2> vertex;
    ARRAY<VECTOR<MATRIX<T,TV::m>,TV::m+2> > non_particle_vertex_jacobian;
    ARRAY<SYMMETRIC_MATRIX<VECTOR<MATRIX<T,TV::m>,TV::m>,TV::m+2> > non_particle_vertex_tensor;
    for(int i=0;i<vertex_dependencies.m;i++){
        int k=vertex_dependencies(i);
        if(k<0){
            TV v;
            VECTOR<MATRIX<T,TV::m>,TV::m+2> dv;
            SYMMETRIC_MATRIX<VECTOR<MATRIX<T,TV::m>,TV::m>,TV::m+2> ddv;
            Add_Plane_Edge_Intersection(non_particle_vertex_nodes(~k),X,v,dv,ddv);    
            vertex(i)=v;
            non_particle_vertex_jacobian.Append(dv);
            non_particle_vertex_tensor.Append(ddv);}
        else
            vertex(i)=X(k);}
    VECTOR<TV,2*TV::m+2> dintegral;
    SYMMETRIC_MATRIX<MATRIX<T,TV::m>,2*TV::m+2> ddintegral;
    Integrate_Levelset(vertex,nodewise_undeformed_phi,pe,dintegral,ddintegral);
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
                        ddf(vdi,npj(m))+=ddintegral(i,j)*non_particle_vertex_jacobian(j)(npj(m));}}}
        else{
            const VECTOR<int,TV::m+2>& npi=non_particle_vertex_nodes(~vdi);
            for(int l=0;l<TV::m+2;l++){
                df(vdi)+=non_particle_vertex_jacobian(i)(npi(l)).Transpose_Times(dintegral(npi(l)));
                for(int j=0;j<vertex_dependencies.m;j++){
                    int vdj=vertex_dependencies(j);
                    if(vdj>=0)
                        ddf(npi(l),vdj)+=ddintegral(i,j)*non_particle_vertex_jacobian(i)(npi(l));
                    else{
                        const VECTOR<int,TV::m+2>& npj=non_particle_vertex_nodes(~vdj);
                        for(int m=0;m<TV::m+2;m++){
                            ddf(npi(l),npj(m))+=ddintegral(i,j)*non_particle_vertex_jacobian(i)(npi(l))*non_particle_vertex_jacobian(j)(npj(m));
                            const VECTOR<MATRIX<T,TV::m>,TV::m> t=non_particle_vertex_tensor(i)(npi(l),npi(m));
                            for(int p=0;p<TV::m;p++) ddf(npi(l),npi(m))+=t(p)*dintegral(i)(p);}}}}}}
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
    collision_body.hierarchy->Update_Boxes(1e-15);
    for(int e=0;e<=collision_body.mesh.elements.m;e++){
        const VECTOR<int,TV::m+1>& nodes1=collision_body.mesh.elements(e);
        ARRAY<int> intersection_list;
        SIMPLEX s(particles.X.Subset(nodes1));
        collision_body.hierarchy->Intersection_List(s.Bounding_Box(),intersection_list);
        for(int n=0;n<intersection_list.m;n++){
            int t=intersection_list(n);
            if(t==e) continue;
            VECTOR<int,TV::m+1> nodes2=collision_body.mesh.elements(t);
            ARRAY<unsigned char> polytope(TV::m+1,false);
            auto v=particles.X.Subset(nodes2);
            ARRAY<HYPER_PLANE> f;
            for(int i=0;i<TV::m+1;i++){
                unsigned char p=1u<<(TV::m+1+i);
                TV n=Normal_To_Face(i,v);
                T s=v(0).Dot(n);
                f.Append({p,n,s});}
            Simplex_Intersection(s,f,polytope);
            if(polytope.m==0)
                continue;
            ARRAY<VECTOR<unsigned char,TV::m+1> > triangulation;
            Triangulate(2*(TV::m+1),polytope,triangulation);
            VECTOR<TV,2*TV::m+2> df;
            MATRIX<MATRIX<T,TV::m>,2*TV::m+2> ddf;
            VECTOR<int,2*TV::m+2> nodes=nodes1.Append_Elements(nodes2);
            for(int i=0;i<triangulation.m;i++)
                Integrate_Simplex(triangulation(i),particles.X.Subset(nodes),undeformed_phi.Subset(nodes1),df,ddf);
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
