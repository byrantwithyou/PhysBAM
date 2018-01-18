//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Math_Tools/cyclic_shift.h>
#include <Core/Matrices/IDENTITY_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/ZERO_MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Forces/SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input,
    T stiffness_coefficient,T friction)
    :BASE(particles_input),stiffness_coefficient(stiffness_coefficient),
    friction(friction)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
~SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION()
{
}
namespace{
template<class CP>
auto Attachment(const CP& c)
{
    if(c.diff_entry.m) return c.diff_entry.Last().Y;
    return c.Y0;
}
template<class CP>
int Element(const CP& c)
{
    if(c.diff_entry.m) return c.diff_entry.Last().e;
    return c.e0;
}

template<class TV,class CP,class TV_INT>
void Evalute_Derivatives(TV& dY,TV& dw,const CP& c,const TV& dZ,ARRAY_VIEW<const TV> dX,const ARRAY<TV_INT>& elements)
{
    dY=dX.Subset(elements(c.e0)).Weighted_Sum(c.w0);
    TV dYp;
    for(int i=0;i<c.diff_entry.m;i++){
        dYp=dY;
        const auto& de=c.diff_entry(i);
        dY=de.dYdI(0)*dY+de.dYdI(1)*dZ;
        for(int i=0;i<TV::m;i++)
            dY+=de.dYdI(2+i)*dX(elements(de.e)(i));}

    if(c.diff_entry.m>0){
        dw=c.dwdI(0)*dYp+c.dwdI(1)*dZ;
        for(int i=0;i<TV::m;i++)
            dw+=c.dwdI(2+i)*dX(elements(c.diff_entry.Last().e)(i));}
    else dw=TV();
}

// return:
// -2 = solution is previous attachment
// -1 = solution point in interior. output: Y
// 0,1,2 = exit edge opposite this vertex. output: Y
template<class T,class TV>
int Find_Next_Triangle(TV& Y,VECTOR<MATRIX<T,TV::m>,5>& dYdI,
    TV& w,VECTOR<MATRIX<T,TV::m>,5>& dwdI,
    const TV& X,const TV& Z,const TV& A,const TV& B,const TV& C,T friction,bool& inside)
{
    typedef DIFF_LAYOUT<T,TV::m,TV::m,TV::m,TV::m,TV::m> LAYOUT;
    auto XX=Diff_From_Var<LAYOUT,0>(X);
    auto ZZ=Diff_From_Var<LAYOUT,1>(Z);
    auto AA=Diff_From_Var<LAYOUT,2>(A);
    auto BB=Diff_From_Var<LAYOUT,3>(B);
    auto CC=Diff_From_Var<LAYOUT,4>(C);
    auto u=AA-CC;
    auto v=BB-CC;
    auto x=XX-CC;
    auto z=ZZ-CC;
    auto N=u.Cross(v);
    auto N_mag=N.Magnitude();
    auto n=N/N_mag;
    auto d=-z.Dot(n);
    inside=d.x>=0;
    auto p=z+d*n;
    auto td=abs(d)*friction; // max dist from solution to projection point
    auto xy2=(x-p).Magnitude_Squared();
    if(xy2.x<sqr(td.x)) return -2; // Input is feasible
    auto a=td/sqrt(xy2);
    auto y=(x-p)*a+p; // y = a * u + b * v;
    auto M00=u.Dot(u);
    auto M01=u.Dot(v);
    auto M11=v.Dot(v);
    auto My0=u.Dot(y);
    auto My1=v.Dot(y);
    auto Mx0=u.Dot(x);
    auto Mx1=v.Dot(x);
    auto DM=M00*M11-M01*M01;
    auto q0=(M11*My0-M01*My1)/DM;
    auto q1=(M00*My1-M01*My0)/DM;
    auto r0=(M11*Mx0-M01*Mx1)/DM;
    auto r1=(M00*Mx1-M01*Mx0)/DM;
    auto YY=y+CC;

    int cs=-1;
    if(q0.x<0){
        cs=0;
        q1=r1+(q1-r1)*(r0/(r0-q0));
        q0=q0*0;
        YY=q1*v+CC;}
    if(q1.x<0){
        cs=1;
        q0=r0+(q0-r0)*(r1/(r1-q1));
        q1=q1*0;
        YY=q0*u+CC;}
    if(q0.x+q1.x>1){
        cs=2;
        auto sw=r0+r1;
        auto sq=q0+q1;
        auto a=(1-sw)/(sq-sw);
        q0=r0+(q0-r0)*a;
        q1=r1+(q1-r1)*a;
        YY=q0*u+q1*v+CC;}
    Y=YY.x;
    Extract<0>(dYdI,YY.dx);
    auto WW=Make_Vector<T>(q0,q1,(T)1-q0-q1);
    w=WW.x;
    Extract<0>(dwdI,WW.dx);
    return cs;
}

// return:
// -2 = solution is previous attachment
// -1 = solution point in interior. output: Y
// 0,1 = exit at this endpoint.
template<class T,class TV>
int Stuck_On_Edge(TV& Y,VECTOR<MATRIX<T,TV::m>,4>& dYdI,T& w,VECTOR<TV,4>& dwdI,
    const TV& X,const TV& Z,const TV& A,const TV& B,T friction)
{
    typedef DIFF_LAYOUT<T,TV::m,TV::m,TV::m,TV::m> LAYOUT;
    auto XX=Diff_From_Var<LAYOUT,0>(X);
    auto ZZ=Diff_From_Var<LAYOUT,1>(Z);
    auto AA=Diff_From_Var<LAYOUT,2>(A);
    auto BB=Diff_From_Var<LAYOUT,3>(B);
    auto u=AA-BB;
    auto x=XX-BB;
    auto z=ZZ-BB;
    auto u2=u.Dot(u);
    auto p=u.Dot(z)/u2*u;
    auto d2=(z-p).Magnitude_Squared();
    auto td2=d2*sqr(friction); // max dist from solution to projection point
    auto xy2=(x-p).Magnitude_Squared();
    if(xy2.x<td2.x) return -2; // Input is feasible
    auto a=sqrt(td2/xy2);
    auto y=(x-p)*a+p;
    auto q=u.Dot(y)/u2;
    if(q.x<0){
        Y=BB.x;
        Extract<0>(dYdI,BB.dx);
        w=0;
        dwdI.Fill(TV());
        return 1;}
    if(q.x>1){
        Y=AA.x;
        Extract<0>(dYdI,AA.dx);
        w=1;
        dwdI.Fill(TV());
        return 0;}
    auto YY=q*u+BB;
    Y=YY.x;
    Extract<0>(dYdI,YY.dx);
    w=q.x;
    Extract<0>(dwdI,q.dx);
    return -1;
}

// return:
// (e,v) = exit in trangle e; on edge (p,v)
// (v,-1) = exit along edge (p,v)
// (-1,-1) = no exit from vertex
template<class T,class TV>
VECTOR<int,2> Handle_Vertex(const TV& Z,const TRIANGULATED_SURFACE<T>& ts,int p,T friction)
{
    // Try leaving vertex via one of the neighboring triangles
    const ARRAY<int>& ie=(*ts.mesh.incident_elements)(p);
    TV O=ts.particles.X(p),z=Z-O;
    for(int i=0;i<ie.m;i++){
        VECTOR<int,3> e(ts.mesh.elements(ie(i)));
        if(e.x!=p){
            cyclic_shift(e);
            if(e.x!=p) cyclic_shift(e);}
        PHYSBAM_ASSERT(e.x==p);
        TV A=ts.particles.X(e.y),u=A-O;
        TV B=ts.particles.X(e.z),v=B-O;
        TV n=u.Cross(v),un=u.Cross(n),vn=v.Cross(n);
        // z = a*u+b*v+c*n;
        T n2=n.Magnitude_Squared();
        T a=z.Dot(vn)/n2;
        T b=-z.Dot(un)/n2;
        if(a>=0 && b>=0) // can exit via this triangle.
            return {ie(i),e.y};}

    // Try leaving vertex via one of the neighboring edges
    const ARRAY<int>& nn=(*ts.mesh.neighbor_nodes)(p);
    for(int i=0;i<nn.m;i++){
        TV A=ts.particles.X(nn(i)),u=A-O;
        if(u.Dot(z)>0)
            return {nn(i),-1};}

    return {-1,-1};
}

// return active
template<class T,class TV,class CP>
void Relax_Attachment(CP& c,const TV& Z,const TRIANGULATED_SURFACE<T>& ts,T friction)
{
    typedef VECTOR<int,TV::m> TV_INT;
    enum STATE {tri_int,tri_edge,edge_int,point} state=tri_int;
    const ARRAY<ARRAY<int> >& adjacent_elements=*ts.mesh.adjacent_elements;

    c.diff_entry.Remove_All();
    int e0=c.e0;
    TV Y=ts.particles.X.Subset(ts.mesh.elements(e0)).Weighted_Sum(c.w0);
    // For tri_edge, vertices from previous element
    // For edge_int, vertices of the edge
    VECTOR<int,2> edge;
    int p=-1;
    typename SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::DIFF_ENTRY de;

    c.Y0=Y;
    c.w=c.w0;
    c.dwdI.Fill(MATRIX<T,3>());

    while(1)
    {
//        LOG::printf("STATE %i\n",state);
        switch(state)
        {
            case tri_int:
            case tri_edge:{
                TV_INT e=ts.mesh.elements(e0);
                VECTOR<TV,TV::m> P(ts.particles.X.Subset(e));
                de.e=e0;
                int ret=Find_Next_Triangle(de.Y,de.dYdI,c.w,c.dwdI,
                    Y,Z,P(0),P(1),P(2),friction,c.active);
                if(ret==-2){/*puts("RET U");*/return;}
                Y=de.Y;
                c.diff_entry.Append(de);
                if(ret==-1){/*puts("RET T");*/return;}
                int rem=e(ret),next_e=-1;
                VECTOR<int,2> new_edge=e.Remove_Index(ret).Sorted();
                const ARRAY<int>& a=adjacent_elements(e0);
                for(int i=0;i<a.m;i++)
                    if(!ts.mesh.elements(a(i)).Contains(rem)){
                        e0=next_e=a(i);
                        break;}
                if(next_e==-1){c.active=false;/*puts("RET B");*/return;}
                if(state==tri_edge && new_edge==edge) state=edge_int;
                else state=tri_edge;
                edge=new_edge;
                break;}

            case edge_int:{
                VECTOR<TV,2> P(ts.particles.X.Subset(edge));
                VECTOR<MATRIX<T,TV::m>,4> dYdE;
                T w;
                VECTOR<TV,4> dwdE;
                int ret=Stuck_On_Edge(de.Y,dYdE,w,dwdE,Y,Z,P(0),P(1),friction);

                // Get both triangles neighboring edge
                TV_INT e[2];
                int tris[2]={-1,-1},next_tri=0;
                de.e=-1;
                const ARRAY<int>& ie=(*ts.mesh.incident_elements)(edge.x);
                for(int i=0;i<ie.m;i++){
                    TV_INT te=ts.mesh.elements(ie(i));
                    if(te.Contains(edge.y)){
                        e[next_tri]=te;
                        tris[next_tri++]=ie(i);
                        if(next_tri==2) break;}}
                de.e=tris[0];

                // Determine inside/outside based on two neighbor triangles.
                VECTOR<TV,3> V0(ts.particles.X.Subset(e[0]));
                VECTOR<TV,3> V1(ts.particles.X.Subset(e[1]));
                TV n0=(V0.y-V0.x).Cross(V0.z-V0.x);
                TV n1=(V1.y-V1.x).Cross(V1.z-V1.x);
                bool i0=n0.Dot(Z-V0.x)<=0;
                bool i1=n1.Dot(Z-V1.x)<=0;
                if(i0==i1) c.active=i0; // Both triangles agree on inside/outside
                else{
                    TV OP=ts.particles.X(e[1].Sum()-edge.Sum());
                    c.active=n0.Dot(OP-V0.x)>0;} // outside if convex

                if(ret==-2){/*puts("RET F");*/return;}

                PHYSBAM_ASSERT(de.e>=0);

                int idx[2]={e[0].Find(edge.x),e[0].Find(edge.y)},miss=3-idx[0]-idx[1];

                Y=de.Y;
                de.dYdI(0)=dYdE(0);
                de.dYdI(1)=dYdE(1);
                de.dYdI(2+idx[0])=dYdE(2);
                de.dYdI(2+idx[1])=dYdE(3);
                de.dYdI(2+miss)=MATRIX<T,3>();
                c.w(idx[0])=w;
                c.w(idx[1])=1-w;
                c.w(miss)=0;
                for(int i=0;i<2;i++)
                    for(int j=0,s=1;j<2;j++,s-=2){
                        c.dwdI(i).Set_Row(idx[j],dwdE(i)*s);
                        c.dwdI(2+idx[i]).Set_Row(idx[j],dwdE(2+i)*s);}
                c.dwdI(2+miss)=MATRIX<T,3>();
                c.diff_entry.Append(de);

                if(ret==-1){/*puts("RET E");*/return;}

                // NOTE: This breaks the dependency chain (!!)

                p=edge(ret);
                state=point;
                break;}

            case point:
                auto ret=Handle_Vertex(Z,ts,p,friction);
                if(ret.x==-1){
                    c.active=ts.Signed_Solid_Angle_Of_Triangle_Web(Y,p)>0;
                    //puts("RET V");
                    int e=c.diff_entry.Last().e;
                    int j=ts.mesh.elements(e).Find(p);
                    c.dwdI.Fill(MATRIX<T,3>());
                    c.w=TV();
                    c.w(j)=1;
                    return;}
                if(ret.y==-1){
                    state=edge_int;
                    edge={p,ret.x};}
                else{
                    state=tri_edge;
                    e0=ret.x;
                    edge={p,ret.y};}
        }
    }
}
template<class T,class TV,class CP>
void Relax_Attachment(CP& c,const TV& Z,const SEGMENTED_CURVE_2D<T>& sc,T friction)
{
    PHYSBAM_NOT_IMPLEMENTED();
}

template<class T,class TV,class CP>
void Visualize(const CP& c,const TV& Z,const TRIANGULATED_SURFACE<T>& ts)
{
    Dump_Surface(ts,TV(.5,.5,.5));
    TV Y=c.Y0;
    Add_Debug_Particle(Y,TV(1,0,0));
    Add_Debug_Particle(Z,TV(1,0,0));
    for(int i=0;i<c.diff_entry.m;i++){
        const auto& de=c.diff_entry(i);
        Add_Debug_Object(VECTOR<TV,2>(Y,de.Y),TV(1,0,1));
        Add_Debug_Particle(de.Y,TV(1,1,0));
        Y=de.Y;}
}

template<class TV,class CP,class SURF>
void Test_Weights(const CP& c,const SURF& ts)
{
    TV Y0=c.Y0,Y1=c.Y0;
    if(c.diff_entry.m>0){
        const auto& le=c.diff_entry.Last();
        Y1=ts.particles.X.Subset(ts.mesh.elements(le.e)).Weighted_Sum(c.w);
        Y0=le.Y;}
    LOG::printf("Y vs w: %P %P   %P\n",Y0,Y1,(Y0-Y1).Magnitude());
}

template<class TV,class TV_INT,class CP,class T>
void Test_Diff(const CP& c0,const CP& c1,
    const TV& dZ,const ARRAY<TV>& dX,const ARRAY<TV_INT>& elements,T eps)
{
    if(c0.diff_entry.m!=c1.diff_entry.m){
        LOG::printf("SKIP Test_Diff; paths differ: len %i != %i\n",
            c0.diff_entry.m,c1.diff_entry.m);
        return;}
    for(int i=0;i<c0.diff_entry.m;i++)
        if(c0.diff_entry(i).e!=c1.diff_entry(i).e){
            LOG::printf("SKIP Test_Diff; paths differ: tri(%i) %i != %i\n",
                i,c0.diff_entry(i).e,c1.diff_entry(i).e);
            return;}

    TV dY0,dY1,dw0,dw1;
    Evalute_Derivatives<TV>(dY0,dw0,c0,dZ,dX,elements);
    Evalute_Derivatives<TV>(dY1,dw1,c1,dZ,dX,elements);
    TV da=(Attachment(c1)-Attachment(c0))/eps;
    TV db=(dY0+dY1)/(2*eps);
    LOG::printf("diff test: %.16P %.16P %P\n",da.Magnitude(),db.Magnitude(),(da-db).Magnitude());

    TV dc=(c1.w-c0.w)/eps;
    TV dd=(dw0+dw1)/(2*eps);
    LOG::printf("diff weight: %.16P %.16P %P\n",dc.Magnitude(),dd.Magnitude(),(dc-dd).Magnitude());
}
}
//#####################################################################
// Function Test_Relax
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Test_Relax(int cp)
{
    T eps=1e-6;
    
    COLLISION_PAIR c1=collision_pairs(cp),c2=collision_pairs(cp);
    auto& ts=*surfaces(c1.s);

    RANDOM_NUMBERS<T> rand;
    ARRAY<TV> dX(particles.X.m),X0(particles.X);
    rand.Fill_Uniform(dX,-eps,eps);

    TV Z=particles.X(c1.p);
    ::Relax_Attachment(c1,Z,ts,friction);

    Test_Weights<TV>(c1,ts);
    LOG::printf("ret %i\n",c1.active);
    particles.X+=dX;
    TV Z2=particles.X(c2.p),dZ=Z2-Z;
    ::Relax_Attachment(c2,Z2,ts,friction);
    particles.X=X0;
    LOG::printf("ret %i\n",c2.active);
    Test_Diff(c1,c2,dZ,dX,ts.mesh.elements,eps);
    TV Y=Attachment(c1);
    // if(c1.active!=(torus.Normal(Y).Dot(Z-Y)<=0)){
    //     LOG::printf("INSIDE TEST DOES NOT MATCH %i %g\n",c1.active,torus.Signed_Distance(Z));
    //     Visualize(c1,Z,ts);
    // }
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active){
            const ARRAY<TV_INT>& elements=surfaces(c.s)->mesh.elements;
            TV_INT e=elements(Element(c));
            TV j=stiffness_coefficient*(particles.X(c.p)-Attachment(c));
            F(c.p)-=j;
            for(int k=0;k<TV::m;k++)
                F(e(k))+=c.w(k)*j;}}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    get_candidates();

    for(int i=0;i<collision_pairs.m;i++)
        Relax_Attachment(i);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active){
            const ARRAY<TV_INT>& elements=surfaces(c.s)->mesh.elements;
            TV_INT e=elements(Element(c));
            TV j=stiffness_coefficient*(particles.X(c.p)-Attachment(c));
            TV dY,dw,dZ=V(c.p);
            Evalute_Derivatives(dY,dw,c,dZ,V,elements);
            TV dj=stiffness_coefficient*(V(c.p)-dY);
            F(c.p)-=dj;
            for(int k=0;k<TV::m;k++)
                F(e(k))+=c.w(k)*dj+dw(k)*j;}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Potential_Energy(const T time) const
{
    T pe=0;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active)
            pe+=(T).5*stiffness_coefficient*(particles.X(c.p)-Attachment(c)).Magnitude_Squared();}
    return pe;
}
//#####################################################################
// Function Relax_Attachment
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Relax_Attachment(int cp)
{
    COLLISION_PAIR& c=collision_pairs(cp);
    TV Z=particles.X(c.p);
    ::Relax_Attachment(c,Z,*surfaces(c.s),friction);
}
//#####################################################################
// Function Update_Attachments_And_Prune_Pairs
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Attachments_And_Prune_Pairs()
{
    int k=0;
    for(int i=0;i<collision_pairs.m;i++){
        COLLISION_PAIR c=collision_pairs(i);
        if(c.active){
            if(c.diff_entry.m>0){
                c.w0=c.w;
                c.e0=c.diff_entry.Last().e;
                c.Y0=c.diff_entry.Last().Y;}
            collision_pairs(k++)=c;}
        else hash.Delete({c.p,c.s});}
    collision_pairs.Resize(k);
}
//#####################################################################
// Function Add_Pair
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Pair(int p,int s,const TV& w0,int e0)
{
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({p,s})) return;
    const auto& ts=*surfaces(s);
    TV Y=ts.particles.X.Subset(ts.mesh.elements(e0)).Weighted_Sum(w0);
    COLLISION_PAIR c={p,s,w0,e0,Y};
    collision_pairs.Append(c);
    hash.Insert({p,s});
}
//#####################################################################
// Function Add_Surface
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Surface(T_SURFACE& surface)
{
    surface.mesh.Initialize_Incident_Elements();
    surface.mesh.Initialize_Neighbor_Nodes();
    surface.mesh.Initialize_Adjacent_Elements();
    int si=surfaces.Append(&surface);
    for(int i=0;i<surface.mesh.elements.m;i++)
        object_from_element.Set(surface.mesh.elements(i).Sorted(),{si,i});
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
CFL_Strain_Rate() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
namespace PhysBAM{
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,2> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,3> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,2> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,3> >;
}
