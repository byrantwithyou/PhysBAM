//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/cyclic_shift.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/RELAX_ATTACHMENT_MESH.h>
using namespace PhysBAM;
// return:
// -2 = solution is previous attachment
// -1 = solution point in interior. output: Y
// 0,1,2 = exit edge opposite this vertex. output: Y
template<class T,class TV>
int Find_Next_Triangle(TV& Y,VECTOR<MATRIX<T,TV::m>,5>& dYdI,
    TV& w,VECTOR<MATRIX<T,TV::m>,5>& dwdI,
    const TV& X,const TV& Z,const TV& A,const TV& B,const TV& C,T friction,bool& inside)
{
    static const T small_number=std::numeric_limits<T>::epsilon()*128;
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
    if(xy2.x<sqr(td.x)+small_number) return -2; // Input is feasible
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
        if(sw.x<=1){
            auto a=(1-sw)/(sq-sw);
            q0=r0+(q0-r0)*a;
            q1=r1+(q1-r1)*a;}
        else{
            q0=q0*0+r0;
            q1=q1*0+r1;}
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
    static const T small_number=std::numeric_limits<T>::epsilon()*128;
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
    if(xy2.x<td2.x+small_number) return -2; // Input is feasible
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
    static const T small_number=std::numeric_limits<T>::epsilon()*128;
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
        if(a>=small_number && b>=small_number) // can exit via this triangle.
            return {ie(i),e.y};}

    // Try leaving vertex via one of the neighboring edges
    const ARRAY<int>& nn=(*ts.mesh.neighbor_nodes)(p);
    for(int i=0;i<nn.m;i++){
        TV A=ts.particles.X(nn(i)),u=A-O;
        if(u.Dot(z)>small_number)
            return {nn(i),-1};}

    return {-1,-1};
}

// return active
template<class T,class TV>
bool Relax_Attachment(RELAX_ATTACHMENT_MESH<TV>& c,int e0,const TV& w0,const TV& Z,
    const TRIANGULATED_SURFACE<T>& ts,int ex_pt,T friction)
{
    typedef VECTOR<int,TV::m> TV_INT;
    enum STATE {tri_int,tri_edge,edge_int,point} state=tri_int;
    const ARRAY<ARRAY<int> >& adjacent_elements=*ts.mesh.adjacent_elements;
    static const T small_number=std::numeric_limits<T>::epsilon()*128;

    c.diff_entry.Remove_All();
    TV Y=ts.particles.X.Subset(ts.mesh.elements(e0)).Weighted_Sum(w0);
    // For tri_edge, vertices from previous element
    // For edge_int, vertices of the edge
    VECTOR<int,2> edge;
    int p=-1;
    typename RELAX_ATTACHMENT_MESH<TV>::DIFF_ENTRY de;

    c.w=w0;
    c.dwdI.Fill(MATRIX<T,3>());

    bool done=false;
    while(!done)
    {
//        LOG::printf("STATE %i\n",state);
        switch(state)
        {
            case tri_int:
            case tri_edge:{
                TV_INT e=ts.mesh.elements(e0);
                if(e.Contains(ex_pt)) return false;
                VECTOR<TV,TV::m> P(ts.particles.X.Subset(e));
                de.e=e0;
                int ret=Find_Next_Triangle(de.Y,de.dYdI,c.w,c.dwdI,
                    Y,Z,P(0),P(1),P(2),friction,c.active);
                if(ret==-2){/*puts("RET U");*/done=true;break;}
                Y=de.Y;
                c.diff_entry.Append(de);
                if(ret==-1){/*puts("RET T");*/done=true;break;}
                bool use_vertex=false;
                for(int i=0;i<TV::m;i++)
                    if(c.w(i)>1-small_number){
                        p=e(i);
                        state=point;
                        use_vertex=true;
                        // BUG: what about de.dYdI and c.dwdI?
                        LOG::printf("CLAMP TO VERTEX\n");
                        break;}
                if(use_vertex) break;
                int rem=e(ret),next_e=-1;
                VECTOR<int,2> new_edge=e.Remove_Index(ret).Sorted();
                const ARRAY<int>& a=adjacent_elements(e0);
                for(int i=0;i<a.m;i++)
                    if(!ts.mesh.elements(a(i)).Contains(rem)){
                        e0=next_e=a(i);
                        break;}
                if(next_e==-1) return false;
                if(state==tri_edge && new_edge==edge) state=edge_int;
                else state=tri_edge;
                edge=new_edge;
                break;}

            case edge_int:{
                if(edge.Contains(ex_pt)) return false;
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

                if(ret==-2){/*puts("RET F");*/done=true;break;}
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

                if(ret==-1){/*puts("RET E");*/done=true;break;}
                // NOTE: This breaks the dependency chain (!!)

                p=edge(ret);
                state=point;
                break;}

            case point:
                if(p==ex_pt) return false;
                auto ret=Handle_Vertex(Z,ts,p,friction);
                if(ret.x==-1){
                    c.active=ts.Signed_Solid_Angle_Of_Triangle_Web(Y,p)>0;
                    //puts("RET V");
                    int e=c.diff_entry.Last().e;
                    int j=ts.mesh.elements(e).Find(p);
                    c.dwdI.Fill(MATRIX<T,3>());
                    c.w=TV();
                    c.w(j)=1;
                    done=true;
                    break;}
                if(ret.y==-1){
                    state=edge_int;
                    edge={p,ret.x};}
                else{
                    state=tri_edge;
                    e0=ret.x;
                    edge={p,ret.y};}
        }
    }
    c.Y=Y;
    c.e=de.e;
    return c.active;
}
template<class T,class TV>
bool Relax_Attachment(RELAX_ATTACHMENT_MESH<TV>& c,int e0,const TV& w0,const TV& Z,
    const SEGMENTED_CURVE_2D<T>& sc,int ex_pt,T friction)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Relax
//#####################################################################
template<class TV> void RELAX_ATTACHMENT_MESH<TV>::
Relax(int e0,const TV& w0,const TV& Z,const T_SURFACE& surface,int ex_pt,T friction)
{
    active=Relax_Attachment(*this,e0,w0,Z,surface,ex_pt,friction);
}
template class RELAX_ATTACHMENT_MESH<VECTOR<float,2> >;
template class RELAX_ATTACHMENT_MESH<VECTOR<float,3> >;
template class RELAX_ATTACHMENT_MESH<VECTOR<double,2> >;
template class RELAX_ATTACHMENT_MESH<VECTOR<double,3> >;
