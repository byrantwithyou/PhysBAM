//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "FLUID_LAYOUT_FEM.h"
#include <list>
#include <map>

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> FLUID_LAYOUT_FEM<VECTOR<T,2> >::
FLUID_LAYOUT_FEM(): area_hidden(*new TRIANGULATED_AREA<T>)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FLUID_LAYOUT_FEM<VECTOR<T,2> >::
~FLUID_LAYOUT_FEM()
{
    delete &area_hidden;
}
//#####################################################################
// Function Generate_Pipe
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Generate_Pipe(PIPE_ID pipe,const PARSE_DATA_FEM<TV,TW>& pd,const CONNECTION& con)
{
    typedef VECTOR<PARTICLE_ID,3> E;
    PARTICLE_ID base=Number_Particles();
    const ARRAY<PARTICLE_ID>* f0=con.Get_Pointer({pd.pipes(pipe).x,pipe}),*f1=con.Get_Pointer({pd.pipes(pipe).y,pipe});
    if(!f0 || !f1) return;
    TV v0=X((*f0)(pd.half_width)),v1=X((*f1)(pd.half_width));
    TV d=v1-v0;
    T l=d.Normalize();
    if(!l) return;
    TV t(-d.y,d.x);
    int width=2*pd.half_width;
    int height=rint(l/pd.unit_length);
    if(height==0) height=1;
    if(height>1)
        Add_Particles((width+1)*(height-1));
    auto pid=[base,width,height,&pd,f0,f1](int i,int j)
    {
        if(i==0) return (*f0)(pd.half_width+j);
        else if(i==height) return (*f1)(pd.half_width+j);
        else return base+(i-1)*(width+1)+j+pd.half_width;
    };

    TV inc_i=pd.unit_length*d,inc_j=pd.unit_length*t;
    TV next=v0+inc_i-pd.half_width*inc_j;
    PARTICLE_ID k=base;
    for(int i=1;i<height;i++){
        TV pt=next;
        for(int j=-pd.half_width;j<=pd.half_width;j++){
            X(k++)=pt;
            pt+=inc_j;}
        next+=inc_i;}

    TRIANGLE_ID canonical_element=Number_Triangles();
    for(int i=0;i<height;i++){
        for(int j=-pd.half_width+1;j<=pd.half_width;j++){
            PARTICLE_ID a=pid(i,j),b=pid(i+1,j),c=pid(i,j-1),d=pid(i+1,j-1);
            Append_Triangle(E(a,c,b));
            Append_Triangle(E(d,b,c));
            elem_data.Append({blocks.m,j>0?3:0});
            elem_data.Append({blocks.m,j>0?0:3});}
        for(int j=-pd.half_width;j<=pd.half_width;j++){
            node_blocks(pid(i+(j<0),j))=blocks.m;
            PARTICLE_ID cp=pid(i+(j>=0),j);
            if((i==0 || i==height-1) && node_blocks(cp)<BLOCK_ID())
                node_blocks(cp)=blocks.m;}
        blocks.Append({pipe});}

    bool irregular_last=abs(l-height*pd.unit_length)>(T)1e-6*l;
    if(irregular_last) blocks.Last().pipe_id=PIPE_ID(-1);

    pipes(pipe)={-inc_j,inc_i,canonical_element};
}
//#####################################################################
// Function Mark_BC
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Mark_BC(const ARRAY<PARTICLE_ID>& pindices,BC_ID bc_id)
{
    particle_bc_map.Set(pindices(0),bc_id);
    for(int i=1;i<pindices.m;i++){
        PARTICLE_ID p0=pindices(i),p1=pindices(i-1);
        if(p0>p1) std::swap(p0,p1);
        particle_bc_map.Set(pindices(i),bc_id);
        bc_map.Set({p0,p1},bc_id);}
}
//#####################################################################
// Function Generate_End
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Generate_End(VERTEX_ID i,PIPE_ID pipe,const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con)
{
    VERTEX_ID j=pd.pipes(pipe).x;
    if(j==i) j=pd.pipes(pipe).y;
    TV p=pd.pts(i).pt;
    TV v=pd.pts(j).pt-pd.pts(i).pt;
    TV u=TV(-v.y,v.x).Normalized();
    PARTICLE_ID base=Add_Particles(2*pd.half_width+1);
    ARRAY<PARTICLE_ID>& indices=con.Get_Or_Insert({i,pipe});
    for(int k=-pd.half_width;k<=pd.half_width;k++){
        PARTICLE_ID pid=base+k+pd.half_width;
        indices.Append(pid);
        X(pid)=p+k*pd.unit_length*u;}
    if(pd.pipes(pipe).x!=i) indices.Reverse();
    Mark_BC(indices,pd.pts(i).bc_id);
}
//#####################################################################
// Function Wedge
//#####################################################################
template<class T> VECTOR<VECTOR<T,2>,3> FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Wedge(const TV& joint,const TV& p0,const TV& p1,int half_width,T unit_length) const
{
    TV e0=(p0-joint).Normalized(),e1=(p1-joint).Normalized();
    TV m=(e0+e1).Normalized();
    T a=e0.Cross(m)(0);
    if(abs(a)<1e-6){
        TV v(-e0(1),e0(0));
        return VECTOR<TV,3>(joint+half_width*unit_length*v,e0,e1);}
    else{
        T s=half_width*unit_length/a;
        return VECTOR<TV,3>(joint+s*m,(s*m.Dot(e0)*e0-s*m).Normalized(),(s*m.Dot(e1)*e1-s*m).Normalized());}
}
//#####################################################################
// Function Sample_Interpolated
//#####################################################################
template<class T> ARRAY<PARTICLE_ID> FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Sample_Interpolated(T s,const ARRAY<PARTICLE_ID>& side0,const ARRAY<PARTICLE_ID>& side1,T unit_length)
{
    ARRAY<T> l0(side0.m),l1(side1.m);
    l0(0)=0;
    for(int j=1;j<side0.m;j++)
        l0(j)=l0(j-1)+(X(side0(j))-X(side0(j-1))).Magnitude();
    l1(0)=0;
    for(int j=1;j<side1.m;j++)
        l1(j)=l1(j-1)+(X(side1(j))-X(side1(j-1))).Magnitude();
    auto point=[this](int j,T t,const ARRAY<PARTICLE_ID>& side)
    {
        if(j>=side.m-1) return X(side(j));
        TV p=X(side(j));
        return p+t*(X(side(j+1))-p);
    };
    auto loc=[point](T lambda,const ARRAY<PARTICLE_ID>& side,const ARRAY<T>& l)
    {
        T dist=lambda*l.Last();
        auto iter=std::lower_bound(l.begin(),l.end(),dist);
        PHYSBAM_ASSERT(iter!=l.end());
        int j=iter-l.begin();
        if(*iter==dist) return point(j,0,side);
        else{
            T cur=*iter,prev=*(iter-1);
            return point(j-1,(dist-prev)/(cur-prev),side);}
    };

    std::list<PAIR<T,TV> > verts;
    TV p0=(1-s)*loc(0,side0,l0)+s*loc(0,side1,l1);
    TV p1=(1-s)*loc(1,side0,l0)+s*loc(1,side1,l1);
    auto begin=verts.insert(verts.end(),{0,p0});
    auto end=verts.insert(verts.end(),{1,p1});
    T max_len=(p1-p0).Magnitude();
    T min_percent=0.2;
    while(max_len>(1+min_percent)*unit_length){
        T lambda=(begin->x+end->x)*0.5;
        TV v=(1-s)*loc(lambda,side0,l0)+s*loc(lambda,side1,l1);
        verts.insert(end,{lambda,v});
        max_len=0;
        for(auto k=verts.begin();k!=verts.end();k++){
            if(k==verts.begin()) continue;
            auto prev=k;
            prev--;
            T d=(k->y-prev->y).Magnitude();
            if(d>max_len){
                begin=prev;
                end=k;
                max_len=d;}}}
    PARTICLE_ID base=Add_Particles(verts.size());
    ARRAY<PARTICLE_ID> indices(verts.size());
    int j=0;
    for(auto iter=verts.begin();iter!=verts.end();iter++){
        X(base+j)=iter->y;
        indices(j)=base+j;
        j++;}
    return indices;
}
//#####################################################################
// Function Merge_Interpolated
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Merge_Interpolated(const ARRAY<PARTICLE_ID>& left,const ARRAY<PARTICLE_ID>& right)
{
    typedef VECTOR<PARTICLE_ID,3> E;
    int i=0,j=0,alt=0;
    auto angle=[this](PARTICLE_ID v0,PARTICLE_ID v1,PARTICLE_ID v2)
    {
        TV u=(X(v2)-X(v1)).Normalized();
        TV v=(X(v0)-X(v1)).Normalized();
        return TV::Angle_Between(u,v);
    };
    ARRAY<E> elems;
    auto assign_node=[this](PARTICLE_ID n,BLOCK_ID bid,bool greedy)
    {
        auto& nb=node_blocks(n);
        if(nb<BLOCK_ID() || greedy)
            nb=bid;
    };
    auto max_angle_index=[this,angle](PARTICLE_ID* p)
    {
        T max_angle=-1;
        int max_index=-1;
        for(int j=0;j<3;j++){
            T a=angle(p[(j+2)%3],p[j],p[(j+1)%3]);
            if(a>max_angle){
                max_angle=a;
                max_index=j;}}
        return max_index;
    };
    assign_node(left(i),blocks.m,false);
    assign_node(right(j),blocks.m,true);
    while(i<left.m-1 || j<right.m-1){
        T a0=0;
        if(i+1<left.m) a0=angle(right(j),left(i+1),left(i));
        T a1=0;
        if(j+1<right.m) a1=angle(right(j),right(j+1),left(i));
        int greed;
        if(j+1>=right.m || (i+1<left.m && abs(a0-a1)<1e-6 && alt==0) || (i+1<left.m && a0>a1)){
            PARTICLE_ID p[]={left(i+1),left(i),right(j)};
            int k=max_angle_index(p);
            Append_Triangle(E(p[k],p[(k+1)%3],p[(k+2)%3]));
            greed=i<(left.m-1)/2?1:2;
            assign_node(left(i+1),blocks.m,i>=left.m/2);
            i++;
            alt=1;}
        else{
            PARTICLE_ID p[]={right(j+1),left(i),right(j)};
            int k=max_angle_index(p);
            Append_Triangle(E(p[k],p[(k+1)%3],p[(k+2)%3]));
            greed=j<(right.m-1)/2?2:1;
            assign_node(right(j+1),blocks.m,j<right.m/2);
            j++;
            alt=0;}
        elem_data.Append({blocks.m,greed});}
}
//#####################################################################
// Function Polyline
//#####################################################################
template<class T> ARRAY<PARTICLE_ID> FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Polyline(const ARRAY<TV>& points,T unit_length)
{
    ARRAY<PARTICLE_ID> verts;
    PARTICLE_ID base;
    for(int j=1;j<points.m;j++){
        TV v=points(j)-points(j-1);
        T l=v.Normalize();
        for(int k=0;(k+1)*unit_length<=l+0.5*unit_length;k++){
            base=Add_Particle();
            X(base)=points(j-1)+v*k*unit_length;
            verts.Append(base);}}
    base=Add_Particle();
    X(base)=points.Last();
    verts.Append(base);
    return verts;
}
//#####################################################################
// Function Arc
//#####################################################################
template<class T> PAIR<ARRAY<PARTICLE_ID>,ARRAY<PARTICLE_ID> > FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Arc(const TV& c,const TV& p0,const TV& p1,int half_width,T unit_length,
    const TV& dir0,T n0,const TV& dir1,T n1)
{
    ARRAY<PARTICLE_ID> side0,side1;
    PARTICLE_ID base;
    if(n0){
        base=Add_Particles(2);
        X(base)=c+dir0*unit_length*n0;
        X(base+1)=p0+dir0*unit_length*n0;
        side0.Append(base);
        side1.Append(base+1);}
    base=Add_Particle();
    X(base)=c;
    side0.Append(base);
    T unit_angle=2*asin((T)1/(half_width*4));
    TV v0=(p0-c).Normalized(),v1=(p1-c).Normalized();
    T total=acos(v0.Dot(v1));
    for(int j=0;(j+1)*unit_angle<=total+0.5*unit_angle;j++){
        T a=-j*unit_angle;
        TV v(cos(a)*v0(0)-sin(a)*v0(1),sin(a)*v0(0)+cos(a)*v0(1));
        base=Add_Particle();
        X(base)=c+v*2*half_width*unit_length;
        side1.Append(base);}
    base=Add_Particle();
    X(base)=p1;
    side1.Append(base);
    if(n1){
        base=Add_Particles(2);
        X(base)=c+dir1*unit_length*n1;
        X(base+1)=p1+dir1*unit_length*n1;
        side0.Append(base);
        side1.Append(base+1);}
    return {side0,side1};
}
//#####################################################################
// Function Corner
//#####################################################################
template<class T> PAIR<ARRAY<PARTICLE_ID>,ARRAY<PARTICLE_ID> > FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Corner(const TV& c,const TV& elbow,const TV& p0,const TV& p1,T unit_length,
    const TV& dir0,T n0,const TV& dir1,T n1)
{
    ARRAY<PARTICLE_ID> side0,side1;
    PARTICLE_ID base;
    if(n0){
        base=Add_Particles(2);
        X(base)=c+dir0*unit_length*n0;
        X(base+1)=p0+dir0*unit_length*n0;
        side0.Append(base);
        side1.Append(base+1);}
    base=Add_Particle();
    X(base)=c;
    side0.Append(base);
    TV start=p0;
    TV v=elbow-start;
    T total=v.Normalize();
    for(int j=0;(j+1)*unit_length<=total+0.5*unit_length;j++){
        base=Add_Particle();
        X(base)=start+v*j*unit_length;
        side1.Append(base);}
    TV end=p1;
    v=end-elbow;
    total=v.Normalize();
    for(int j=0;(j+1)*unit_length<=total+0.5*unit_length;j++){
        base=Add_Particle();
        X(base)=elbow+v*j*unit_length;
        side1.Append(base);}
    base=Add_Particle();
    X(base)=end;
    side1.Append(base);
    if(n1){
        base=Add_Particles(2);
        X(base)=c+dir1*unit_length*n1;
        X(base+1)=p1+dir1*unit_length*n1;
        side0.Append(base);
        side1.Append(base+1);}
    return {side0,side1};
}
//#####################################################################
// Function Weld
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Weld(int n,const ARRAY<PARTICLE_ID>& side0,const ARRAY<PARTICLE_ID>& side1,T unit_length,ARRAY<PARTICLE_ID>& f0,ARRAY<PARTICLE_ID>& f1)
{
    f0.Append(side0(0));
    f1.Append(side0.Last());
    ARRAY<PARTICLE_ID> prev=side0;
    for(int j=1;j<n;j++){
        T s=(T)j/n;
        ARRAY<PARTICLE_ID> cur=Sample_Interpolated(s,side0,side1,unit_length);
        Merge_Interpolated(cur,prev);
        blocks.Append({PIPE_ID(-1)});
        f0.Append(cur(0));
        f1.Append(cur.Last());
        prev=cur;}
    Merge_Interpolated(side1,prev);
    blocks.Append({PIPE_ID(-1)});
    f0.Append(side1(0));
    f1.Append(side1.Last());
}
//#####################################################################
// Function Generate_2_Joint
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Generate_2_Joint(VERTEX_ID i,const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con,JOINT_TYPE jt)
{
    typedef VECTOR<PARTICLE_ID,3> E;
    const ARRAY<PIPE_ID>& pipes=pd.pts(i).joints;
    PHYSBAM_ASSERT(pipes.m==2);
    PIPE_ID p0=pipes(0),p1=pipes(1);
    VERTEX_ID j0=pd.pipes(p0)(0),j1=pd.pipes(p1)(0);
    if(j0==i) j0=pd.pipes(p0)(1);
    if(j1==i) j1=pd.pipes(p1)(1);
    TV joint=pd.pts(i).pt;
    TV pipe_dir[2]={(pd.pts(j0).pt-joint).Normalized(),(pd.pts(j1).pt-joint).Normalized()};
    if(pipe_dir[0].Cross(pipe_dir[1])(0)<0){
        std::swap(pipe_dir[0],pipe_dir[1]);
        std::swap(p0,p1);
        std::swap(j0,j1);}

    VECTOR<TV,3> w=Wedge(joint,pd.pts(j0).pt,pd.pts(j1).pt,pd.half_width,pd.unit_length);
    T edge=2*pd.half_width*pd.unit_length;
    PAIR<ARRAY<PARTICLE_ID>,ARRAY<PARTICLE_ID> > sides;
    if(jt==corner_joint)
        sides=Corner(w(0),w(0)+2*(joint-w(0)),w(0)+w(1)*edge,w(0)+w(2)*edge,pd.unit_length,
            pipe_dir[0],0.5,pipe_dir[1],0.5);
    else
        sides=Arc(w(0),w(0)+w(1)*edge,w(0)+w(2)*edge,pd.half_width,pd.unit_length,
            pipe_dir[0],0.5,pipe_dir[1],0.5);
    ARRAY<PARTICLE_ID> g0,g1;
    Weld(2*pd.half_width,sides.x,sides.y,pd.unit_length,g0,g1);

    if(pd.pipes(p0).x==i) g0.Reverse();
    if(pd.pipes(p1).x!=i) g1.Reverse();
    con.Set({i,p0},g0);
    con.Set({i,p1},g1);
}
//#####################################################################
// Function Generate_3_Joint_SmallMin
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Generate_3_Joint_SmallMin(VERTEX_ID i,const VECTOR<VERTEX_ID,3>& ends,const VECTOR<PIPE_ID,3>& pipes,
    const VECTOR<T,3>& angles,const VECTOR<VECTOR<TV,3>,3>& tri,
    const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con)
{
    TV pipe_dir[3];
    T min_angle=2*pi;
    int min_angle_idx=-1;
    TV c=pd.pts(i).pt;
    for(int j=0;j<3;j++){
        pipe_dir[j]=(pd.pts(ends(j)).pt-c).Normalized();
        if(angles(j)<min_angle){
            min_angle=angles(j);
            min_angle_idx=j;}}

    int end0=(min_angle_idx+1)%3,end1=(min_angle_idx+2)%3,end2=min_angle_idx;
    int merging_end[2]={end1,end0};
    if(angles(end0)>angles(end1)){
        end2=end0;
        end0=end1;
        end1=min_angle_idx;
        merging_end[true]=end1;
        merging_end[false]=end0;}
    bool reverse[3];
    reverse[end2]=pd.pipes(pipes(end2)).x!=i;
    reverse[merging_end[true]]=pd.pipes(pipes(merging_end[true])).x!=i;
    reverse[merging_end[false]]=pd.pipes(pipes(merging_end[false])).x==i;
    T edge=2*pd.half_width*pd.unit_length;
    VECTOR<TV,3> w=tri(min_angle_idx);
    ARRAY<TV> merging_points={w(0)+w(1)*edge,w(0),w(0)+w(2)*edge};

    TV u=tri(end0)(0)-c;
    TV v=(tri(min_angle_idx)(0)-c).Normalized();
    TV q=c+(v.Dot(u)*2*v-u);
    v=(pd.pts(ends(merging_end[false])).pt-c).Normalized();
    TV p=c+(v.Dot(u)*2*v-u);

    TV elbow=tri((merging_end[true]+1)%3)(0);
    PAIR<ARRAY<PARTICLE_ID>,ARRAY<PARTICLE_ID> > sides=Corner(tri(end0)(0),elbow,p,q,pd.unit_length,
        pipe_dir[merging_end[false]],0.5,TV(),0);
    if(merging_end[true]==end0) std::swap(sides.x,sides.y);
    ARRAY<PARTICLE_ID> f[3],left;
    Weld(2*pd.half_width,sides.x,sides.y,pd.unit_length,f[merging_end[false]],left);

    ARRAY<PARTICLE_ID> merging_side,g[2];
    for(int k=1;k<3;k++){
        v=(merging_points(k)-merging_points(k-1)).Normalized();
        for(int j=0;j<2*pd.half_width;j++){
            PARTICLE_ID base=Add_Particle();
            X(base)=merging_points(k-1)+v*j*pd.unit_length;
            g[k-1].Append(base);
            merging_side.Append(base);}}
    g[0].Append(g[1](0));
    PARTICLE_ID base=Add_Particle();
    X(base)=merging_points(2);
    g[1].Append(base);
    merging_side.Append(base);
    f[(merging_end[false]+1)%3]=g[0];
    f[(merging_end[false]+2)%3]=g[1];

    int res=(X(merging_side(0))-X(left(0))).Magnitude()/pd.unit_length;
    ARRAY<PARTICLE_ID> dummy;
    Weld(res+1,merging_side,left,pd.unit_length,dummy,dummy);

    if(reverse[end0]) f[end0].Reverse();
    if(reverse[end1]) f[end1].Reverse();
    if(reverse[end2]) f[end2].Reverse();
    con.Set({i,pipes(end0)},f[end0]);
    con.Set({i,pipes(end1)},f[end1]);
    con.Set({i,pipes(end2)},f[end2]);
}
//#####################################################################
// Function Generate_3_Joint_LargeMin
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Generate_3_Joint_LargeMin(VERTEX_ID i,const VECTOR<VERTEX_ID,3>& ends,const VECTOR<PIPE_ID,3>& pipes,
    const VECTOR<T,3>& angles,const VECTOR<VECTOR<TV,3>,3>& tri,
    const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con)
{
    TV c=pd.pts(i).pt;
    TV pipe_dir[3];
    T max_angle=-1;
    int end0=-1;
    for(int j=0;j<3;j++){
        pipe_dir[j]=(pd.pts(ends(j)).pt-c).Normalized();
        if(angles(j)>max_angle){
            end0=j;
            max_angle=angles(j);}}
    int end1=(end0+1)%3,end2=(end0+2)%3;

    TV p[2]={tri(end0)(0),tri(end2)(0)};
    TV v=p[1]-p[0];
    T proj=v.Dot(pipe_dir[end0]);
    if(proj>=0) p[0]+=pipe_dir[end0]*proj;
    else p[1]-=pipe_dir[end2]*proj;
    TV q[2]={tri(end0)(0),tri(end1)(0)};
    v=q[1]-q[0];
    proj=v.Dot(pipe_dir[end1]);
    if(proj>=0) q[0]+=pipe_dir[end1]*proj;
    else q[1]-=pipe_dir[end1]*proj;
    T ext=2*pd.unit_length;
    TV r[2]={tri(end2)(0)+pipe_dir[end2]*ext,tri(end1)(0)+pipe_dir[end2]*ext};
    v=r[1]-r[0];
    proj=v.Dot(pipe_dir[end2]);
    if(proj>=0) r[0]+=pipe_dir[end2]*proj;
    else r[1]-=pipe_dir[end2]*proj;

    ARRAY<PARTICLE_ID> side0=Polyline({r[0],p[1]},pd.unit_length);
    ARRAY<PARTICLE_ID> side1=Polyline({r[1],q[1]},pd.unit_length);
    ARRAY<PARTICLE_ID> f0,f1;
    Weld(2*pd.half_width,side0,side1,pd.unit_length,f0,f1);

    ARRAY<PARTICLE_ID> bottom=Polyline({p[0],tri(end0)(0),q[0]},pd.unit_length);
    ARRAY<PARTICLE_ID> g0,g1;
    Weld(2*pd.half_width,bottom,f1,pd.unit_length,g0,g1);

    if(pd.pipes(pipes(end0)).x==i) g0.Reverse();
    if(pd.pipes(pipes(end1)).x!=i) g1.Reverse();
    if(pd.pipes(pipes(end2)).x==i) f0.Reverse();
    con.Set({i,pipes(end0)},g0);
    con.Set({i,pipes(end1)},g1);
    con.Set({i,pipes(end2)},f0);
}
//#####################################################################
// Function Generate_3_Joint
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Generate_3_Joint(VERTEX_ID i,const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con)
{
    Dump_Input(pd);
    const ARRAY<PIPE_ID>& pipes=pd.pts(i).joints;
    PHYSBAM_ASSERT(pipes.m==3);
    VECTOR<VERTEX_ID,3> ends;
    VECTOR<PIPE_ID,3> ccw_pipes;
    for(int j=0;j<3;j++){
        ccw_pipes(j)=pipes(j);
        ends(j)=pd.pipes(pipes(j)).x;
        if(ends(j)==i) ends(j)=pd.pipes(pipes(j)).y;}
    T signed_area=(pd.pts(ends(1)).pt-pd.pts(ends(0)).pt).Cross(pd.pts(ends(2)).pt-pd.pts(ends(1)).pt)(0);
    if(signed_area<0){
        std::swap(ccw_pipes(1),ccw_pipes(2));
        std::swap(ends(1),ends(2));}

    TV c=pd.pts(i).pt;
    VECTOR<VECTOR<TV,3>,3> tri;
    VECTOR<T,3> angles;
    T min_angle=2*pi;
    for(int j=0;j<3;j++){
        int j1=(j+1)%3;
        tri(j)=Wedge(c,pd.pts(ends(j)).pt,pd.pts(ends(j1)).pt,pd.half_width,pd.unit_length);
        TV u=(pd.pts(ends(j)).pt-c).Normalized();
        TV v=(pd.pts(ends(j1)).pt-c).Normalized();
        angles(j)=acos(u.Dot(v));
        if(u.Cross(v)(0)<0) angles(j)=2*pi-angles(j);
        if(angles(j)<min_angle) min_angle=angles(j);}

    for(int j=0;j<3;j++){
        std::string s=LOG::sprintf("%f",angles(j)*180/pi);
        Add_Debug_Text(tri(j)(0),s,VECTOR<T,3>(1,1,1));}

    if(min_angle<pi/4) Generate_3_Joint_SmallMin(i,ends,ccw_pipes,angles,tri,pd,con);
    else Generate_3_Joint_LargeMin(i,ends,ccw_pipes,angles,tri,pd,con);
}
//#####################################################################
// Function Generate_Joint
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Generate_Joint(VERTEX_ID i,const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con)
{
    const ARRAY<PIPE_ID>& p=pd.pts(i).joints;
    if(p.m==1)
        Generate_End(i,p(0),pd,con);
    else if(p.m==2)
        Generate_2_Joint(i,pd,con,pd.pts(i).joint_type);
    else if(p.m==3)
        Generate_3_Joint(i,pd,con);
    else PHYSBAM_FATAL_ERROR("joint type not handled");
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Compute(const PARSE_DATA_FEM<TV,TW>& pd)
{
    CONNECTION con;
    for(VERTEX_ID i(0);i<pd.pts.m;i++){
        Generate_Joint(i,pd,con);}
    pipes.Resize(pd.pipes.m);
    for(PIPE_ID i(0);i<pd.pipes.m;i++){
        Generate_Pipe(i,pd,con);}
    area_hidden.mesh.Set_Number_Nodes(area_hidden.particles.number);
    Allocate_Dofs(pd);
    area_hidden.mesh.Initialize_Incident_Elements();
    area_hidden.mesh.Initialize_Element_Edges();
}
//#####################################################################
// Function Allocate_Dofs
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Allocate_Dofs(const PARSE_DATA_FEM<TV,TW>& pd)
{
    area_hidden.mesh.Initialize_Segment_Mesh();
    area_hidden.mesh.Initialize_Edge_Triangles();
    vel_edge_dofs.Resize(Number_Edges(),init_all,DOF_ID(-2));
    edge_blocks.Resize(Number_Edges(),init_all,BLOCK_ID(-1));

    for(EDGE_ID i(0);i<Number_Edges();i++)
        if(Edge_Triangles(i).m==1){
            auto key=Edge(i).Sorted();
            if(!bc_map.Contains(key)){
                if(BC_ID* pbc=particle_bc_map.Get_Pointer(key.x)) *pbc=pd.wall_bc;
                else particle_bc_map.Set(key.x,pd.wall_bc);
                if(BC_ID* pbc=particle_bc_map.Get_Pointer(key.y)) *pbc=pd.wall_bc;
                else particle_bc_map.Set(key.y,pd.wall_bc);
                bc_map.Set({key.x,key.y},pd.wall_bc);}}

    ARRAY<int,BLOCK_ID> block_vel_dofs(blocks.m),block_pressure_dofs(blocks.m);
    for(PARTICLE_ID i(0);i<node_blocks.m;i++){
        block_pressure_dofs(node_blocks(i))++;
        BC_ID* bc_idx=particle_bc_map.Get_Pointer(i);
        if(!bc_idx || pd.bc(*bc_idx).type!=dirichlet_v)
            block_vel_dofs(node_blocks(i))++;}
    for(EDGE_ID ei(0);ei<Number_Edges();ei++){
        BC_ID* bc_idx=bc_map.Get_Pointer(Edge(ei).Sorted());
        if(!bc_idx || pd.bc(*bc_idx).type!=dirichlet_v){
            auto neighbor_tris=Edge_Triangles(ei);
            PHYSBAM_ASSERT(neighbor_tris.m==1 || neighbor_tris.m==2);
            BLOCK_ID bid=elem_data(neighbor_tris(0)).block_id;
            if(neighbor_tris.m==2 && elem_data(neighbor_tris(1)).greed>elem_data(neighbor_tris(0)).greed)
                bid=elem_data(neighbor_tris(1)).block_id;
            edge_blocks(ei)=bid;
            block_vel_dofs(bid)++;}}
    for(BLOCK_ID i(1);i<blocks.m;i++){
        block_vel_dofs(i)+=block_vel_dofs(i-1);
        block_pressure_dofs(i)+=block_pressure_dofs(i-1);}
    int num_pressure_dofs=block_pressure_dofs.Last();
    int num_vel_dofs=block_vel_dofs.Last()*TV::m;
    num_dofs=DOF_ID(num_pressure_dofs+num_vel_dofs);

    for(BLOCK_ID i(0);i<blocks.m;i++) blocks(i).num_dofs=0;
    dof_map.Resize(num_dofs,use_init,{-1,-1});
    pressure_dofs.Resize(Number_Particles(),use_init,DOF_ID(-1));
    for(PARTICLE_ID i(0);i<Number_Particles();i++){
        BLOCK_ID bid=node_blocks(i);
        pressure_dofs(i)=DOF_ID(num_vel_dofs+--block_pressure_dofs(bid));
        dof_map(pressure_dofs(i))={Value(bid),blocks(bid).num_dofs++};}

    vel_node_dofs.Resize(Number_Particles(),use_init,DOF_ID(-2));
    for(PARTICLE_ID i(0);i<Number_Particles();i++){
        BC_ID* bc_idx=particle_bc_map.Get_Pointer(i);
        if(bc_idx && pd.bc(*bc_idx).type==dirichlet_v) continue;
        BLOCK_ID bid=node_blocks(i);
        vel_node_dofs(i)=DOF_ID(--block_vel_dofs(bid)*TV::m);
        for(int a=0;a<TV::m;a++)
            dof_map(vel_node_dofs(i)+a)={Value(bid),blocks(bid).num_dofs++};}

    vel_edge_dofs.Resize(Number_Edges(),use_init,DOF_ID(-2));
    for(EDGE_ID i(0);i<Number_Edges();i++){
        BC_ID* bc_idx=bc_map.Get_Pointer(Edge(i).Sorted());
        if(bc_idx && pd.bc(*bc_idx).type==dirichlet_v) continue;
        BLOCK_ID bid=edge_blocks(i);
        vel_edge_dofs(i)=DOF_ID(--block_vel_dofs(bid)*TV::m);
        for(int a=0;a<TV::m;a++)
            dof_map(vel_edge_dofs(i)+a)={Value(bid),blocks(bid).num_dofs++};}
}
//#####################################################################
// Function Print_Statistics
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Print_Statistics() const
{
    area_hidden.mesh.Initialize_Adjacent_Elements();
    LOG::printf("orientation consistent: %d\n",area_hidden.mesh.Orientations_Consistent());

    PARTICLE_HIERARCHY<TV> ph(area_hidden.particles.X);
    VECTOR<PARTICLE_ID,2> pair=Triangle(TRIANGLE_ID()).Remove_Index(2);
    T min_dist=(X(pair.x)-X(pair.y)).Magnitude();
    for(PARTICLE_ID j(0);j<Number_Particles();j++){
        ARRAY<int> q;
        ph.Intersection_List(X(j),q,min_dist);
        for(int i=0;i<q.m;i++){
            if(PARTICLE_ID(q(i))==j) continue;
            T d=(X(PARTICLE_ID(q(i)))-X(j)).Magnitude();
            if(d<min_dist){
                pair.x=j;
                pair.y=PARTICLE_ID(q(i));
                min_dist=d;}}}
    LOG::printf("min dist: %f\n",min_dist);
}
//#####################################################################
// Function Dump_Mesh
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Dump_Mesh() const
{
    for(TRIANGLE_ID i(0);i<Number_Triangles();i++){
        auto tri=Triangle(i);
        for(int j=0;j<3;j++){
            PARTICLE_ID v0=tri(j),v1=tri((j+1)%3);
            Add_Debug_Object<TV,2>(VECTOR<TV,2>(X(v0),X(v1)),VECTOR<T,3>(0.5,0.5,0.5));}}
}
//#####################################################################
// Function Dump_Edge_Blocks
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Dump_Edge_Blocks() const
{
    Dump_Mesh();
    VECTOR<T,3> colors[]={VECTOR<T,3>(1,0,1),VECTOR<T,3>(0,1,1)};
    for(EDGE_ID i(0);i<Number_Edges();i++){
        auto e=Edge(i);
        TV p=0.5*(X(e.x)+X(e.y));
        std::string s=LOG::sprintf("%i",edge_blocks(i));
        Add_Debug_Text(p,s,colors[Value(edge_blocks(i))%2]);}
}
//#####################################################################
// Function Dump_Node_Blocks
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Dump_Node_Blocks() const
{
    Dump_Mesh();
    VECTOR<T,3> colors[]={VECTOR<T,3>(1,0,1),VECTOR<T,3>(0,1,1)};
    for(PARTICLE_ID i(0);i<node_blocks.m;i++){
        if(node_blocks(i)<BLOCK_ID()) continue;
        std::string s=LOG::sprintf("%i",node_blocks(i));
        Add_Debug_Text(X(i),s,colors[Value(node_blocks(i))%2]);}
}
//#####################################################################
// Function Dump_Dofs
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Dump_Dofs() const
{
    VECTOR<T,3> colors[]={VECTOR<T,3>(1,0,1),VECTOR<T,3>(0,1,1)};
    Dump_Mesh();
    for(PARTICLE_ID i(0);i<Number_Particles();i++){
        if(vel_node_dofs(i)<DOF_ID()) continue;
        std::string s=LOG::sprintf("%i",vel_node_dofs(i));
        Add_Debug_Text(X(i),s,colors[Value(node_blocks(i))%2]);}
    for(EDGE_ID i(0);i<Number_Edges();i++){
        if(vel_edge_dofs(i)<DOF_ID()) continue;
        auto e=Edge(i);
        TV p=0.5*(X(e.x)+X(e.y));
        std::string s=LOG::sprintf("%i",vel_edge_dofs(i));
        BLOCK_ID bid=edge_blocks(i);
        Add_Debug_Text(p,s,colors[Value(bid)%2]);}
    Flush_Frame("vel dofs");

    Dump_Mesh();
    for(PARTICLE_ID i(0);i<Number_Particles();i++){
        if(pressure_dofs(i)<DOF_ID()) continue;
        std::string s=LOG::sprintf("%i",pressure_dofs(i));
        Add_Debug_Text(X(i),s,colors[Value(node_blocks(i))%2]);}
    Flush_Frame("pressure dofs");
}
//#####################################################################
// Function Dump_Layout
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Dump_Layout() const
{
    Dump_Mesh();
    VECTOR<T,3> reg(0,1,0),irreg(1,0,0);
    for(TRIANGLE_ID i(0);i<Number_Triangles();i++){
        std::string s=LOG::sprintf("%i",elem_data(i).block_id);
        Add_Debug_Text(area_hidden.particles.X.Subset(area_hidden.mesh.elements(Value(i))).Average(),s,
            blocks(elem_data(i).block_id).pipe_id>=PIPE_ID()?reg:irreg);}
}
//#####################################################################
// Function Dump_Input
//#####################################################################
template<class T> template<class TW> void FLUID_LAYOUT_FEM<VECTOR<T,2> >::
Dump_Input(const PARSE_DATA_FEM<TV,TW>& pd) const
{
    for(auto& i:pd.pts)
        Add_Debug_Particle(i.pt,VECTOR<T,3>(0,0,1));
    for(auto& i:pd.pipes)
        Add_Debug_Object<TV,2>(VECTOR<TV,2>(pd.pts(i.x).pt,pd.pts(i.y).pt),VECTOR<T,3>(0,0,1));
}
template struct FLUID_LAYOUT_FEM<VECTOR<double,2> >;
template void FLUID_LAYOUT_FEM<VECTOR<double,2> >::Compute<VECTOR<double,2> >(PARSE_DATA_FEM<VECTOR<double,2>,VECTOR<double,2> > const&);
template void FLUID_LAYOUT_FEM<VECTOR<double,2> >::Dump_Input<VECTOR<double,2> >(PARSE_DATA_FEM<VECTOR<double,2>,VECTOR<double,2> > const&) const;
template void FLUID_LAYOUT_FEM<VECTOR<double,2> >::Generate_Joint<VECTOR<double,3> >(VERTEX_ID,
    PARSE_DATA_FEM<VECTOR<double,2>,VECTOR<double,3> > const&,HASHTABLE<PAIR<VERTEX_ID,PIPE_ID>,ARRAY<PARTICLE_ID,int> >&);
template void FLUID_LAYOUT_FEM<VECTOR<double,2> >::Generate_Pipe<VECTOR<double,3> >(PIPE_ID,
    PARSE_DATA_FEM<VECTOR<double,2>,VECTOR<double,3> > const&,HASHTABLE<PAIR<VERTEX_ID,PIPE_ID>,ARRAY<PARTICLE_ID,int> > const&);
}
