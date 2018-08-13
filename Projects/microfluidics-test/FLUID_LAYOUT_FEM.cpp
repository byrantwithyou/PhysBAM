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
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <list>
#include <map>
#include "FLUID_LAYOUT_FEM.h"

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_LAYOUT_FEM<TV>::
FLUID_LAYOUT_FEM(): area(*new TRIANGULATED_AREA<T>)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_LAYOUT_FEM<TV>::
~FLUID_LAYOUT_FEM()
{
    delete &area;
}
//#####################################################################
// Function Generate_Pipe
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Pipe(int pipe,const PARSE_DATA_FEM<TV>& pd,const CONNECTION& con)
{
    typedef VECTOR<int,3> E;
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    int base=particles.number;
    const ARRAY<CONNECTION_DATA>* f0=con.Get_Pointer({pd.pipes(pipe).x,pipe}),*f1=con.Get_Pointer({pd.pipes(pipe).y,pipe});
    if(!f0 || !f1) return;
    TV v0=particles.X((*f0)(pd.half_width).pid),v1=particles.X((*f1)(pd.half_width).pid);
    TV d=v1-v0;
    T l=d.Normalize();
    TV t(-d.y,d.x);
    int width=2*pd.half_width+1;
    int height=(int)(l/pd.unit_length)+1;
    T min_percent=0.5;
    bool irregular_last=false;
    T rem=l-(height-1)*pd.unit_length;
    if(rem>0){
        irregular_last=true;
        if(rem>min_percent*pd.unit_length)
            height++;}
    if(height>2)
        particles.Add_Elements(width*(height-2));
    auto pid=[base,width,height,&pd,f0,f1](int i,int j)
    {
        if(i==0) return (*f0)(pd.half_width+j).pid;
        else if(i==height-1) return (*f1)(pd.half_width+j).pid;
        else return base+(i-1)*width+j+pd.half_width;
    };
    auto is_collapsed=[&pd,height,f0,f1](int i,int j)
    {
        if(i==0) return (*f0)(pd.half_width+j).collapsed;
        else if(i==height-1) return (*f1)(pd.half_width+j).collapsed;
        else return false;
    };
    auto assign_edge=[this](int n0,int n1,int bid,bool overwrite)
    {
        if(n0>n1) std::swap(n0,n1);
        int* b=edge_blocks.Get_Pointer({n0,n1});
        if(b && overwrite) *b=bid;
        else edge_blocks.Set({n0,n1},bid);
    };
    auto assign_node=[this](int n,int bid,bool overwrite)
    {
        if(node_blocks_assigned(n)==false || overwrite){
            node_blocks(n)=bid;
            node_blocks_assigned(n)=true;}
    };
    TV next=pd.unit_length*d;;
    bool regular=true;
    for(int i=1;i<height;i++){
        int num_left_node=pd.half_width+1,num_right_node=pd.half_width+1;
        int num_left_edge=pd.half_width,num_right_edge=pd.half_width;
        for(int j=-pd.half_width;j<=pd.half_width;j++){
            if(i!=height-1)
                particles.X(pid(i,j))=v0+next+j*pd.unit_length*t;
            assign_node(pid(i-1,j),blocks.m,num_left_node--<=0);
            assign_node(pid(i,j),blocks.m,num_right_node-->0);
            if(j!=pd.half_width){
                regular=regular && !is_collapsed(i,j) && !is_collapsed(i-1,j+1) && !is_collapsed(i-1,j);
                area.mesh.elements.Append(E(pid(i,j),pid(i-1,j+1),pid(i-1,j)));
                assign_edge(pid(i-1,j),pid(i,j),blocks.m,false);
                assign_edge(pid(i-1,j+1),pid(i,j),blocks.m,false);
                assign_edge(pid(i-1,j+1),pid(i-1,j),blocks.m,num_left_edge--<=0);
                elem_data.Append({blocks.m});}
            if(j!=-pd.half_width){
                regular=regular && !is_collapsed(i,j) && !is_collapsed(i-1,j) && !is_collapsed(i,j-1);
                area.mesh.elements.Append(E(pid(i,j),pid(i-1,j),pid(i,j-1)));
                assign_edge(pid(i-1,j),pid(i,j),blocks.m,false);
                assign_edge(pid(i-1,j),pid(i,j-1),blocks.m,false);
                assign_edge(pid(i,j),pid(i,j-1),blocks.m,num_right_edge-->0);
                elem_data.Append({blocks.m});}}
        if(i!=0) blocks.Append({regular});
        if(i==height-2){
            if(irregular_last) regular=false;
            next=l*d;}
        else {
            regular=true;
            next=(i+1)*pd.unit_length*d;}}
}
//#####################################################################
// Function Mark_BC
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Mark_BC(const ARRAY<CONNECTION_DATA>& pindices,BC_TYPE bc_type)
{
    int bc_index=bc.Append({bc_type});
    for(int i=0;i<pindices.m;i++)
        bc_map.Set(pindices(i).pid,bc_index);
}
//#####################################################################
// Function Generate_End
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_End(int i,int pipe,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    int j=pd.pipes(pipe).x;
    if(j==i) j=pd.pipes(pipe).y;
    TV p=pd.pts(i).pt;
    TV v=pd.pts(j).pt-pd.pts(i).pt;
    TV u=TV(-v.y,v.x).Normalized();
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    int base=particles.Add_Elements(2*pd.half_width+1);
    ARRAY<CONNECTION_DATA>& indices=con.Get_Or_Insert({i,pipe});
    for(int k=-pd.half_width;k<=pd.half_width;k++){
        int pid=base+k+pd.half_width;
        indices.Append({pid,false});
        particles.X(pid)=p+k*pd.unit_length*u;}
    if(pd.pipes(pipe).x!=i) indices.Reverse();
    if(pd.pts(i).bc_type!=nobc) Mark_BC(indices,pd.pts(i).bc_type);
}
//#####################################################################
// Function Wedge
//#####################################################################
template<class TV> VECTOR<TV,3> FLUID_LAYOUT_FEM<TV>::
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
template<class TV> ARRAY<int> FLUID_LAYOUT_FEM<TV>::
Sample_Interpolated(T s,const ARRAY<int>& side0,const ARRAY<int>& side1,T unit_length)
{
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    ARRAY<T> l0(side0.m),l1(side1.m);
    l0(0)=0;
    for(int j=1;j<side0.m;j++)
        l0(j)=l0(j-1)+(particles.X(side0(j))-particles.X(side0(j-1))).Magnitude();
    l1(0)=0;
    for(int j=1;j<side1.m;j++)
        l1(j)=l1(j-1)+(particles.X(side1(j))-particles.X(side1(j-1))).Magnitude();
    auto point=[&particles](int j,T t,const ARRAY<int>& side)
    {
        if(j>=side.m-1) return particles.X(side(j));
        TV p=particles.X(side(j));
        return p+t*(particles.X(side(j+1))-p);
    };
    auto loc=[point](T lambda,const ARRAY<int>& side,const ARRAY<T>& l)
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
    int base=particles.Add_Elements(verts.size());
    ARRAY<int> indices(verts.size());
    int j=0;
    for(auto iter=verts.begin();iter!=verts.end();iter++){
        particles.X(base+j)=iter->y;
        indices(j)=base+j;
        j++;}
    return indices;
}
//#####################################################################
// Function Merge_Interpolated
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Merge_Interpolated(const ARRAY<int>& left,const ARRAY<int>& right)
{
    typedef VECTOR<int,3> E;
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    int i=0,j=0,alt=0;
    auto angle=[&particles](int v0,int v1,int v2)
    {
        TV u=(particles.X(v2)-particles.X(v1)).Normalized();
        TV v=(particles.X(v0)-particles.X(v1)).Normalized();
        return acos(u.Dot(v));
    };
    ARRAY<E> elems;
    auto assign_edge=[this](int n0,int n1,int bid,bool overwrite)
    {
        if(n0>n1) std::swap(n0,n1);
        int* b=edge_blocks.Get_Pointer({n0,n1});
        if(b && overwrite) *b=bid;
        else edge_blocks.Set({n0,n1},bid);
    };
    auto assign_node=[this](int n,int bid,bool overwrite)
    {
        if(node_blocks_assigned(n)==false || overwrite){
            node_blocks(n)=bid;
            node_blocks_assigned(n)=true;}
    };
    int num_left_edge=(left.m-1)/2,num_right_edge=(right.m-1)/2;
    int num_left_node=left.m/2,num_right_node=right.m/2;
    assign_node(left(i),blocks.m,num_left_node--<=0);
    assign_node(right(j),blocks.m,num_right_node-->0);
    while(i<left.m-1 || j<right.m-1){
        T a0=0;
        if(i+1<left.m) a0=angle(right(j),left(i+1),left(i));
        T a1=0;
        if(j+1<right.m) a1=angle(right(j),right(j+1),left(i));
        assign_edge(left(i),right(j),blocks.m,false);
        if(j+1>=right.m || (i+1<left.m && abs(a0-a1)<1e-6 && alt==0) || (i+1<left.m && a0>a1)){
            area.mesh.elements.Append(E(left(i+1),left(i),right(j)));
            assign_edge(left(i+1),right(j),blocks.m,false);
            assign_edge(left(i+1),left(i),blocks.m,num_left_edge--<=0);
            assign_node(left(i+1),blocks.m,num_left_node--<=0);
            i++;
            alt=1;}
        else{
            area.mesh.elements.Append(E(right(j+1),left(i),right(j)));
            assign_edge(right(j+1),left(i),blocks.m,false);
            assign_edge(right(j+1),right(j),blocks.m,num_right_edge-->0);
            assign_node(right(j+1),blocks.m,num_right_node-->0);
            j++;
            alt=0;}
        elem_data.Append({blocks.m});}
}
//#####################################################################
// Function Polyline
//#####################################################################
template<class TV> ARRAY<int> FLUID_LAYOUT_FEM<TV>::
Polyline(const ARRAY<TV>& points,T unit_length)
{
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    ARRAY<int> verts;
    int base;
    for(int j=1;j<points.m;j++){
        TV v=points(j)-points(j-1);
        T l=v.Normalize();
        for(int k=0;(k+1)*unit_length<=l+0.5*unit_length;k++){
            base=particles.Add_Element();
            particles.X(base)=points(j-1)+v*k*unit_length;
            verts.Append(base);}}
    base=particles.Add_Element();
    particles.X(base)=points.Last();
    verts.Append(base);
    return verts;
}
//#####################################################################
// Function Arc
//#####################################################################
template<class TV> PAIR<ARRAY<int>,ARRAY<int> > FLUID_LAYOUT_FEM<TV>::
Arc(const TV& c,const TV& p0,const TV& p1,int half_width,T unit_length,
    const TV& dir0,T n0,const TV& dir1,T n1)
{
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    ARRAY<int> side0,side1;
    int base;
    if(n0){
        base=particles.Add_Elements(2);
        particles.X(base)=c+dir0*unit_length*n0;
        particles.X(base+1)=p0+dir0*unit_length*n0;
        side0.Append(base);
        side1.Append(base+1);}
    base=particles.Add_Element();
    particles.X(base)=c;
    side0.Append(base);
    T unit_angle=2*asin((T)1/(half_width*4));
    TV v0=(p0-c).Normalized(),v1=(p1-c).Normalized();
    T total=acos(v0.Dot(v1));
    for(int j=0;(j+1)*unit_angle<=total+0.5*unit_angle;j++){
        T a=-j*unit_angle;
        TV v(cos(a)*v0(0)-sin(a)*v0(1),sin(a)*v0(0)+cos(a)*v0(1));
        base=particles.Add_Element();
        particles.X(base)=c+v*2*half_width*unit_length;
        side1.Append(base);}
    base=particles.Add_Element();
    particles.X(base)=p1;
    side1.Append(base);
    if(n1){
        base=particles.Add_Elements(2);
        particles.X(base)=c+dir1*unit_length*n1;
        particles.X(base+1)=p1+dir1*unit_length*n1;
        side0.Append(base);
        side1.Append(base+1);}
    return {side0,side1};
}
//#####################################################################
// Function Corner
//#####################################################################
template<class TV> PAIR<ARRAY<int>,ARRAY<int> > FLUID_LAYOUT_FEM<TV>::
Corner(const TV& c,const TV& elbow,const TV& p0,const TV& p1,T unit_length,
    const TV& dir0,T n0,const TV& dir1,T n1)
{
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    ARRAY<int> side0,side1;
    int base;
    if(n0){
        base=particles.Add_Elements(2);
        particles.X(base)=c+dir0*unit_length*n0;
        particles.X(base+1)=p0+dir0*unit_length*n0;
        side0.Append(base);
        side1.Append(base+1);}
    base=particles.Add_Element();
    particles.X(base)=c;
    side0.Append(base);
    TV start=p0;
    TV v=elbow-start;
    T total=v.Normalize();
    for(int j=0;(j+1)*unit_length<=total+0.5*unit_length;j++){
        base=particles.Add_Element();
        particles.X(base)=start+v*j*unit_length;
        side1.Append(base);}
    TV end=p1;
    v=end-elbow;
    total=v.Normalize();
    for(int j=0;(j+1)*unit_length<=total+0.5*unit_length;j++){
        base=particles.Add_Element();
        particles.X(base)=elbow+v*j*unit_length;
        side1.Append(base);}
    base=particles.Add_Element();
    particles.X(base)=end;
    side1.Append(base);
    if(n1){
        base=particles.Add_Elements(2);
        particles.X(base)=c+dir1*unit_length*n1;
        particles.X(base+1)=p1+dir1*unit_length*n1;
        side0.Append(base);
        side1.Append(base+1);}
    return {side0,side1};
}
//#####################################################################
// Function Weld
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Weld(int n,const ARRAY<int>& side0,const ARRAY<int>& side1,T unit_length,ARRAY<int>& f0,ARRAY<int>& f1)
{
    f0.Append(side0(0));
    f1.Append(side0.Last());
    ARRAY<int> prev=side0;
    for(int j=1;j<n;j++){
        T s=(T)j/n;
        ARRAY<int> cur=Sample_Interpolated(s,side0,side1,unit_length);
        Merge_Interpolated(cur,prev);
        blocks.Append({false});
        f0.Append(cur(0));
        f1.Append(cur.Last());
        prev=cur;}
    Merge_Interpolated(side1,prev);
    blocks.Append({false});
    f0.Append(side1(0));
    f1.Append(side1.Last());
}
//#####################################################################
// Function Generate_2_Joint
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_2_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con,JOINT_TYPE jt)
{
    typedef VECTOR<int,3> E;
    const ARRAY<int>* pipes=pd.joints.Get_Pointer(i);
    PHYSBAM_ASSERT(pipes);
    PHYSBAM_ASSERT(pipes->m==2);
    int p0=(*pipes)(0),p1=(*pipes)(1);
    int j0=pd.pipes(p0)(0),j1=pd.pipes(p1)(0);
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
    PAIR<ARRAY<int>,ARRAY<int> > sides;
    if(jt==corner_joint)
        sides=Corner(w(0),w(0)+2*(joint-w(0)),w(0)+w(1)*edge,w(0)+w(2)*edge,pd.unit_length,
            pipe_dir[0],0.5,pipe_dir[1],0.5);
    else
        sides=Arc(w(0),w(0)+w(1)*edge,w(0)+w(2)*edge,pd.half_width,pd.unit_length,
            pipe_dir[0],0.5,pipe_dir[1],0.5);
    ARRAY<int> g0,g1;
    Weld(2*pd.half_width,sides.x,sides.y,pd.unit_length,g0,g1);

    ARRAY<CONNECTION_DATA> f0,f1;
    for(int j=0;j<g0.m;j++){
        f0.Append({g0(j),false});
        f1.Append({g1(j),false});}
    if(pd.pipes(p0).x==i) f0.Reverse();
    if(pd.pipes(p1).x!=i) f1.Reverse();
    con.Set({i,p0},f0);
    con.Set({i,p1},f1);
}
//#####################################################################
// Function Generate_3_Joint_SmallMin
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_3_Joint_SmallMin(int i,const VECTOR<int,3>& ends,const VECTOR<int,3>& pipes,
    const VECTOR<T,3>& angles,const VECTOR<VECTOR<TV,3>,3>& tri,
    const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
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
    PAIR<ARRAY<int>,ARRAY<int> > sides=Corner(tri(end0)(0),elbow,p,q,pd.unit_length,
        pipe_dir[merging_end[false]],0.5,TV(),0);
    if(merging_end[true]==end0) std::swap(sides.x,sides.y);
    ARRAY<int> f[3],left;
    Weld(2*pd.half_width,sides.x,sides.y,pd.unit_length,f[merging_end[false]],left);

    ARRAY<int> merging_side,g[2];
    for(int k=1;k<3;k++){
        v=merging_points(k)-merging_points(k-1);
        T h=v.Normalize()/(2*pd.half_width);
        for(int j=0;j<2*pd.half_width;j++){
            int base=particles.Add_Element();
            particles.X(base)=merging_points(k-1)+v*j*pd.unit_length;
            g[k-1].Append(base);
            merging_side.Append(base);}}
    g[0].Append(g[1](0));
    int base=particles.Add_Element();
    particles.X(base)=merging_points(2);
    g[1].Append(base);
    merging_side.Append(base);
    f[(merging_end[false]+1)%3]=g[0];
    f[(merging_end[false]+2)%3]=g[1];

    int res=(particles.X(merging_side(0))-particles.X(left(0))).Magnitude()/pd.unit_length;
    ARRAY<int> dummy;
    Weld(res+1,merging_side,left,pd.unit_length,dummy,dummy);

    if(reverse[end0]) f[end0].Reverse();
    if(reverse[end1]) f[end1].Reverse();
    if(reverse[end2]) f[end2].Reverse();
    ARRAY<CONNECTION_DATA> tmp[3];
    for(int j=0;j<2*pd.half_width+1;j++){
        tmp[end0].Append({f[end0](j),false});
        tmp[end1].Append({f[end1](j),false});
        tmp[end2].Append({f[end2](j),false});}
    con.Set({i,pipes(end0)},tmp[end0]);
    con.Set({i,pipes(end1)},tmp[end1]);
    con.Set({i,pipes(end2)},tmp[end2]);
}
//#####################################################################
// Function Generate_3_Joint_LargeMin
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_3_Joint_LargeMin(int i,const VECTOR<int,3>& ends,const VECTOR<int,3>& pipes,
    const VECTOR<T,3>& angles,const VECTOR<VECTOR<TV,3>,3>& tri,
    const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
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

    ARRAY<int> side0=Polyline({r[0],p[1]},pd.unit_length);
    ARRAY<int> side1=Polyline({r[1],q[1]},pd.unit_length);
    ARRAY<int> f0,f1;
    Weld(2*pd.half_width,side0,side1,pd.unit_length,f0,f1);

    ARRAY<int> bottom=Polyline({p[0],tri(end0)(0),q[0]},pd.unit_length);
    ARRAY<int> g0,g1;
    Weld(2*pd.half_width,bottom,f1,pd.unit_length,g0,g1);

    if(pd.pipes(pipes(end0)).x==i) g0.Reverse();
    if(pd.pipes(pipes(end1)).x!=i) g1.Reverse();
    if(pd.pipes(pipes(end2)).x==i) f0.Reverse();
    ARRAY<CONNECTION_DATA> tmp[3];
    for(int j=0;j<2*pd.half_width+1;j++){
        tmp[end0].Append({g0(j),false});
        tmp[end1].Append({g1(j),false});
        tmp[end2].Append({f0(j),false});}
    con.Set({i,pipes(end0)},tmp[end0]);
    con.Set({i,pipes(end1)},tmp[end1]);
    con.Set({i,pipes(end2)},tmp[end2]);
}
//#####################################################################
// Function Generate_3_Joint
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_3_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    Dump_Input(pd);
    const ARRAY<int>* pipes=pd.joints.Get_Pointer(i);
    PHYSBAM_ASSERT(pipes && pipes->m==3);
    VECTOR<int,3> ends,ccw_pipes;
    for(int j=0;j<3;j++){
        ccw_pipes(j)=(*pipes)(j);
        ends(j)=pd.pipes((*pipes)(j)).x;
        if(ends(j)==i) ends(j)=pd.pipes((*pipes)(j)).y;}
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
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    const ARRAY<int>* p=pd.joints.Get_Pointer(i);
    PHYSBAM_ASSERT(p);
    if(p->m==1)
        Generate_End(i,(*p)(0),pd,con);
    else if(p->m==2)
        Generate_2_Joint(i,pd,con,pd.pts(i).joint_type);
    else if(p->m==3)
        Generate_3_Joint(i,pd,con);
    else PHYSBAM_FATAL_ERROR("joint type not handled");
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Compute(const PARSE_DATA_FEM<TV>& pd)
{
    area.particles.Add_Array("node_blocks",&node_blocks);
    area.particles.Add_Array("node_blocks_assigned",&node_blocks_assigned);
    CONNECTION con;
    for(int i=0;i<pd.pts.m;i++){
        Generate_Joint(i,pd,con);}
    for(int i=0;i<pd.pipes.m;i++){
        Generate_Pipe(i,pd,con);}
    area.mesh.Set_Number_Nodes(area.particles.number);
}
//#####################################################################
// Function Print_Statistics
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Print_Statistics() const
{
    area.mesh.Initialize_Adjacent_Elements();
    LOG::printf("orientation consistent: %d\n",area.mesh.Orientations_Consistent());

    const GEOMETRY_PARTICLES<TV>& particles=area.particles;
    PARTICLE_HIERARCHY<TV> ph(particles.X);
    PAIR<int,int> pair(area.mesh.elements(0)(0),area.mesh.elements(0)(1));
    T min_dist=(particles.X(pair.x)-particles.X(pair.y)).Magnitude();
    for(int j=0;j<particles.number;j++){
        ARRAY<int> q;
        ph.Intersection_List(particles.X(j),q,min_dist);
        for(int i=0;i<q.m;i++){
            if(q(i)==j) continue;
            T d=(particles.X(q(i))-particles.X(j)).Magnitude();
            if(d<min_dist){
                pair.x=j;
                pair.y=q(i);
                min_dist=d;}}}
    LOG::printf("min dist: %f\n",min_dist);
}
//#####################################################################
// Function Dump_Mesh
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Dump_Mesh() const
{
    const GEOMETRY_PARTICLES<TV>& particles=area.particles;
    //for(int i=0;i<particles.number;i++)
    //    Add_Debug_Particle(particles.X(i),VECTOR<T,3>(1,1,1));
    VECTOR<T,3> default_edge(0.5,0.5,0.5),dirichlet_bc(1,0,1),traction_bc(0,1,1);
    for(int i=0;i<area.mesh.elements.m;i++){
        //Add_Debug_Object(VECTOR<TV,3>(particles.X.Subset(area.mesh.elements(i))),VECTOR<T,3>(1,1,1));
        //VECTOR<TV,3> tri(particles.X.Subset(area.mesh.elements(i)));
        //T a=(tri(2)-tri(1)).Cross(tri(0)-tri(1))(0);
        //Add_Debug_Object(tri,a<0?VECTOR<T,3>(0,1,1):default_edge);
        for(int j=0;j<3;j++){
            VECTOR<T,3> color=default_edge;
            int v0=area.mesh.elements(i)(j),v1=area.mesh.elements(i)((j+1)%3);
            if(const int* bc0=bc_map.Get_Pointer(v0)){
                const int* bc1=bc_map.Get_Pointer(v1);
                if(bc1 && *bc1==*bc0)
                    color=bc(*bc0).bc_type==dirichlet_v?dirichlet_bc:traction_bc;}
            Add_Debug_Object<TV,2>(VECTOR<TV,2>(particles.X(v0),particles.X(v1)),color);}}
}
//#####################################################################
// Function Dump_Edge_Blocks
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Dump_Edge_Blocks() const
{
    const GEOMETRY_PARTICLES<TV>& particles=area.particles;
    Dump_Mesh();
    ARRAY<VECTOR<int,2> > keys;
    ARRAY<int> data;
    edge_blocks.Get_Keys(keys);
    edge_blocks.Get_Data(data);
    VECTOR<T,3> colors[]={VECTOR<T,3>(1,0,1),VECTOR<T,3>(0,1,1)};
    for(int i=0;i<keys.m;i++){
        VECTOR<TV,2> edge(particles.X.Subset(keys(i)));
        TV p=0.5*(edge(0)+edge(1));
        std::string s=LOG::sprintf("%i",data(i));
        Add_Debug_Text(p,s,colors[data(i)%2]);}
}
//#####################################################################
// Function Dump_Node_Blocks
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Dump_Node_Blocks() const
{
    Dump_Mesh();
    VECTOR<T,3> colors[]={VECTOR<T,3>(1,0,1),VECTOR<T,3>(0,1,1)};
    for(int i=0;i<node_blocks.m;i++){
        if(node_blocks(i)<0) continue;
        std::string s=LOG::sprintf("%i",node_blocks(i));
        Add_Debug_Text(area.particles.X(i),s,colors[node_blocks(i)%2]);}
}
//#####################################################################
// Function Dump_Layout
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Dump_Layout() const
{
    const GEOMETRY_PARTICLES<TV>& particles=area.particles;
    Dump_Mesh();
    VECTOR<T,3> reg(0,1,0),irreg(1,0,0);
    for(int i=0;i<area.mesh.elements.m;i++){
        std::string s=LOG::sprintf("%i",elem_data(i).block_id);
        Add_Debug_Text(particles.X.Subset(area.mesh.elements(i)).Average(),s,
            blocks(elem_data(i).block_id).regular?reg:irreg);}
}
//#####################################################################
// Function Dump_Input
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Dump_Input(const PARSE_DATA_FEM<TV>& pd) const
{
    for(auto& i:pd.pts)
        Add_Debug_Particle(i.pt,VECTOR<T,3>(0,0,1));
    for(auto& i:pd.pipes)
        Add_Debug_Object<TV,2>(VECTOR<TV,2>(pd.pts(i.x).pt,pd.pts(i.y).pt),VECTOR<T,3>(0,0,1));
}
template struct FLUID_LAYOUT_FEM<VECTOR<double,2> >;
}
