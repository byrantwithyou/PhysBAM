//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
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
    TV next=pd.unit_length*d;;
    bool regular=true;
    for(int i=1;i<height;i++){
        for(int j=-pd.half_width;j<=pd.half_width;j++){
            if(i!=height-1)
                particles.X(pid(i,j))=v0+next+j*pd.unit_length*t;
            if(j!=pd.half_width){
                regular=regular && !is_collapsed(i,j) && !is_collapsed(i-1,j+1) && !is_collapsed(i-1,j);
                area.mesh.elements.Append(E(pid(i,j),pid(i-1,j+1),pid(i-1,j)));
                elem_data.Append({blocks.m});}
            if(j!=-pd.half_width){
                regular=regular && !is_collapsed(i,j) && !is_collapsed(i-1,j) && !is_collapsed(i,j-1);
                area.mesh.elements.Append(E(pid(i,j),pid(i-1,j),pid(i,j-1)));
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
        return VECTOR<TV,3>(joint+unit_length*v,e0,e1);}
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
    auto loc=[&particles,point](T lambda,const ARRAY<int>& side,const ARRAY<T>& l)
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
    const GEOMETRY_PARTICLES<TV>& particles=area.particles;
    int i=0,j=0,alt=0;
    auto angle=[&particles](int v0,int v1,int v2)
    {
        TV u=(particles.X(v2)-particles.X(v1)).Normalized();
        TV v=(particles.X(v0)-particles.X(v1)).Normalized();
        return acos(u.Dot(v));
    };
    while(i<left.m-1 || j<right.m-1){
        T a0=0;
        if(i+1<left.m) a0=angle(right(j),left(i+1),left(i));
        T a1=0;
        if(j+1<right.m) a1=angle(right(j),right(j+1),left(i));
        if(j+1>=right.m || (i+1<left.m && abs(a0-a1)<1e-6 && alt==0) || (i+1<left.m && a0>a1)){
            area.mesh.elements.Append(E(left(i+1),left(i),right(j)));
            i++;
            alt=1;}
        else{
            area.mesh.elements.Append(E(right(j+1),left(i),right(j)));
            j++;
            alt=0;}
        elem_data.Append({blocks.m});}
}
//#####################################################################
// Function Generate_Arc
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Arc(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
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

    VECTOR<TV,3> arc=Wedge(joint,pd.pts(j0).pt,pd.pts(j1).pt,pd.half_width,pd.unit_length);
    ARRAY<CONNECTION_DATA> f0,f1;
    ARRAY<int> side0,side1;
    int base=particles.Add_Elements(5);
    particles.X(base)=arc(0)+pipe_dir[0]*pd.unit_length*0.5;
    particles.X(base+1)=arc(0);
    particles.X(base+2)=arc(0)+pipe_dir[1]*pd.unit_length*0.5;
    particles.X(base+3)=arc(0)+arc(1)*2*pd.half_width*pd.unit_length+pipe_dir[0]*pd.unit_length*0.5;
    particles.X(base+4)=arc(0)+arc(2)*2*pd.half_width*pd.unit_length+pipe_dir[1]*pd.unit_length*0.5;
    side0.Append(base);
    side0.Append(base+1);
    side0.Append(base+2);
    side1.Append(base+3);
    int side1_end=base+4;
    f0.Append({base,false});
    f1.Append({base+2,false});

    T unit_angle=2*asin((T)1/(pd.half_width*4));
    T total=acos(arc(1).Dot(arc(2)));
    for(int j=0;(j+1)*unit_angle<=total;j++){
        base=particles.Add_Element();
        T a=-j*unit_angle;
        TV v(cos(a)*arc(1)(0)-sin(a)*arc(1)(1),sin(a)*arc(1)(0)+cos(a)*arc(1)(1));
        particles.X(base)=arc(0)+v*2*pd.half_width*pd.unit_length;
        side1.Append(base);}
    base=particles.Add_Element();
    particles.X(base)=arc(0)+arc(2)*2*pd.half_width*pd.unit_length;
    side1.Append(base);
    side1.Append(side1_end);

    ARRAY<int> prev=side0,verts;
    for(int j=1;j<2*pd.half_width;j++){
        T s=(T)j/(2*pd.half_width);
        verts=Sample_Interpolated(s,side0,side1,pd.unit_length);
        Merge_Interpolated(verts,prev);
        blocks.Append({false});
        prev=verts;
        f0.Append({verts(0),false});
        f1.Append({verts.Last(),false});}
    Merge_Interpolated(side1,prev);
    blocks.Append({false});
    f0.Append({side1(0),false});
    f1.Append({side1.Last(),false});
    if(pd.pipes(p0).x==i) f0.Reverse();
    if(pd.pipes(p1).x!=i) f1.Reverse();
    con.Set({i,p0},f0);
    con.Set({i,p1},f1);
}
//#####################################################################
// Function Generate_Corner
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Corner(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
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

    VECTOR<TV,3> arc=Wedge(joint,pd.pts(j0).pt,pd.pts(j1).pt,pd.half_width,pd.unit_length);
    ARRAY<CONNECTION_DATA> f0,f1;
    ARRAY<int> side0,side1;
    int base=particles.Add_Elements(5);
    particles.X(base)=arc(0)+pipe_dir[0]*pd.unit_length*0.5;
    particles.X(base+1)=arc(0);
    particles.X(base+2)=arc(0)+pipe_dir[1]*pd.unit_length*0.5;
    particles.X(base+3)=arc(0)+arc(1)*2*pd.half_width*pd.unit_length+pipe_dir[0]*pd.unit_length*0.5;
    particles.X(base+4)=arc(0)+arc(2)*2*pd.half_width*pd.unit_length+pipe_dir[1]*pd.unit_length*0.5;
    side0.Append(base);
    side0.Append(base+1);
    side0.Append(base+2);
    side1.Append(base+3);
    int side1_end=base+4;
    f0.Append({base,false});
    f1.Append({base+2,false});

    TV elbow=arc(0)+(joint-arc(0))*2;
    TV start=arc(0)+arc(1)*2*pd.half_width*pd.unit_length;
    TV v=elbow-start;
    T total=v.Normalize();
    for(int j=0;(j+1)*pd.unit_length<=total;j++){
        base=particles.Add_Element();
        particles.X(base)=start+v*j*pd.unit_length;
        side1.Append(base);}
    TV end=arc(0)+arc(2)*2*pd.half_width*pd.unit_length;
    v=end-elbow;
    total=v.Normalize();
    for(int j=0;(j+1)*pd.unit_length<=total;j++){
        base=particles.Add_Element();
        particles.X(base)=elbow+v*j*pd.unit_length;
        side1.Append(base);}
    base=particles.Add_Element();
    particles.X(base)=end;
    side1.Append(base);
    side1.Append(side1_end);

    ARRAY<int> prev=side0,verts;
    for(int j=1;j<2*pd.half_width;j++){
        T s=(T)j/(2*pd.half_width);
        verts=Sample_Interpolated(s,side0,side1,pd.unit_length);
        Merge_Interpolated(verts,prev);
        blocks.Append({false});
        prev=verts;
        f0.Append({verts(0),false});
        f1.Append({verts.Last(),false});}
    Merge_Interpolated(side1,prev);
    blocks.Append({false});
    f0.Append({side1(0),false});
    f1.Append({side1.Last(),false});
    if(pd.pipes(p0).x==i) f0.Reverse();
    if(pd.pipes(p1).x!=i) f1.Reverse();
    con.Set({i,p0},f0);
    con.Set({i,p1},f1);
}
//#####################################################################
// Function March_Corner
//#####################################################################
template<class TV> ARRAY<int> FLUID_LAYOUT_FEM<TV>::
March_Corner(const TV& start_point,int p1,const ARRAY<int>& side,T unit_length)
{
    typedef VECTOR<int,3> E;
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    TV d=particles.X(p1)-start_point;
    T l=d.Normalize();
    int height=(int)(l/unit_length)+1;
    T min_percent=0.3;
    T min_edge_ratio=0.3;
    T edge=(particles.X(side.Last())-particles.X(p1)).Magnitude();
    T rem=l-(height-1)*unit_length;
    bool collapse_last=true;
    if(rem>1e-2*unit_length){
        height++;
        if(rem>min_percent*unit_length && rem/edge>min_edge_ratio)
            collapse_last=false;}
    else collapse_last=false;

    ARRAY<TV> points;
    for(int k=height-2;k>=0;k--){
        if(k==height-2 && collapse_last) continue;
        points.Append(start_point+k*unit_length*d);}
    points.Reverse();

    int base=particles.Add_Elements(points.m);
    ARRAY<int> new_side;
    for(int i=0;i<points.m;i++){
        int pid=base+i;
        particles.X(pid)=points(i);
        new_side.Append(pid);}
    new_side.Append(p1);

    auto next=[height,p1,base,collapse_last](int i)
    {
        int last=collapse_last?height-3:height-2;
        if(i==last) return p1;else return base+i+1;
    };
    auto order=[d,&particles,side](const TV& a,int si)
    {
        T l=a.Dot(d)-(particles.X(side(si))-particles.X(side(0))).Dot(d);
        if(abs(l)<1e-6) return (T)0;
        else return l;
    };
    int p=new_side(0),j=0,i=0,alt=0;
    while(p!=p1){
        int pnext=next(i);
        TV w=particles.X(pnext)-start_point;
        if(j==side.m-1 || order(w,j+1)<0 || (order(w,j+1)==0 && alt==0)){
            area.mesh.elements.Append(E(pnext,p,side(j)));
            p=pnext;
            i++;
            alt=1;}
        else{
            area.mesh.elements.Append(E(side(j+1),p,side(j)));
            j++;
            alt=0;}
        elem_data.Append({blocks.m});}
    if(j==side.m-2){
        area.mesh.elements.Append(E(side(j+1),p,side(j)));
        elem_data.Append({blocks.m});}
    return new_side;
}
//#####################################################################
// Function Generate_Triangle_Junction
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Triangle_Junction(int i,const VECTOR<int,3>& ends,const VECTOR<int,3>& pipes,
    const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    typedef VECTOR<int,3> E;
    GEOMETRY_PARTICLES<TV>& particles=area.particles;
    TV c=pd.pts(i).pt;
    VECTOR<TV,3> tri;
    for(int j=0;j<3;j++){
        int j1=(j+1)%3;
        tri(j)=Wedge(c,pd.pts(ends(j)).pt,pd.pts(ends(j1)).pt,pd.half_width,pd.unit_length)(0);}
    VECTOR<int,3> tri_ids;
    int base=particles.Add_Elements(3);
    for(int j=0;j<3;j++){
        particles.X(base+j)=tri(j);
        tri_ids(j)=base+j;}

    ARRAY<int> tri_f[3];
    for(int j=0;j<3;j++){
        int j0=(j+2)%3;
        TV pipe_dir=(pd.pts(ends(j)).pt-c).Normalized();
        int start=tri_ids(j);
        TV edge_dir=tri(j0)-tri(j);
        TV proj_dir=edge_dir-edge_dir.Dot(pipe_dir)*pipe_dir;
        if(proj_dir.Cross(edge_dir)(0)>0){
            start=tri_ids(j0);
            proj_dir=-proj_dir;
            edge_dir=-edge_dir;}
        int last=start==tri_ids(j0)?tri_ids(j):tri_ids(j0);

        ARRAY<CONNECTION_DATA> f;
        f.Append({start,false});
        tri_f[j].Append(start);
        ARRAY<int> side;
        side.Append(start);
        edge_dir/=(2*pd.half_width);
        proj_dir.Normalize();
        for(int k=0;k<2*pd.half_width;k++){
            TV start_point=particles.X(start)+proj_dir*pd.unit_length*(k+1);
            int s;
            if(k!=2*pd.half_width-1){
                s=particles.Add_Element();
                particles.X(s)=particles.X(start)+edge_dir*(k+1);}
            else s=last;
            tri_f[j].Append(s);
            int elem_num=area.mesh.elements.m;
            side=March_Corner(start_point,s,side,pd.unit_length);
            f.Append({side(0),particles.X(side(0))==start_point});
            if(area.mesh.elements.m>elem_num)
                blocks.Append({false});}
        if((pd.pipes(pipes(j)).x==i && start==tri_ids(j)) ||
            (pd.pipes(pipes(j)).x!=i && start==tri_ids(j0)))
            f.Reverse();
        if(start==tri_ids(j)) tri_f[j].Reverse();
        con.Set({i,pipes(j)},f);}

    ARRAY<int> bottom=tri_f[2];
    TV bottom_dir=particles.X(bottom(1))-particles.X(bottom(0));
    for(int k=2*pd.half_width;k>0;k--){
        int left_index=tri_f[1](k-1);
        TV left=particles.X(left_index);
        ARRAY<int> top;
        top.Append(left_index);
        int base=particles.number;
        if(k>2)
            base=particles.Add_Elements(k-2);
        for(int l=0;l<k-2;++l){
            particles.X(base+l)=left+bottom_dir*(l+1);
            top.Append(base+l);}
        int right_index=tri_f[0](2*pd.half_width-k+1);
        if(right_index!=left_index)
            top.Append(right_index);
        for(int l=1;l<bottom.m;l++){
            area.mesh.elements.Append(E(bottom(l-1),bottom(l),top(l-1)));
            elem_data.Append({blocks.m});}
        for(int l=1;l<top.m;l++){
            area.mesh.elements.Append(E(top(l),top(l-1),bottom(l)));
            elem_data.Append({blocks.m});}
        bottom=top;
        blocks.Append({false});}
}
//#####################################################################
// Function Generate_3_Joint
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_3_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    const ARRAY<int>* pipes=pd.joints.Get_Pointer(i);
    PHYSBAM_ASSERT(pipes && pipes->m==3);
    VECTOR<int,3> ends,ccw_pipes;
    for(int j=0;j<3;j++){
        ccw_pipes(j)=(*pipes)(j);
        ends(j)=pd.pipes((*pipes)(j)).x;
        if(ends(j)==i)
            ends(j)=pd.pipes((*pipes)(j)).y;}
    T signed_area=(pd.pts(ends(1)).pt-pd.pts(ends(0)).pt).Cross(pd.pts(ends(2)).pt-pd.pts(ends(1)).pt)(0);
    if(signed_area<0){
        std::swap(ccw_pipes(1),ccw_pipes(2));
        std::swap(ends(1),ends(2));}

    Generate_Triangle_Junction(i,ends,ccw_pipes,pd,con);
}
//#####################################################################
// Function Generate_Joint
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    const ARRAY<int>* p=pd.joints.Get_Pointer(i);
    PHYSBAM_ASSERT(p);
    if(pd.pts(i).joint_type==default_joint){
        if(p->m==1){
            Generate_End(i,(*p)(0),pd,con);
            return;}
        if(p->m==2){
            Generate_Arc(i,pd,con);
            return;}
        if(p->m==3){
            Generate_3_Joint(i,pd,con);
            return;}}
    if(pd.pts(i).joint_type==corner_joint){
        Generate_Corner(i,pd,con);
        return;}

    PHYSBAM_FATAL_ERROR("joint type not handled");
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Compute(const PARSE_DATA_FEM<TV>& pd)
{
    CONNECTION con;
    for(int i=0;i<pd.pts.m;i++){
        Generate_Joint(i,pd,con);}
    for(int i=0;i<pd.pipes.m;i++){
        Generate_Pipe(i,pd,con);}
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
    VECTOR<T,3> white(1,1,1),dirichlet_bc(1,0,1),traction_bc(0,1,1);
    for(int i=0;i<area.mesh.elements.m;i++){
        //Add_Debug_Object(VECTOR<TV,3>(particles.X.Subset(area.mesh.elements(i))),VECTOR<T,3>(1,1,1));
        for(int j=0;j<3;j++){
            VECTOR<T,3> color=white;
            int v0=area.mesh.elements(i)(j),v1=area.mesh.elements(i)((j+1)%3);
            if(const int* bc0=bc_map.Get_Pointer(v0)){
                const int* bc1=bc_map.Get_Pointer(v1);
                if(bc1 && *bc1==*bc0)
                    color=bc(*bc0).bc_type==dirichlet_v?dirichlet_bc:traction_bc;}
            Add_Debug_Object<TV,2>(VECTOR<TV,2>(particles.X(v0),particles.X(v1)),color);}
    }
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
