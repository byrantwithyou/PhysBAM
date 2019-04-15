//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/ROTATION.h>
#include "BLOCK_MESHING_ITERATOR.h"
#include "COMPONENT_JOINT.h"
namespace PhysBAM{

//#####################################################################
// Function Make_Component
//#####################################################################
template<class T> PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > COMPONENT_JOINT<T>::
Make_Component(int d,T width,const ARRAY<T>& angles)
{
    PHYSBAM_ASSERT(angles.m>=1);
    JOINT_KEY<T> key={d,width,angles};
    auto it=canonical_joints.insert({key,{}});
    if(!it.second) return it.first->second;

    if(angles.m==1)
        it.first->second=Make_Joint_2(d,width,angles);
    else if(angles.m==2)
        it.first->second=Make_Joint_3(d,width,angles);
    else PHYSBAM_FATAL_ERROR("joint type not supported");
    return it.first->second;
}
//#####################################################################
// Function Make_Joint_2
//#####################################################################
template<class T> PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > COMPONENT_JOINT<T>::
Make_Joint_2(int d,T width,const ARRAY<T>& angles)
{
    T sep=target_length;
    TV vert;
    T ext,angle;
    std::tie(vert,ext,angle)=Elbow_Pit(angles(0),width);
    PAIR<ARRAY<TV>,ARRAY<TV> > sides=Arc(vert,angle,width,sep,sep);
    if(angle>0) std::swap(sides.x,sides.y);

    CANONICAL_COMPONENT<T>* cc=new CANONICAL_COMPONENT<T>;
    cc->blocks.Resize(CC_BLOCK_ID(d-1));
    CC_IRREG_ID ic0=cc->irregular_connections.Append({CC_BLOCK_ID(~0)});
    CC_IRREG_ID ic1=cc->irregular_connections.Append({CC_BLOCK_ID(~1)});
    for(BLOCK_MESHING_ITERATOR<TV> it(sides.x,sides.y,d-1,target_length);it.Valid();it.Next())
    {
        CANONICAL_BLOCK<T>* cb=new CANONICAL_BLOCK<T>;
        it.Build(cb->X,cb->E,cb->S);
        ARRAY<CC_BLOCK_CONNECTION,CON_ID> con;
        Joint_Connection(0,it,cb,
            &cc->irregular_connections(ic0),
            &cc->irregular_connections(ic1),con,it.k==1?CON_ID(0):CON_ID(1));
        if(it.k==0)
        {
            for(int j=0;j<it.X0.m;j++) cb->bc_v.Append(j);
            for(int j=0;j<it.First_Diagonal_Edge();j++) cb->bc_e.Append(j);
        }
        if(it.k==it.nseg-1)
        {
            for(int j=it.X0.m;j<cb->X.m;j++) cb->bc_v.Append(j);
            for(int j=it.Last_Diagonal_Edge()+1;j<cb->S.m;j++) cb->bc_e.Append(j);
        }
        cc->blocks(CC_BLOCK_ID(it.k))={cb,{},con,{{ic0,d-2-it.k},{ic1,it.k}}};
    }
    cc->irregular_connections(ic0).edge_on.Reverse();
    for(auto& e:cc->irregular_connections(ic0).edge_on)
        std::swap(e.v0,e.v1);
    return {cc,{ext+sep,ext+sep}};
}
//#####################################################################
// Function Make_Joint_3
//#####################################################################
template<class T> PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > COMPONENT_JOINT<T>::
Make_Joint_3(int d,T width,const ARRAY<T>& angles)
{
    T minimum_angle=min(2*(T)pi-angles.Sum(),angles.Min());
    if(minimum_angle<pi/4) return Make_Joint_3_Small(d,width,angles);
    else return Make_Joint_3_Average(d,width,angles);
}
//#####################################################################
// Function Make_Joint_3_Small
//#####################################################################
template<class T> PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > COMPONENT_JOINT<T>::
Make_Joint_3_Small(int d,T width,const ARRAY<T>& angles)
{
    T sep=target_length;
    ARRAY<T> a(angles);
    a.Append(2*pi-angles.Sum());
    int k=0;
    T min_abs_angle=a(k),tot=0;
    VECTOR<TV,3> dirs,w;
    for(int i=0;i<3;i++)
    {
        if(a(i)<min_abs_angle)
        {
            k=i;
            min_abs_angle=a(i);
        }
        dirs(i)=ROTATION<TV>::From_Angle(tot).Rotated_X_Axis();
        w(i)=ROTATION<TV>::From_Angle(tot).Rotate(Elbow_Pit_Oriented(a(i),width));
        tot+=a(i);
    }

    VECTOR<TV,2> e=Extrude(w((k+2)%3),w((k+1)%3),dirs((k+2)%3));
    e+=dirs((k+2)%3)*sep;
    ARRAY<TV> l0={e(0),w((k+2)%3)},l1={e(1),w((k+1)%3)};
    if(a((k+1)%3)>a((k+2)%3)+pi/10)
        l1.Append(w((k+1)%3)+dirs((k+1)%3)*sep);
    if(a((k+2)%3)>a((k+1)%3)+pi/10)
        l0.Append(w((k+2)%3)+dirs(k)*sep);
    l0.Reverse();
    l1.Reverse();
    ARRAY<TV> s0=Polyline(l0,target_length),s1=Polyline(l1,target_length);
    CANONICAL_COMPONENT<T>* cc=new CANONICAL_COMPONENT<T>;
    ARRAY<TV> t1;
    t1.Append(s0(0));
    CC_IRREG_ID index=cc->irregular_connections.Add_End();
    CC_IRREG_ID ick2=cc->irregular_connections.Append({CC_BLOCK_ID(~((k+2)%3))});
    for(BLOCK_MESHING_ITERATOR<TV> it(s0,s1,d-1,target_length);it.Valid();it.Next())
    {
        CANONICAL_BLOCK<T>* cb=new CANONICAL_BLOCK<T>;
        it.Build(cb->X,cb->E,cb->S);
        ARRAY<CC_BLOCK_CONNECTION,CON_ID> con;
        Joint_Connection(0,it,cb,
            &cc->irregular_connections(index),
            &cc->irregular_connections(ick2),con,it.k==1?CON_ID(0):CON_ID(1));
        for(int i=0;it.k==0 && i<it.X0.m;i++) cb->bc_v.Append(i);
        for(int i=0;it.k==0 && i<it.First_Diagonal_Edge();i++) cb->bc_e.Append(i);
        for(int i=it.X0.m;it.k==it.nseg-1 && i<cb->X.m;i++) cb->bc_v.Append(i);
        for(int i=it.Last_Diagonal_Edge()+1;it.k==it.nseg-1 && i<cb->S.m;i++) cb->bc_e.Append(i);
        cc->blocks.Append({cb,{},con,{{index,d-2-it.k},{ick2,it.k}}});
        t1.Append(it.X1(0));
    }
    t1.Reverse();
    VECTOR<TV,2> g0=Extrude(w((k+2)%3),w(k),dirs(k)),g1=Extrude(w((k+1)%3),w(k),dirs((k+1)%3));
    ARRAY<TV> t0=Polyline({g1(0),w(k),g0(0)},width/(d-1));
    int nseg=rint(w(k).Magnitude()/target_length);
    int offset=Value(cc->blocks.m);
    CANONICAL_BLOCK<T>* sep_cb=new CANONICAL_BLOCK<T>;
    for(BLOCK_MESHING_ITERATOR<TV> it(t0,t1,nseg,target_length);it.Valid();it.Next())
    {
        CANONICAL_BLOCK<T>* cb=new CANONICAL_BLOCK<T>;
        it.Build(cb->X,cb->E,cb->S);
        ARRAY<CC_BLOCK_CONNECTION,CON_ID> con;
        if(it.k==0)
        {
            sep_cb->X.Append_Elements(it.X0);
            cb->cross_sections.Append({{0,it.X0.m},{0,it.First_Diagonal_Edge()},false});
            con.Append({CC_BLOCK_ID(offset+nseg),CON_ID(0)});
            cb->bc_v.Append(d-1);
        }
        Joint_Connection(offset,it,cb,nullptr,nullptr,con,CON_ID(1));
        if(it.k==nseg-1)
        {
            cb->cross_sections.Append({{it.X0.m,cb->X.m},{it.Last_Diagonal_Edge()+1,cb->S.m},true});
            con.Append({index});
            cc->irregular_connections(index).regular=CC_BLOCK_ID(offset+it.k);
            cc->irregular_connections(index).con_id=CON_ID(1);
        }
        cb->bc_v={0,it.X0.m-1,it.X0.m,cb->X.m-1};
        cb->bc_e={it.First_Diagonal_Edge(),it.Last_Diagonal_Edge()};
        cc->blocks.Append({cb,{},con,{}});
    }

    PHYSBAM_ASSERT(sep_cb->X.m==2*d-1);
    sep_cb->S.Resize(2*(4*d-3));
    for(int j=0;j<sep_cb->X.m-1;j++) sep_cb->S(j)={j,j+1};
    sep_cb->cross_sections.Append({{0,sep_cb->X.m},{0,sep_cb->X.m-1},true});
    for(int i=0;i<2;i++)
    {
        int offset_x=i*(d-1),s=2*d-1,offset_edge=s-1+i*(3*(d-1)+1);
        for(int j=0;j<d;j++) sep_cb->X.Append(sep_cb->X(offset_x+j)+dirs((k+1-i)%3)*sep);
        for(int j=0;j<d-1;j++)
        {
            int p=offset_x+j,q=p+s+i;
            sep_cb->E.Append({q,p,q+1});
            sep_cb->E.Append({p+1,q+1,p});

            sep_cb->S(offset_edge+j)={q,q+1};
            sep_cb->S(offset_edge+j+(d-1))={q+1,p};
            sep_cb->S(offset_edge+j+2*(d-1))={p,q};
        }
        sep_cb->S(offset_edge+3*(d-1))={offset_x+d-1,offset_x+d-1+i+s};
        sep_cb->bc_v.Append_Elements(ARRAY<int>{offset_x,offset_x+i+s,offset_x+d-1,offset_x+d-1+i+s});
        sep_cb->bc_e.Append_Elements(ARRAY<int>{offset_edge+2*(d-1),offset_edge+3*(d-1)});
        sep_cb->cross_sections.Append({{offset_x+s+i,offset_x+d+s+i},{offset_edge,offset_edge+(d-1)},false});
    }

    ARRAY<CC_BLOCK_CONNECTION,CON_ID> sep_con;
    sep_con.Append({CC_BLOCK_ID(offset),CON_ID(0)});
    sep_con.Append({CC_BLOCK_ID(~((k+1)%3)),CON_ID(1)});
    sep_con.Append({CC_BLOCK_ID(~k),CON_ID(1)});
    cc->blocks.Append({sep_cb,{},sep_con,{}});

    cc->irregular_connections(index).edge_on.Reverse();
    for(auto& e:cc->irregular_connections(index).edge_on)
        std::swap(e.v0,e.v1);
    ARRAY<T> ext={g1.Average().Magnitude()+sep,g0.Average().Magnitude()+sep,e.Average().Magnitude()};
    return {cc,{ext((3-k)%3),ext((4-k)%3),ext((5-k)%3)}};
}
//#####################################################################
// Function Make_Joint_3_Average
//#####################################################################
template<class T> PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > COMPONENT_JOINT<T>::
Make_Joint_3_Average(int d,T width,const ARRAY<T>& angles)
{
    T sep=2*target_length;
    ARRAY<T> a(angles);
    a.Append(2*pi-angles.Sum());
    int k=0;
    T max_abs_angle=a(k),tot=0;
    VECTOR<TV,3> dirs,w;
    for(int i=0;i<3;i++)
    {
        if(a(i)>max_abs_angle)
        {
            k=i;
            max_abs_angle=a(i);
        }
        dirs(i)=ROTATION<TV>::From_Angle(tot).Rotated_X_Axis();
        w(i)=ROTATION<TV>::From_Angle(tot).Rotate(Elbow_Pit_Oriented(a(i),width));
        tot+=a(i);
    }

    CANONICAL_COMPONENT<T>* cc=new CANONICAL_COMPONENT<T>;
    VECTOR<TV,2> e=Extrude(w((k+2)%3),w((k+1)%3),dirs((k+2)%3));
    e+=dirs((k+2)%3)*sep;
    ARRAY<TV> s0=Polyline({w((k+2)%3),e(0)},target_length),s1=Polyline({w((k+1)%3),e(1)},target_length),b;
    CC_IRREG_ID index=cc->irregular_connections.Add_End();
    CC_IRREG_ID ick2=cc->irregular_connections.Append({CC_BLOCK_ID(~((k+2)%3))});
    for(BLOCK_MESHING_ITERATOR<TV> it(s0,s1,d-1,target_length);it.Valid();it.Next())
    {
        CANONICAL_BLOCK<T>* cb=new CANONICAL_BLOCK<T>;
        it.Build(cb->X,cb->E,cb->S);
        ARRAY<CC_BLOCK_CONNECTION,CON_ID> con;
        Joint_Connection(0,it,cb,
            &cc->irregular_connections(index),
            &cc->irregular_connections(ick2),con,it.k==1?CON_ID(0):CON_ID(1));
        for(int i=0;it.k==0 && i<it.X0.m;i++) cb->bc_v.Append(i);
        for(int i=0;it.k==0 && i<it.First_Diagonal_Edge();i++) cb->bc_e.Append(i);
        for(int i=it.X0.m;it.k==it.nseg-1 && i<cb->X.m;i++) cb->bc_v.Append(i);
        for(int i=it.Last_Diagonal_Edge()+1;it.k==it.nseg-1 && i<cb->S.m;i++) cb->bc_e.Append(i);
        cc->blocks.Append({cb,{},con,{{index,d-2-it.k},{ick2,it.k}}});
        b.Append(it.X1(0));
    }
    VECTOR<TV,2> e0=Extrude(w(k),w((k+1)%3),dirs((k+1)%3)),e1=Extrude(w(k),w((k+2)%3),dirs(k));
    e0+=dirs((k+1)%3)*sep;
    e1+=dirs(k)*sep;
    ARRAY<TV> t1=Polyline({e1(0),w(k),e0(0)},target_length);
    ARRAY<TV> t0;
    t0.Append_Elements(Polyline({e1(1),w((k+2)%3)},target_length));
    int n0=t0.m;
    t0.Append_Elements(b);
    t0.Append_Elements(Polyline({t0.Pop_Value(),e0(1)},target_length));
    CC_IRREG_ID ick1=cc->irregular_connections.Append({CC_BLOCK_ID(~((k+1)%3))});
    CC_IRREG_ID ick=cc->irregular_connections.Append({CC_BLOCK_ID(~k)});
    int offset=Value(cc->blocks.m);
    for(BLOCK_MESHING_ITERATOR<TV> it(t0,t1,d-1,target_length);it.Valid();it.Next())
    {
        CANONICAL_BLOCK<T>* cb=new CANONICAL_BLOCK<T>;
        it.Build(cb->X,cb->E,cb->S);
        ARRAY<CC_BLOCK_CONNECTION,CON_ID> con;
        if(it.k==0)
        {
            cb->cross_sections.Append({{n0-1,n0+b.m},{n0-1,n0+b.m-1},false});
            con.Append({index});
            cc->irregular_connections(index).regular=CC_BLOCK_ID(offset+it.k);
            for(int i=0;i<it.X0.m;i++)
            {
                if(i<n0 || i>=n0+b.m-1) cb->bc_v.Append(i);
                if((i>0 && i<n0) || i>(n0+b.m-1)) cb->bc_e.Append(i-1);
            }
        }
        Joint_Connection(offset,it,cb,
            &cc->irregular_connections(ick),
            &cc->irregular_connections(ick1),con,CON_ID(1));
        for(int i=it.X0.m;it.k==it.nseg-1 && i<cb->X.m;i++) cb->bc_v.Append(i);
        for(int i=it.Last_Diagonal_Edge()+1;it.k==it.nseg-1 && i<cb->S.m;i++) cb->bc_e.Append(i);
        cc->blocks.Append({cb,{},con,{{ick,d-2-it.k},{ick1,it.k}}});
    }
    cc->irregular_connections(index).edge_on.Reverse();
    for(auto& e:cc->irregular_connections(index).edge_on)
        std::swap(e.v0,e.v1);
    cc->irregular_connections(ick).edge_on.Reverse();
    for(auto& e:cc->irregular_connections(ick).edge_on)
        std::swap(e.v0,e.v1);
    ARRAY<T> ext={e1.Average().Magnitude(),e0.Average().Magnitude(),e.Average().Magnitude()};
    return {cc,{ext((3-k)%3),ext((4-k)%3),ext((5-k)%3)}};
}
//#####################################################################
// Function Elbow_Pit
//#####################################################################
template<class T> auto COMPONENT_JOINT<T>::
Elbow_Pit(T angle,T width) const -> std::tuple<TV,T,T>
{
    if(angle>pi) angle-=2*pi;
    TV m=ROTATION<TV>::From_Angle(angle/2).Rotated_X_Axis();
    T l=1/sin(abs(angle/2))*width/2;
    T w=std::tan(pi/2-abs(angle/2))*width/2;
    return std::make_tuple(l*m,w,angle);
}
//#####################################################################
// Function Elbow_Pit_Oriented
//#####################################################################
template<class T> auto COMPONENT_JOINT<T>::
Elbow_Pit_Oriented(T angle,T width) const -> TV
{
    TV m=ROTATION<TV>::From_Angle(angle/2).Rotated_X_Axis();
    T s=width/2/sin(angle/2);
    return s*m;
}
//#####################################################################
// Function Extrude
//#####################################################################
template<class T> auto COMPONENT_JOINT<T>::
Extrude(const TV& v0,const TV& v1,const TV& n) const -> VECTOR<TV,2>
{
    VECTOR<TV,2> e(v0,v1);
    T d=(e(1)-e(0)).Dot(n);
    if(d>=0) e(0)+=n*d;
    else e(1)-=n*d;
    return e;
}
//#####################################################################
// Function Arc
//#####################################################################
template<class T> auto COMPONENT_JOINT<T>::
Arc(const TV& c,T angle,T len_arm,T ext0,T ext1) const -> PAIR<ARRAY<TV>,ARRAY<TV> >
{
    PHYSBAM_ASSERT(abs(angle)<=pi);
    T arc_angle=angle>0?angle-pi:angle+pi;
    TV arm0(c(0),-c(1));
    TV arm1=c+ROTATION<TV>::From_Angle(arc_angle).Rotate(arm0-c);
    ARRAY<TV> inner,outer;
    if(ext0)
    {
        inner.Append(c+TV(ext0,0));
        outer.Append(arm0+TV(ext0,0));
    }
    inner.Append(c);

    T da=sign(arc_angle)*2*asin(target_length/2*len_arm),offset=0,a=arc_angle;
    TV v=(arm0-c).Normalized();
    while(abs(a)>abs(da)*1.5)
    {
        outer.Append(c+ROTATION<TV>::From_Angle(offset).Rotate(v)*len_arm);
        a-=da;
        offset+=da;
    }
    outer.Append(c+ROTATION<TV>::From_Angle(arc_angle).Rotate(v)*len_arm);

    if(ext1)
    {
        TV d=ROTATION<TV>::From_Angle(angle).Rotated_X_Axis();
        inner.Append(c+d*ext1);
        outer.Append(arm1+d*ext1);
    }
    return {inner,outer};
}
//#####################################################################
// Function Polyline
//#####################################################################
template<class T> auto COMPONENT_JOINT<T>::
Polyline(const ARRAY<TV>& points,T dx) const -> ARRAY<TV>
{
    ARRAY<TV> verts;
    for(int j=1;j<points.m;j++)
    {
        TV v=points(j)-points(j-1);
        T l=v.Normalize();
        int n=std::round(l/dx);
        T dx=l/n;
        for(int k=0;k<n;k++) verts.Append(points(j-1)+v*k*dx);
    }
    verts.Append(points.Last());
    return verts;
}
//#####################################################################
// Function Joint_Connection
//#####################################################################
template<class T> void COMPONENT_JOINT<T>::
Joint_Connection(int offset,BLOCK_MESHING_ITERATOR<TV>& it,
    CANONICAL_BLOCK<T>* cb,CC_IRREGULAR_CONNECTION* ic0,CC_IRREGULAR_CONNECTION* ic1,
    ARRAY<CC_BLOCK_CONNECTION,CON_ID>& con,CON_ID prev) const
{
    CC_BLOCK_ID self(offset+it.k);
    if(it.k!=0)
    {
        cb->cross_sections.Append({{0,it.X0.m},{0,it.First_Diagonal_Edge()},false});
        con.Append({CC_BLOCK_ID(offset+it.k-1),prev});
    }
    if(it.k!=it.nseg-1)
    {
        cb->cross_sections.Append({{it.X0.m,cb->X.m},{it.Last_Diagonal_Edge()+1,cb->S.m},true});
        con.Append({CC_BLOCK_ID(offset+it.k+1),CON_ID(0)});
    }
    if(ic0) ic0->edge_on.Append({self,it.First_Diagonal_Edge(),0,it.X0.m});
    if(ic1) ic1->edge_on.Append({self,it.Last_Diagonal_Edge(),it.X0.m-1,cb->X.m-1});
}

template class COMPONENT_JOINT<float>;
template class COMPONENT_JOINT<double>;
}
