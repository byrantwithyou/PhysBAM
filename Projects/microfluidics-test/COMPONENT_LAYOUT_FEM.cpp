//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/ROTATION.h>
#include <fstream>
#include <list>
#include "CACHED_ELIMINATION_MATRIX.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include <tuple>
namespace PhysBAM{

template<class T>
bool Canonical_Direction(VECTOR<T,2> u)
{
    auto tol=COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::comp_tol;
    if(u.x>tol) return true;
    if(u.x<-tol) return false;
    if(u.x>tol) return true;
    return false;
}

// c <cross-section-name> <num-elements> <width>
// v <vertex-name> <vertex-location-2d>
// j <cross-section-name> <num-pipes> <origin-vertex> [<vertex-name> <connection-name>]*
// p <cross-section-name> <connection-name> <connection-name>
// g <cross-section-name> <cross-section-name> <vertex-name> <vertex-name> <distance> <length> <connection-name> <connection-name>
//#####################################################################
// Function Parse_Input
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Parse_Input(const std::string& pipe_file)
{
    std::ifstream fin(pipe_file);
    std::string line;

    std::map<std::string,VERTEX_ID> pts_index;
    TV pt;
    std::string name,name2,name3,name4,name5;
    int i0;
    T t0,t1;
    char c;
    TV v0;
    HASHTABLE<std::string,CROSS_SECTION_TYPE_ID> cross_section_hash;
    HASHTABLE<std::string,TV> vertices;
    HASHTABLE<std::string,VERTEX_DATA> connection_points;
    xforms.Append(MATRIX<T,2>()+1);
    
    while(getline(fin,line))
    {
        std::stringstream ss(line);
        ss>>c;
        if(isspace(c) || c=='#') continue;
        switch(c)
        {
            case 'c':
                ss>>name>>i0>>t0;
                cross_section_hash.Set(name,Get_Cross_Section_ID({i0,t0}));
                break;

            case 'v':
                ss>>name>>v0;
                vertices.Set(name,{v0});
                break;

            case 'j':
                {
                    ss>>name2>>i0>>name3;
                    JOINT_KEY key;
                    key.type=cross_section_hash.Get(name2);
                    TV O=vertices.Get(name3);
                    ARRAY<TV> dirs;
                    ARRAY<std::string> names;
                    for(int i=0;i<i0;i++)
                    {
                        ss>>name>>name2;
                        TV dir=(vertices.Get(name)-O).Normalized();
                        if(i) key.angles.Append(TV::Oriented_Angle_Between(dirs.Last(),dir));
                        dirs.Append(dir);
                        names.Append(name2);
                    }

                    auto cj=Make_Canonical_Joint(key);
                    XFORM xf={Compute_Xform(dirs(0)),O};
                    ARRAY<VERTEX_DATA> vd(i0);
                    Emit_Component_Blocks(cj.x,xf,vd);
                    for(int i=0;i<i0;i++)
                    {
                        vd(i).X=O+dirs(i)*cj.y(i);
                        connection_points.Set(names(i),vd(i));
                    }
                }
                break;

            case 'p':
                {
                    ss>>name>>name2>>name3;
                    ARRAY<VERTEX_DATA> vd(2);
                    vd(0)=connection_points.Get(name2);
                    vd(1)=connection_points.Get(name3);
                    connection_points.Delete(name2);
                    connection_points.Delete(name3);
                    PIPE_KEY key;
                    key.type=cross_section_hash.Get(name);
                    TV dir=vd(1).X-vd(0).X;
                    if(!Canonical_Direction(dir))
                    {
                        dir=-dir;
                        std::swap(vd(0),vd(1));
                    }
                    key.length=dir.Normalize();
                    XFORM xf={Compute_Xform(dir),vd(0).X};
                    auto cc=Make_Canonical_Pipe(key);
                    Emit_Component_Blocks(cc,xf,vd);
                }
                break;

            case 'g':
                {
                    PIPE_CHANGE_KEY key;
                    ss>>name>>name2;
                    key.type[0]=cross_section_hash.Get(name);
                    key.type[1]=cross_section_hash.Get(name2);
                    ss>>name>>name2>>t0>>t1;
                    TV A=vertices.Get(name);
                    TV B=vertices.Get(name2);
                    TV dir=(B-A).Normalized();
                    if(!Canonical_Direction(dir))
                    {
                        dir=-dir;
                        std::swap(A,B);
                    }
                    TV C=A+dir*t0;
                    TV D=C+dir*t1;
                    key.length=t1;
                    auto cc=Make_Canonical_Pipe_Change(key);
                    XFORM xf={Compute_Xform(dir),A};
                    ss>>name>>name2;

                    ARRAY<VERTEX_DATA> vd(i0);
                    Emit_Component_Blocks(cc,xf,vd);
                    vd(0).X=C;
                    vd(1).X=D;
                    connection_points.Set(name,vd(0));
                    connection_points.Set(name2,vd(1));
                }
                break;
            default:
                LOG::printf("PARSE FAIL: %c %s\n",c,ss.str());
        }
    }
}
//#####################################################################
// Function Set_Connector
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Set_Connector(VERTEX_DATA& vd,BLOCK_ID id,int con_id)
{
    if(vd.con.regular>=BLOCK_ID())
    {
        blocks(vd.con.regular).connections.Append({id,con_id});
        blocks(id).connections.Append({vd.con.regular,vd.con.con_id});
    }
    else
    {
        vd.con.regular=id;
        vd.con.con_id=con_id;
        if(vd.con.edge_on_v.m)
            irregular_connections.Append(vd.con);
    }
}
//#####################################################################
// Function Emit_Component_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Emit_Component_Blocks(const CANONICAL_COMPONENT* cc,const XFORM& xf,ARRAY<VERTEX_DATA>& vd)
{
    int offset=Value(blocks.m);
    for(BLOCK_ID b(0);b<cc->blocks.m;b++)
    {
        auto bb=cc->blocks(b);
        for(int c=0;c<bb.connections.m;c++)
        {
            auto& cc=bb.connections(c);
            if(cc.id>=BLOCK_ID()) cc.id+=offset;
            else Set_Connector(vd(~Value(cc.id)),b,c);
        }
        bb.xform=Compose_Xform(xf,bb.xform);
        blocks.Append(bb);
    }

    for(auto a:cc->irregular_connections)
    {
        for(auto& b:a.edge_on_v)
            b.x+=offset;

        if(a.regular>=BLOCK_ID())
        {
            a.regular+=offset;
            irregular_connections.Append(a);
        }
        else vd(~Value(a.regular)).con=a;
    }
}
//#####################################################################
// Function Compute_Xform
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Xform(const TV& dir) -> XFORM_ID
{
    auto it=xforms_lookup.find(dir);
    if(it!=xforms_lookup.end()) return it->second;
    MATRIX<T,TV::m> M(dir,dir.Orthogonal_Vector());
    XFORM_ID id=xforms.Append(M);
    xforms_lookup[dir]=id;
    return id;
}
//#####################################################################
// Function Compose_Xform
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compose_Xform(const XFORM& a,const XFORM& b) -> XFORM
{
    if(a.id==XFORM_ID()) return {b.id,a.b+b.b};
    if(b.id==XFORM_ID()) return {a.id,a.b+xforms(a.id)*b.b};

    auto pr=xform_comp_table.Insert({a.id,b.id},XFORM_ID());
    if(pr.y)
    {
        *pr.x=xforms.Append(xforms(a.id)*xforms(b.id));
        xform_comp_table.Set({a.id,b.id},*pr.x);
    }
    return {*pr.x,a.b+xforms(a.id)*b.b};
}
//#####################################################################
// Function Make_Canonical_Pipe_Block
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Make_Canonical_Pipe_Block(const PIPE_KEY& key) -> CANONICAL_BLOCK_ID
{
    auto it=canonical_pipe_blocks.insert({key,{}});
    if(!it.second) return it.first->second;

    it.first->second=canonical_blocks.Add_End();
    auto& cb=canonical_blocks(it.first->second);
    const auto& cst=cross_section_types(key.type);

    int n=cst.num_dofs;
    PHYSBAM_ASSERT(n%2);
    cb.X.Resize(2*n);
    cb.S.Resize(4*(n-1)+1);
    for(int i=0;i<n;i++)
    {
        T y=cst.width/(n-1)*i-cst.width/2;
        cb.X(i)=TV(0,y);
        cb.X(i+n)=TV(key.length,y);
    }

    for(int i=0;i<n-1;i++)
    {
        cb.E.Append({i,i+n,i+1});
        cb.E.Append({i+n,i+n+1,i+1});
        cb.S(i)={i,i+1};
        cb.S(i+(n-1))={i+n,i+n+1};
        cb.S(i+2*(n-1))={i+n,i+1};
        cb.S(i+3*(n-1))={i,i+n};
    }
    cb.S(4*(n-1))={n-1,2*n-1};

    cb.cross_sections.Append({{0,n},{0,n-1},false});
    cb.cross_sections.Append({{n,2*n},{n-1,2*(n-1)},true});
    return it.first->second;
}
//#####################################################################
// Function Make_Canonical_Pipe
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Make_Canonical_Pipe(const PIPE_KEY& key) -> CANONICAL_COMPONENT*
{
    auto it=canonical_pipes.insert({key,{}});
    if(!it.second) return it.first->second;

    CANONICAL_COMPONENT* cc=new CANONICAL_COMPONENT;
    it.first->second=cc;

    T length=key.length;
    T offset=0;
    if(length>target_length*(T)1.5)
    {
        CANONICAL_BLOCK_ID id=Make_Canonical_Pipe_Block({key.type,target_length});
        while(length>target_length*(T)1.5)
        {
            cc->blocks.Append(
                {
                    id,
                    {XFORM_ID(),TV(offset,0)},
                    {{cc->blocks.m-1,1},{cc->blocks.m+1,0}}
                });
            length-=target_length;
            offset+=target_length;
        }
    }
    CANONICAL_BLOCK_ID id=Make_Canonical_Pipe_Block({key.type,length});
    cc->blocks.Append(
        {
            id,
            {XFORM_ID(),TV(offset,0)},
            {{cc->blocks.m-1,1},{BLOCK_ID(~1),0}}
        });
    return cc;
}
//#####################################################################
// Function Vertex
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Vertex(T angle,T width) const -> std::tuple<TV,T,T>
{
    PHYSBAM_ASSERT(abs(angle)<2*pi);
    if(angle>pi) angle-=2*pi;
    else if(angle<-pi) angle+=2*pi;
    TV m=ROTATION<TV>::From_Angle(angle/2).Rotated_X_Axis();
    T l=1/sin(abs(angle/2))*width/2;
    T w=tan(pi/2-abs(angle/2))*width/2;
    return std::make_tuple(l*m,w,angle);
}
//#####################################################################
// Function Arc
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
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
// Function Merge_Interpolated
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Merge_Interpolated(const ARRAY<TV>& X,int n0,int n1) const -> ARRAY<IV3>
{
    auto angle=[X](int v0,int v1,int v2)
    {
        TV u=(X(v2)-X(v1)).Normalized();
        TV v=(X(v0)-X(v1)).Normalized();
        return TV::Angle_Between(u,v);
    };
    auto max_angle_index=[angle](int* p)
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

    ARRAY<IV3> elems;
    int i=0,j=n0,alt=0;
    while(i<n0-1 || j<n0+n1-1){
        T a0=0;
        if(i+1<n0) a0=angle(j,i+1,i);
        T a1=0;
        if(j+1<n0+n1) a1=angle(j,j+1,i);
        if(j+1>=n0+n1 || (i+1<n0 && abs(a0-a1)<1e-6 && alt==0) || (i+1<n0 && a0>a1))
        {
            int p[]={i+1,i,j};
            int k=max_angle_index(p);
            elems.Append(IV3(p[k],p[(k+1)%3],p[(k+2)%3]));
            i++;
            alt=1;
        }
        else
        {
            int p[]={j+1,i,j};
            int k=max_angle_index(p);
            elems.Append(IV3(p[k],p[(k+1)%3],p[(k+2)%3]));
            j++;
            alt=0;
        }
    }
    return elems;
}
//#####################################################################
// Function Fill
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill(int nseg,const ARRAY<TV>& inner,const ARRAY<TV>& outer) const
    -> ARRAY<std::tuple<ARRAY<TV>,ARRAY<IV3>,int> >
{
    ARRAY<std::tuple<ARRAY<TV>,ARRAY<IV3>,int> > ret(nseg); // X,E,#verts in the first layer
    std::get<0>(ret(0))=inner;
    for(int j=1;j<nseg;j++)
    {
        T s=(T)j/nseg;
        ARRAY<TV> inter=Interpolated(s,inner,outer);
        int n0=std::get<0>(ret(j-1)).m,n1=inter.m;
        std::get<0>(ret(j-1)).Append_Elements(inter);
        std::get<1>(ret(j-1))=Merge_Interpolated(std::get<0>(ret(j-1)),n0,n1);
        std::get<0>(ret(j)).Append_Elements(inter);
    }
    int n0=std::get<0>(ret(nseg-1)).m;
    std::get<0>(ret(nseg-1)).Append_Elements(outer);
    std::get<1>(ret(nseg-1))=Merge_Interpolated(std::get<0>(ret(nseg-1)),n0,outer.m);
    return ret;
}
//#####################################################################
// Function Interpolate
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Interpolated(T s,const ARRAY<TV>& side0,const ARRAY<TV>& side1) const -> ARRAY<TV>
{
    ARRAY<T> l0(side0.m),l1(side1.m);
    l0(0)=0;
    for(int j=1;j<side0.m;j++)
        l0(j)=l0(j-1)+(side0(j)-side0(j-1)).Magnitude();
    l1(0)=0;
    for(int j=1;j<side1.m;j++)
        l1(j)=l1(j-1)+(side1(j)-side1(j-1)).Magnitude();
    auto point=[](int j,T t,const ARRAY<TV>& side)
    {
        if(j>=side.m-1) return side(j);
        TV p=side(j);
        return p+t*(side(j+1)-p);
    };
    auto loc=[point](T u,const ARRAY<TV>& side,const ARRAY<T>& l)
    {
        T dist=u*l.Last();
        auto iter=std::lower_bound(l.begin(),l.end(),dist);
        PHYSBAM_ASSERT(iter!=l.end());
        int j=iter-l.begin();
        if(*iter==dist) return point(j,0,side);
        else
        {
            T cur=*iter,prev=*(iter-1);
            return point(j-1,(dist-prev)/(cur-prev),side);
        }
    };

    std::list<PAIR<T,TV> > verts;
    TV p0=(1-s)*loc(0,side0,l0)+s*loc(0,side1,l1);
    TV p1=(1-s)*loc(1,side0,l0)+s*loc(1,side1,l1);
    auto begin=verts.insert(verts.end(),{0,p0});
    auto end=verts.insert(verts.end(),{1,p1});
    T max_len=(p1-p0).Magnitude();
    while(max_len>1.5*target_length)
    {
        T u=(begin->x+end->x)*0.5;
        TV v=(1-s)*loc(u,side0,l0)+s*loc(u,side1,l1);
        verts.insert(end,{u,v});
        max_len=0;
        for(auto k=verts.begin();k!=verts.end();k++)
        {
            if(k==verts.begin()) continue;
            auto prev=k;
            prev--;
            T d=(k->y-prev->y).Magnitude();
            if(d>max_len)
            {
                begin=prev;
                end=k;
                max_len=d;
            }
        }
    }
    ARRAY<PARTICLE_ID> indices(verts.size());
    ARRAY<TV> ret;
    for(auto iter=verts.begin();iter!=verts.end();iter++)
        ret.Append(iter->y);
    return ret;
}
//#####################################################################
// Function Make_Canonical_Joint_2
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Make_Canonical_Joint_2(const JOINT_KEY& key) -> PAIR<CANONICAL_COMPONENT*,ARRAY<T> >
{
    const CROSS_SECTION_TYPE& cst=cross_section_types(key.type);
    T width=cst.width;
    TV vert;
    T ext,angle;
    std::tie(vert,ext,angle)=Vertex(key.angles(0),width);
    PAIR<ARRAY<TV>,ARRAY<TV> > sides=Arc(vert,angle,width,0,0);
    // TODO fill in edge dofs S
    auto g=Fill(cst.num_dofs-1,sides.x,sides.y);

    CANONICAL_COMPONENT* cc=new CANONICAL_COMPONENT;
    cc->blocks.Resize(BLOCK_ID(g.m));
    for(int i=0;i<g.m;i++)
    {
        CANONICAL_BLOCK_ID id=canonical_blocks.Add_End();
        CANONICAL_BLOCK& cb=canonical_blocks.Last();
        cb.X=std::move(std::get<0>(g(i)));
        cb.E=std::move(std::get<1>(g(i)));
        // TODO create cross sections, and build connections here.
        // int n0=std::get<2>(g(i)),n1=cb.X.m;
        cc->blocks.Append({id,{XFORM_ID(),TV()},{},{}});
    }
    return {cc,{ext,ext}};
}
//#####################################################################
// Function Make_Canonical_Joint
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Make_Canonical_Joint(const JOINT_KEY& key) -> PAIR<CANONICAL_COMPONENT*,ARRAY<T> >
{
    auto it=canonical_joints.insert({key,{}});
    if(!it.second) return it.first->second;

    if(key.angles.m==0)
    {
        CANONICAL_BLOCK_ID id=canonical_blocks.Add_End();
        CANONICAL_COMPONENT* cc=new CANONICAL_COMPONENT;
        it.first->second={cc,{0}};
        cc->blocks.Append(
            {
                id,
                {XFORM_ID(),TV()},
                {{BLOCK_ID(),1}}
            });
    }
    else if(key.angles.m==1)
    {
        it.first->second=Make_Canonical_Joint_2(key);
    }
    else PHYSBAM_FATAL_ERROR("joint type not supported");

    return it.first->second;
}
//#####################################################################
// Function Make_Canonical_Pipe_Change
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Make_Canonical_Pipe_Change(const PIPE_CHANGE_KEY& key) -> CANONICAL_COMPONENT*
{
    auto it=canonical_changes.insert({key,{}});
    if(!it.second) return it.first->second;

    int num_sec=rint(key.length/target_length);
    num_sec+=!num_sec;

    CANONICAL_COMPONENT* cc=new CANONICAL_COMPONENT;
    it.first->second=cc;

    auto c0=cross_section_types(key.type[0]);
    auto c1=cross_section_types(key.type[1]);
    T w0=c0.width,w1=c1.width,dw=w1-w0,wid=key.length/num_sec;
    int d0=c0.num_dofs,d1=c1.num_dofs,dd=d1-d0;
    for(int i=0;i<num_sec;i++)
    {
        T ow=w0+dw*i/num_sec;
        T nw=w0+dw*(i+1)/num_sec;
        int od=d0+dd*i/num_sec;
        int nd=d0+dd*(i+1)/num_sec;
        auto oid=Get_Cross_Section_ID({od,ow});
        auto nid=Get_Cross_Section_ID({nd,nw});
        PIPE_CHANGE_KEY k={{oid,nid},wid};
        CANONICAL_BLOCK_ID id=Make_Canonical_Change_Block(k);
        cc->blocks.Append(
            {
                id,
                {XFORM_ID(),TV(i*wid,0)},
                {{cc->blocks.m-1,1},{cc->blocks.m+1,0}}
            });
    }
    cc->blocks.Last().connections(1).id=BLOCK_ID(~1);
    return cc;
}
template<class F> // func(a,b) means val(a)<val(b)
void Cross_Section_Topology(ARRAY<VECTOR<int,3> >& E,ARRAY<VECTOR<int,2> >& S,
    F func,int n0,int n1)
{
    for(int i=0;i<n0;i++) S.Append({i,i+1});
    for(int i=0;i<n1;i++) S.Append({i+n0,i+n0+1});
    int a=0,b=n0,c=n0+n1,ma=a+n0/2,mb=b+n1/2;
    while(a<n0 || b<c)
    {
        S.Append({a,b});
        // Ensure that the midpoints are topologically connected.
        bool need_a=(b>=c || (a==ma && b<mb));
        bool need_b=(a>=b || (b==mb && a<ma));
        if(need_a || (!need_b && func(a+1,b+1)))
        {
            E.Append({a,b,a+1});
            a++;
        }
        else
        {
            E.Append({a,b,b+1});
            b++;
        }
    }
    S.Append({a,b});
}
//#####################################################################
// Function Make_Canonical_Change_Block
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Make_Canonical_Change_Block(const PIPE_CHANGE_KEY& key) -> CANONICAL_BLOCK_ID
{
    auto it=canonical_change_blocks.insert({key,{}});
    if(!it.second) return it.first->second;

    it.first->second=canonical_blocks.Add_End();
    auto& cb=canonical_blocks(it.first->second);
    const auto& cst0=cross_section_types(key.type[0]);
    const auto& cst1=cross_section_types(key.type[1]);

    int n0=cst0.num_dofs;
    int n1=cst1.num_dofs;
    PHYSBAM_ASSERT(n0%2);
    PHYSBAM_ASSERT(n1%2);
    cb.X.Resize(n0+n1);
    for(int i=0;i<n0;i++)
        cb.X(i)=TV(0,cst0.width/n0*i-cst0.width/2);
    for(int i=0;i<n1;i++)
        cb.X(i+n0)=TV(key.length,cst1.width/n1*i-cst1.width/2);

    Cross_Section_Topology(cb.E,cb.S,
        [&cb](int a,int b){return cb.X(a+1).y<=cb.X(b+1).y;},n0,n1);

    cb.cross_sections.Append({{0,n0},{0,n0-1},false});
    cb.cross_sections.Append({{n0,2*n0},{n0-1,n0+n1-2},true});
    return it.first->second;
}
//#####################################################################
// Function Get_Cross_Section_ID
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Get_Cross_Section_ID(const CROSS_SECTION_TYPE& cs) -> CROSS_SECTION_TYPE_ID
{
    auto it=cross_section_type_lookup.insert({cs,{}});
    if(it.second) it.first->second=cross_section_types.Append(cs);
    return it.first->second;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute()
{
    Update_Masters();
    Merge_Blocks();

    canonical_block_matrices.Resize(canonical_blocks.m);
    for(CANONICAL_BLOCK_ID i(0);i<canonical_blocks.m;i++)
        Fill_Canonical_Block_Matrix(canonical_block_matrices(i),canonical_blocks(i));
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Update_Masters()
{
    for(BLOCK_ID b0(0);b0<blocks.m;b0++)
    {
        for(int c0=0;c0<blocks(b0).connections.m;c0++)
        {
            auto& n0=blocks(b0).connections(c0);
            BLOCK_ID b1=n0.id;
            int c1=n0.con_id;
            if(c1<0) // irregular
            {
                n0.master=true;
                continue;
            }
            if(b0>b1) continue; // smaller index does the logic
            auto& n1=blocks(b1).connections(c1);
            bool m=true;
            if(c0==0)
            {
                if(c1!=0) m=true;
                else
                {
                    if(blocks(b0).connections.m==2) m=true;
                    else if(blocks(b1).connections.m==2) m=false;
                }
            }
            else
            {
                if(!c1) m=false;
                else if(blocks(b0).connections.m==2) m=true;
                else if(blocks(b1).connections.m==2) m=false;
            }
            n0.master=m;
            n1.master=!m;
        }
    }
}
// F(dof,first_owns)
template<class F>
void Visit_Cross_Section_Dofs(INTERVAL<int> iv,bool uf,bool mast,F func)
{
    int c=iv.min_corner,n=iv.Size();
    if(uf)
    {
        int m=(iv.min_corner+iv.max_corner)/2+(n&mast);
        for(int i=0;i<m;i++) func(c+i,true);
        for(int i=m;i<n;i++) func(c+i,false);
    }
    else
    {
        int m=(iv.min_corner+iv.max_corner)/2+(n&(1-mast));
        for(int i=0;i<m;i++) func(c+i,false);
        for(int i=m;i<n;i++) func(c+i,true);
    }
}
// F(dof0,dof1,first_owns)
template<class F>
void Visit_Cross_Section_Dofs(INTERVAL<int> i0,INTERVAL<int> i1,bool uf0,bool uf1,bool m0,F func)
{
    int a=i1.min_corner,b=1,c=i0.min_corner,n=i0.Size();
    if(uf0!=uf1)
    {
        a=i1.max_corner-1;
        b=-1;
    }
    if(uf0)
    {
        int m=(i0.min_corner+i0.max_corner)/2+(n&m0);
        for(int i=0;i<m;i++) func(c+i,a+b*i,true);
        for(int i=m;i<n;i++) func(c+i,a+b*i,false);
    }
    else
    {
        int m=(i0.min_corner+i0.max_corner)/2+(n&(1-m0));
        for(int i=0;i<m;i++) func(c+i,a+b*i,false);
        for(int i=m;i<n;i++) func(c+i,a+b*i,true);
    }
}
template<class IC,class CS,class FV,class FE>
void Visit_Irregular_Cross_Section_Dofs(const IC& ic,const CS& cs,FV func_v,FE func_e)
{
    Visit_Cross_Section_Dofs(cs.v,{0,cs.v.Size()},cs.own_first,true,true,
        [&ic,&func_v](int a,int b,bool o)
        {
            const auto& p=ic.edge_on_v(b);
            func_v(a,p.x,p.y,o);
        });
    Visit_Cross_Section_Dofs(cs.e,{0,cs.e.Size()},cs.own_first,true,true,
        [&ic,&func_e](int a,int b,bool o)
        {
            const auto& p=ic.edge_on_e(b);
            func_e(a,p.x,p.y,o);
        });
}
template<class CS,class FV,class FE>
void Visit_Regular_Cross_Section_Dofs(const CS& cs0,const CS& cs1,bool m0,FV func_v,FE func_e)
{
    Visit_Cross_Section_Dofs(cs0.v,cs1.v,cs0.own_first,cs1.own_first,m0,func_v);
    Visit_Cross_Section_Dofs(cs0.e,cs1.e,cs0.own_first,cs1.own_first,m0,func_e);
}
//#####################################################################
// Function Separates_Dofs
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Separates_Dofs(BLOCK_ID b)
{
    CANONICAL_BLOCK_ID id=blocks(b).block;
    int mask=0;
    for(int c=0;c<blocks(b).connections.m;c++)
        if(blocks(b).connections(c).master)
            mask|=1<<c;
    auto pr=separates_dofs.Insert({id,mask},0);
    if(!pr.y) return *pr.x;
    
    int valid=0;
    const CANONICAL_BLOCK& cb=canonical_blocks(id);
    ARRAY<int> owner(cb.X.m,use_init,-1);
    for(int i=0;i<cb.cross_sections.m;i++)
    {
        const auto& cs=cb.cross_sections(i);
        int master=(mask>>i)&1;
        INTERVAL<int> iv=cs.v;
        int m=(iv.min_corner+iv.max_corner)/2;
        if(cs.own_first) iv.max_corner=m+(iv.Size()&master);
        else iv.min_corner=m+(iv.Size()&(1-master));
        for(int j:iv) owner(j)=i;
    }
    for(auto e:cb.E)
    {
        IV3 f(owner.Subset(e));
        int o=-1;
        for(auto i:f)
        {
            if(i==-1 || i==o) continue;
            if(o==-1)
            {
                o=i;
                continue;
            }
            valid|=1<<o;
            valid|=1<<i;
        }
    }
    *pr.x=valid;
    return valid;
}
int Merge_Dofs(ARRAY<int>& index_map,int num_dofs0,int num_dofs1,
    INTERVAL<int> i0,INTERVAL<int> i1,bool uf0,bool uf1)
{
    index_map.Resize(num_dofs1,init_all,-1);

    Visit_Cross_Section_Dofs(i0,i1,uf0,uf1,true,[&index_map](int a,int b,bool){index_map(b)=a;});

    int k=num_dofs0;
    for(int i=0;i<index_map.m;i++)
        if(index_map(i)<0)
            index_map(i)=k++;
    return k;
}
PAIR<INTERVAL<int>,bool> Map_Interval(INTERVAL<int> iv,const ARRAY<int>& index_map)
{
    int a=index_map(iv.min_corner);
    int b=index_map(iv.max_corner-1);
    if(a<=b) return {{a,b+1},false};
    return {{b,a+1},true};
}
template<class CS>
CS Map_Cross_Section(CS cs,const ARRAY<int>& index_v_map,const ARRAY<int>& index_e_map)
{
    auto pv=Map_Interval(cs.v,index_v_map);
    auto pe=Map_Interval(cs.e,index_v_map);
    assert(pv.y==pe.y);
    cs.v=pv.x;
    cs.e=pe.x;
    if(pv.y) cs.own_first=!cs.own_first;
    return cs;
}
//#####################################################################
// Function Merge_Blocks
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Merge_Canonical_Blocks(CANONICAL_BLOCK_ID id0,int con_id0,XFORM xf0,
    CANONICAL_BLOCK_ID id1,int con_id1,XFORM xf1) -> PAIR<CANONICAL_BLOCK_ID,ARRAY<int> >*
{
    auto pr=merge_canonical_blocks.Insert(std::make_tuple(id0,con_id0,id1,con_id1),{});
    if(pr.y) return pr.x;

    const CANONICAL_BLOCK& cb0=canonical_blocks(id0);
    const CANONICAL_BLOCK& cb1=canonical_blocks(id1);
    const CROSS_SECTION& cs0=cb0.cross_sections(con_id0);
    const CROSS_SECTION& cs1=cb1.cross_sections(con_id1);
    
    ARRAY<int> index_v_map,index_e_map;
    int num_v_dofs=Merge_Dofs(index_v_map,cb0.X.m,cb1.X.m,cs0.v,cs1.v,
        cs0.own_first,cs1.own_first);
    int num_e_dofs=Merge_Dofs(index_e_map,cb0.X.m,cb1.X.m,cs0.e,cs1.e,
        cs0.own_first,cs1.own_first);

    pr.x->x=canonical_blocks.Add_End();
    CANONICAL_BLOCK& cb=canonical_blocks(pr.x->x);
    for(int i=0;i<cb0.cross_sections.m;i++)
        if(i!=con_id0)
            cb.cross_sections.Append(cb0.cross_sections(i));
    for(int i=0;i<cb1.cross_sections.m;i++)
        if(i!=con_id1)
            cb.cross_sections.Append(
                Map_Cross_Section(cb1.cross_sections(i),index_v_map,index_e_map));

    const MATRIX<T,TV::m>& M0=xforms(xf0.id);
    const MATRIX<T,TV::m>& M1=xforms(xf1.id);
    MATRIX<T,TV::m> M=M0*M1.Inverse();
    TV B=xf0.b-M*xf1.b;

    cb.X=cb0.X;
    cb.X.Resize(num_v_dofs);
    for(int i=0;i<cb1.X.m;i++)
    {
        TV X=M*cb1.X(i)+B;
        int j=index_v_map(i);
        if(j>=cb0.X.m) cb.X(j)=X;
        else assert((cb.X(j)-X).Magnitude()<(T)1e-6);
    }

    cb.S=cb0.S;
    cb.S.Resize(num_e_dofs);
    for(int i=0;i<cb1.E.m;i++)
    {
        int j=index_e_map(i);
        IV s(index_v_map.Subset(cb1.E(i)));
        s.Sort();
        if(j>=cb0.S.m) cb.S(j)=s;
        else assert(cb.S(j)==s);
    }

    cb.E=cb0.E;
    ARRAY<IV3> E=cb1.E;
    E.Flattened()=index_v_map.Subset(E.Flattened());
    cb.E.Append_Elements(E);

    return pr.x;
}
//#####################################################################
// Function Merge_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Merge_Blocks(BLOCK_ID id,int con_id)
{
    PHYSBAM_ASSERT(con_id>=0);
    PHYSBAM_ASSERT(id>=BLOCK_ID());
    BLOCK& bl=blocks(id);
    BLOCK_ID id2=bl.connections(con_id).id;
    int con_id2=bl.connections(con_id).con_id;
    BLOCK& bl2=blocks(id2);

    auto pr=Merge_Canonical_Blocks(bl.block,con_id,bl.xform,bl2.block,con_id2,bl2.xform);
    bl.block=pr->x;
    
    for(int c=con_id;c<bl.connections.m;c++)
    {
        BLOCK_ID b=bl.connections(c).id;
        int d=bl.connections(c).con_id;
        if(d<0) // edge-on
            irregular_connections(~d).con_id=c-1;
        else
            blocks(b).connections(d).con_id=c-1;
    }
    bl.connections.Pop();

    for(int c=0;c<bl2.connections.m;c++)
    {
        if(c==con_id2) continue;
        BLOCK_ID b=bl2.connections(c).id;
        int d=bl2.connections(c).con_id;
        if(d<0) // edge-on
        {
            irregular_connections(~d).regular=id;
            irregular_connections(~d).con_id=bl.connections.m;
        }
        else
        {
            blocks(b).connections(d).id=id;
            blocks(b).connections(d).con_id=bl.connections.m;
        }
        bl.connections.Append(bl2.connections(c));
    }

    for(int i:bl2.edge_on)
    {
        auto& ic=irregular_connections(i);
        for(auto& j:ic.edge_on_v)
            if(j.x==id2)
            {
                j.x=id;
                j.y=pr->y(j.y);
            }
    }

    bl2={CANONICAL_BLOCK_ID(-1)};
}
//#####################################################################
// Function Merge_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Merge_Blocks()
{
    for(BLOCK_ID b(0);b<blocks.m;b++)
        if(int mask=Separates_Dofs(b))
        {
            int besti=-1,best=INT_MAX;
            for(int i=0;i<blocks(b).connections.m;i++)
            {
                if(!(mask&(1<<i))) continue;
                int c=Approx_Dof_Count(blocks(b).connections(i).id);
                if(c<best) // keep blocks small
                {
                    best=c;
                    besti=i;
                }
            }
            Merge_Blocks(b,besti);
            b--; // repeat the check on this block
        }
}
//#####################################################################
// Function Approx_Dof_Count
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Approx_Dof_Count(BLOCK_ID b)
{
    const auto& bl=blocks(b);
    const auto& cb=canonical_blocks(bl.block);
    int num=cb.X.m;
    for(auto& cs:cb.cross_sections)
        num-=cs.v.Size()/2;
    return num;
}
//#####################################################################
// Function Fill_Canonical_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Canonical_Block_Matrix(BLOCK_MATRIX<T>& mat,const CANONICAL_BLOCK& cb)
{
    // TODO
}
//#####################################################################
// Function Copy_Matrix_Data
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Copy_Matrix_Data(BLOCK_MATRIX<T>& A,const BLOCK_MATRIX<T>& B,
    const ARRAY<IV>& va,const ARRAY<IV>& ea,
    const ARRAY<IV>& vb,const ARRAY<IV>& eb) const
{
    // TODO transforms
    for(auto p:va)
    {
        for(auto r:vb)
        {
            A.Add_vv(p.x,r.x,B.Get_vv(p.y,r.y));
            A.Add_vp(p.x,r.x,B.Get_vp(p.y,r.y));
            A.Add_pv(p.x,r.x,B.Get_pv(p.y,r.y));
        }
        for(auto r:eb)
        {
            A.Add_ve(p.x,r.x,B.Get_ve(p.y,r.y));
            A.Add_pe(p.x,r.x,B.Get_pe(p.y,r.y));
        }
    }
    for(auto p:ea)
    {
        for(auto r:vb)
        {
            A.Add_ev(p.x,r.x,B.Get_ev(p.y,r.y));
            A.Add_ep(p.x,r.x,B.Get_ep(p.y,r.y));
        }
        for(auto r:eb)
        {
            A.Add_ee(p.x,r.x,B.Get_ee(p.y,r.y));
        }
    }
}
//#####################################################################
// Function Fill_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Block_Matrix(BLOCK_ID b)
{
    int ib=reference_block(b).y;
    auto& mat=matrix_block_list(ib);
    Init_Block_Matrix(mat,b,b);
    const auto& rd=reference_block_data(reference_block(b).y);
    const auto& bl=blocks(b);
    const auto& cb=canonical_blocks(bl.block);
    mat.Resize();

    // Copy from self
    ARRAY<IV> v,e;
    for(int i=0;i<rd.dof_map_v.m;i++)
        if(rd.dof_map_v(i)>=0)
            v.Append({rd.dof_map_v(i),i});
    for(int i=0;i<rd.dof_map_e.m;i++)
        if(rd.dof_map_e(i)>=0)
            e.Append({rd.dof_map_e(i),i});
    Copy_Matrix_Data(mat,canonical_block_matrices(bl.block),v,e,v,e);

    for(int cc=0;cc<bl.connections.m;cc++)
    {
        auto& c=bl.connections(cc);
        if(c.con_id>=0)
        {
            // Copy from regular neighbors
            ARRAY<IV> v,e;
            const auto& cb2=canonical_blocks(blocks(c.id).block);
            Visit_Regular_Cross_Section_Dofs(cb.cross_sections(cc),
                cb2.cross_sections(c.con_id),c.master,
                [&](int a,int b,bool o)
                {
                    a=rd.dof_map_v(a);
                    assert((a>=0)==o);
                    if(o) v.Append({a,b});
                },
                [&](int a,int b,bool o)
                {
                    a=rd.dof_map_e(a);
                    assert((a>=0)==o);
                    if(o) e.Append({a,b});
                });
            Copy_Matrix_Data(mat,canonical_block_matrices(blocks(c.id).block),v,e,v,e);
        }
        else
        {
            // Copy from irregular neighbors
            HASHTABLE<BLOCK_ID,VECTOR<ARRAY<IV>,2> > prs;
            auto& ic=irregular_connections(~c.con_id);
            Visit_Irregular_Cross_Section_Dofs(ic,
                canonical_blocks(blocks(ic.regular).block).cross_sections(ic.con_id),
                [&](int a,BLOCK_ID b,int c,bool o)
                {
                    a=rd.dof_map_v(a);
                    assert((a>=0)==o);
                    if(o) prs.Get_Or_Insert(b).x.Append({a,c});
                },
                [&](int a,BLOCK_ID b,int c,bool o)
                {
                    a=rd.dof_map_e(a);
                    assert((a>=0)==o);
                    if(o) prs.Get_Or_Insert(b).y.Append({a,c});
                });
            for(auto&p:prs)
                Copy_Matrix_Data(mat,
                    canonical_block_matrices(blocks(p.key).block),
                    p.data.x,p.data.y,p.data.x,p.data.y);
        }
    }
}
//#####################################################################
// Function Fill_Connection_Matrix
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Connection_Matrix(BLOCK_ID b0,int con_id0,BLOCK_ID b1,int con_id1)
{
    auto key=std::make_tuple(blocks(b0).block,con_id0,blocks(b1).block,con_id1);
    auto pr=regular_connection_matrix_blocks.Insert(key,0);
    if(!pr.y) return *pr.x;
    *pr.x=matrix_block_list.Add_End();
    auto& mat=matrix_block_list(*pr.x);

    Init_Block_Matrix(mat,b0,b1);
    PHYSBAM_ASSERT(blocks(b0).connections(con_id0).master);
    mat.Resize();
    const auto& cb0=canonical_blocks(blocks(b0).block);
    const auto& cb1=canonical_blocks(blocks(b1).block);
    const auto& rd0=reference_block_data(reference_block(b0).y);
    const auto& rd1=reference_block_data(reference_block(b1).y);
    ARRAY<IV> va0,ea0,vb0,eb0,va1,ea1,vb1,eb1;
    Visit_Regular_Cross_Section_Dofs(cb0.cross_sections(con_id0),
        cb1.cross_sections(con_id1),true, // 0 is master
        [&](int a,int b,bool o)
        {
            if(o)
            {
                va0.Append({rd0.dof_map_v(a),a});
                va1.Append({rd0.dof_map_v(a),b});
            }
            else
            {
                vb0.Append({rd1.dof_map_v(b),a});
                vb1.Append({rd1.dof_map_v(b),b});
            }
        },
        [&](int a,int b,bool o)
        {
            if(o)
            {
                ea0.Append({rd0.dof_map_e(a),a});
                ea1.Append({rd0.dof_map_e(a),b});
            }
            else
            {
                eb0.Append({rd1.dof_map_e(b),a});
                eb1.Append({rd1.dof_map_e(b),b});
            }
        });
    Copy_Matrix_Data(mat,canonical_block_matrices(blocks(b0).block),va0,ea0,vb0,eb0);
    Copy_Matrix_Data(mat,canonical_block_matrices(blocks(b1).block),va1,ea1,vb1,eb1);
    return *pr.x;
}
//#####################################################################
// Function Fill_Irregular_Connection_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Irregular_Connection_Matrix(IRREGULAR_CONNECTION& ic)
{
    // Copy from irregular neighbors
    HASHTABLE<BLOCK_ID,VECTOR<ARRAY<IV>,8> > prs;
    const auto& rd0=reference_block_data(reference_block(ic.regular).y);
    Visit_Irregular_Cross_Section_Dofs(ic,
        canonical_blocks(blocks(ic.regular).block).cross_sections(ic.con_id),
        [&](int a,BLOCK_ID b,int c,bool o)
        {
            auto& z=prs.Get_Or_Insert(b);
            const auto& rd1=reference_block_data(reference_block(b).y);
            if(o)
            {
                z(0).Append({rd0.dof_map_v(a),a});
                z(1).Append({rd0.dof_map_v(a),c});
            }
            else
            {
                z(2).Append({rd1.dof_map_v(c),a});
                z(3).Append({rd1.dof_map_v(c),c});
            }
        },
        [&](int a,BLOCK_ID b,int c,bool o)
        {
            auto& z=prs.Get_Or_Insert(b);
            const auto& rd1=reference_block_data(reference_block(b).y);
            if(o)
            {
                z(4).Append({rd0.dof_map_e(a),a});
                z(5).Append({rd0.dof_map_e(a),c});
            }
            else
            {
                z(6).Append({rd1.dof_map_e(c),a});
                z(7).Append({rd1.dof_map_e(c),c});
            }
        });

    HASHTABLE<BLOCK_ID,int> h;
    
    for(auto&p:prs)
    {
        auto& z=p.data;
        auto& A=canonical_block_matrices(blocks(ic.regular).block);
        auto& B=canonical_block_matrices(blocks(p.key).block);
        int m_id=matrix_block_list.Add_End();
        h.Set(p.key,m_id);
        auto& M=matrix_block_list(m_id);
        Init_Block_Matrix(M,ic.regular,p.key);
        Copy_Matrix_Data(M,A,z(0),z(4),z(2),z(6));
        Copy_Matrix_Data(M,B,z(1),z(5),z(3),z(7));
    }
    
    int ir=irregular_reference_block_data.Add_End();
    ic.block_data=ir;
    auto& d=irregular_reference_block_data(ir);
    for(auto p:ic.edge_on_v)
    {
        int& z=h.Get(p.x);
        d.add_block.Append(z);
        z=-1;
    }
}
//#####################################################################
// Function Compute_Matrix_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Matrix_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem)
{
    cem.Begin_Fill_Blocks();
    cem.rows.Resize(Value(blocks.m));
    cem.rhs.Resize(Value(blocks.m));

    HASHTABLE<BLOCK_ID,int> lookup_matrix_by_block;

    Compute_Reference_Blocks();
    Compute_Reference_Irregular_Connections();

    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        int i=reference_block(b).y;
        if(reference_block(b).x==b)
        {
            Compute_Dof_Remapping(b);
            Fill_Block_Matrix(b);
        }
        cem.Add_Block_Matrix_Entry(Value(b),Value(b),i);
    }

    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        auto& bl=blocks(b);

        for(int c=0;c<bl.connections.m;c++)
        {
            int con_id=bl.connections(c).con_id;
            if(con_id>=0)
            {
                if(bl.connections(c).master)
                {
                    BLOCK_ID b2=bl.connections(c).id;
                    int i=Fill_Connection_Matrix(b,c,b2,con_id);
                    cem.Add_Block_Matrix_Entry(Value(b),Value(b2),i);
                }
            }
        }
    }

    for(int i=0;i<irregular_connections.m;i++)
    {
        auto& ic=irregular_connections(i);
        int ref=ic.ref_ic;
        if(ref==i) Fill_Irregular_Connection_Matrix(ic);
        int bd=irregular_connections(ref).block_data;
        auto& d=irregular_reference_block_data(bd);
        for(int j=0;j<ic.edge_on_v.m;j++)
            if(d.add_block(j)>=0)
                cem.Add_Block_Matrix_Entry(Value(ic.regular),Value(ic.edge_on_v(j).x),d.add_block(j));
    }

    cem.End_Fill_Blocks();
}
//#####################################################################
// Function Compute_Block_Hash
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Reference_Blocks()
{
    reference_block.Resize(blocks.m,no_init);

    typedef ARRAY<std::tuple<CANONICAL_BLOCK_ID,int,bool> > REG_CON;
    typedef ARRAY<std::tuple<int,CANONICAL_BLOCK_ID,int> > IRREG_CON_V;
    typedef ARRAY<std::tuple<int,CANONICAL_BLOCK_ID,int> > IRREG_CON_E;
    typedef std::tuple<CANONICAL_BLOCK_ID,REG_CON,IRREG_CON_V,IRREG_CON_E> KEY;

    HASHTABLE<KEY,BLOCK_ID> h;

    int index=0;
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        REG_CON reg_con;
        IRREG_CON_V irreg_con_v;
        IRREG_CON_E irreg_con_e;
        
        auto& bl=blocks(b);
        for(int cc=0;cc<bl.connections.m;cc++)
        {
            auto& c=bl.connections(cc);
            if(c.con_id>=0)
                reg_con.Append(std::make_tuple(blocks(c.id).block,c.con_id,c.master));
            else
            {
                reg_con.Append(std::make_tuple(CANONICAL_BLOCK_ID(-1),-1,true));
                
                auto& ic=irregular_connections(~c.con_id);
                Visit_Irregular_Cross_Section_Dofs(ic,
                    canonical_blocks(blocks(ic.regular).block).cross_sections(ic.con_id),
                    [&](int a,BLOCK_ID b,int c,bool o)
                    {
                        if(o) irreg_con_v.Append(std::make_tuple(a,blocks(b).block,c));
                    },
                    [&](int a,BLOCK_ID b,int c,bool o)
                    {
                        if(o) irreg_con_e.Append(std::make_tuple(a,blocks(b).block,c));
                    });
            }
        }
        irreg_con_v.Sort();
        irreg_con_e.Sort();

        auto pr=h.Insert(std::make_tuple(bl.block,reg_con,irreg_con_v,irreg_con_e),{});
        if(pr.y)
            reference_block(b)={b,index++};
        else reference_block(b)=reference_block(*pr.x);
    }

    reference_block_data.Resize(index);
}
//#####################################################################
// Function Compute_Block_Hash
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Reference_Irregular_Connections()
{
    typedef std::tuple<CANONICAL_BLOCK_ID,int,ARRAY<PAIR<CANONICAL_BLOCK_ID,int> >,ARRAY<PAIR<CANONICAL_BLOCK_ID,int> > > KEY;
    HASHTABLE<KEY,int> ref;

    for(int i=0;i<irregular_connections.m;i++)
    {
        auto& ic=irregular_connections(i);
        ARRAY<PAIR<CANONICAL_BLOCK_ID,int> > av,ae;
        for(auto p:ic.edge_on_v) av.Append({blocks(p.x).block,p.y});
        for(auto p:ic.edge_on_e) ae.Append({blocks(p.x).block,p.y});
        auto pr=ref.Insert(std::make_tuple(blocks(ic.regular).block,ic.con_id,av,ae),i);
        ic.ref_ic=*pr.x;
    }
}
//#####################################################################
// Function Compute_Connection_Hash
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Connection_Hash(BLOCK_ID b0,int con_id0,BLOCK_ID b1,int con_id1)
{
    auto id0=blocks(b0).block,id1=blocks(b1).block;
    return Hash(std::make_tuple(id0,con_id0,id1,con_id1));
}
//#####################################################################
// Function Compute_Dof_Remapping
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Dof_Remapping(BLOCK_ID b)
{
    auto& rd=reference_block_data(reference_block(b).y);
    const auto& bl=blocks(b);
    const auto& cb=canonical_blocks(bl.block);
    rd.dof_map_v.Resize(cb.X.m,use_init,1);
    rd.dof_map_e.Resize(cb.S.m,use_init,1);

    for(int cc=0;cc<bl.connections.m;cc++)
    {
        auto& c=bl.connections(cc);
        if(c.con_id>=0)
        {
            const auto& cb2=canonical_blocks(blocks(c.id).block);
            Visit_Regular_Cross_Section_Dofs(cb.cross_sections(cc),
                cb2.cross_sections(c.con_id),c.master,
                [&](int a,int b,bool o)
                {
                    if(!o) rd.dof_map_v(a)=0;
                },
                [&](int a,int b,bool o)
                {
                    if(!o) rd.dof_map_e(a)=0;
                });
        }
        else
        {
            auto& ic=irregular_connections(~c.con_id);
            Visit_Irregular_Cross_Section_Dofs(ic,
                canonical_blocks(blocks(ic.regular).block).cross_sections(ic.con_id),
                [&](int a,BLOCK_ID b,int c,bool o)
                {
                    if(!o) rd.dof_map_v(a)=0;
                },
                [&](int a,BLOCK_ID b,int c,bool o)
                {
                    if(!o) rd.dof_map_e(a)=0;
                });
        }
    }
    int iv=0,ie=0;
    for(int i=0;i<rd.dof_map_v.m;i++)
        rd.dof_map_v(i)=rd.dof_map_v(i)?iv++:-1;
    for(int i=0;i<rd.dof_map_e.m;i++)
        rd.dof_map_e(i)=rd.dof_map_e(i)?ie++:-1;
    rd.num_dofs_v=iv;
    rd.num_dofs_e=ie;
}
//#####################################################################
// Function Init_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Init_Block_Matrix(BLOCK_MATRIX<T>& M,BLOCK_ID a,BLOCK_ID b) const
{
    const auto& c=reference_block_data(reference_block(a).y);
    const auto& d=reference_block_data(reference_block(b).y);
    M.nv_r=c.num_dofs_v;
    M.ne_r=c.num_dofs_e;
    M.nv_c=d.num_dofs_v;
    M.ne_c=d.num_dofs_e;
    M.Resize();
}
template class COMPONENT_LAYOUT_FEM<VECTOR<double,2> >;
}
