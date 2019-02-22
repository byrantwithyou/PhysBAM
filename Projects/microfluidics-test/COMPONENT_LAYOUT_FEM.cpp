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
        if(vd.con.edge_on.m)
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
        for(auto& b:a.edge_on)
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
    }

    int h=n/2;
    cb.cross_sections.Append({{h+1,n},{0,h},h,key.type});
    cb.cross_sections.Append({{n,h+n},{h+n+1,2*n},h+n,key.type});
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

    T da=sign(arc_angle)*2*asin(target_length/2*len_arm),offset=0;
    TV v=(arm0-c).Normalized();
    while(abs(arc_angle)>abs(da)*1.5)
    {
        outer.Append(c+ROTATION<TV>::From_Angle(offset).Rotate(v)*len_arm);
        arc_angle-=da;
        offset+=da;
    }

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
    -> ARRAY<PAIR<ARRAY<TV>,ARRAY<IV3> > >
{
    ARRAY<PAIR<ARRAY<TV>,ARRAY<IV3> > > ret(nseg);
    ret(0).x=inner;
    for(int j=1;j<nseg;j++)
    {
        T s=(T)j/nseg;
        ARRAY<TV> inter=Interpolated(s,inner,outer);
        int n0=ret(j-1).x.m,n1=inter.m;
        ret(j-1).x.Append_Elements(inter);
        ret(j).x.Append_Elements(inter);
        ret(j-1).y=Merge_Interpolated(ret(j-1).x,n0,n1);
    }
    int n0=ret(nseg-1).x.m;
    ret(nseg-1).x.Append_Elements(outer);
    ret(nseg-1).y=Merge_Interpolated(ret(nseg-1).x,n0,outer.m);
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
    Fill(cst.num_dofs-1,sides.x,sides.y);
    CANONICAL_COMPONENT* cc=new CANONICAL_COMPONENT;
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
void Cross_Section_Topology(ARRAY<VECTOR<int,3> >& E,F func,int n0,int n1)
{
    int a=0,b=n0,c=n0+n1,ma=a+n0/2,mb=b+n1/2;
    while(a<n0 || b<c)
    {
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

    Cross_Section_Topology(cb.E,[&cb](int a,int b){return cb.X(a+1).y<=cb.X(b+1).y;},n0,n1);

    int h0=n0/2;
    int h1=n1/2;
    cb.cross_sections.Append({{h0+1,n0},{0,h0},h0,key.type[0]});
    cb.cross_sections.Append({{n0,h1+n0},{h1+n0+1,2*n0},h1+n0,key.type[1]});
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
        for(int j=cs.used.min_corner;j<cs.used.max_corner;j++)
            owner(j)=i;
        if(!(mask&(1<<i)) && cb.cross_sections(i).master_index>=0)
            owner(cb.cross_sections(i).master_index)=i;
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
    
    ARRAY<int> index_map(cb1.X.m,use_init,-1);
    const CROSS_SECTION& cs0=cb0.cross_sections(con_id0);
    const CROSS_SECTION& cs1=cb1.cross_sections(con_id1);
    for(int i=cs1.used.min_corner;i<cs1.used.max_corner;i++)
        index_map(i)=i-cs1.used.min_corner+cs0.owned.max_corner;
    for(int i=cs1.owned.min_corner;i<cs1.owned.max_corner;i++)
        index_map(i)=i-cs1.owned.min_corner+cs0.used.max_corner;
    if(cs1.master_index>=0)
        index_map(cs1.master_index)=cs0.master_index;

    int k=cb0.X.m;
    for(int i=0;i<index_map.m;i++)
        if(index_map(i)<0)
            index_map(i)=k++;

    pr.x->x=canonical_blocks.Add_End();
    CANONICAL_BLOCK& cb=canonical_blocks(pr.x->x);
    for(int i=0;i<cb0.cross_sections.m;i++)
        if(i!=con_id0)
            cb.cross_sections.Append(cb0.cross_sections(i));
    for(int i=0;i<cb1.cross_sections.m;i++)
        if(i!=con_id1)
        {
            CROSS_SECTION cs=cb1.cross_sections(i);
            cs.owned.min_corner=index_map(cs.owned.min_corner);
            cs.owned.max_corner=index_map(cs.owned.max_corner);
            cs.used.min_corner=index_map(cs.used.min_corner);
            cs.used.max_corner=index_map(cs.used.max_corner);
            if(cs.master_index>=0) cs.master_index=index_map(cs.master_index);
            cb.cross_sections.Append(cs);
        }

    const MATRIX<T,TV::m>& M0=xforms(xf0.id);
    const MATRIX<T,TV::m>& M1=xforms(xf1.id);
    MATRIX<T,TV::m> M=M0*M1.Inverse();
    TV B=xf0.b-M*xf1.b;

    cb.X=cb0.X;
    cb.X.Resize(k);
    for(int i=0;i<cb1.X.m;i++)
    {
        TV X=M*cb1.X(i)+B;
        int j=index_map(i);
        if(j>=cb0.X.m) cb.X(j)=X;
        else assert((cb.X(j)-X).Magnitude()<(T)1e-6);
    }

    cb.E=cb0.E;
    ARRAY<IV3> E=cb1.E;
    E.Flattened()=index_map.Subset(E.Flattened());
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
        for(auto& j:ic.edge_on)
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
        num-=cs.used.Size();
    return num;
}
//#####################################################################
// Function Fill_Canonical_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Canonical_Block_Matrix(MATRIX_MXN<T>& mat,const CANONICAL_BLOCK& cb)
{
    // TODO
}
//#####################################################################
// Function Fill_Canonical_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Block_Matrix(MATRIX_MXN<T>& mat,BLOCK_ID b)
{
    // TODO
}
//#####################################################################
// Function Fill_Canonical_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Connection_Matrix(MATRIX_MXN<T>& mat,BLOCK_ID b0,int con_id0,BLOCK_ID b1,int con_id1)
{
    // TODO
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

    HASHTABLE<int,int> lookup_matrix_by_hash;

    // Diagonal blocks
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        auto& bl=blocks(b);
        int h=Compute_Block_Hash(b);
        auto pr=lookup_matrix_by_hash.Insert(h,{});
        if(pr.y)
        {
            *pr.x=cem.Create_Matrix_Block(true);
            Fill_Block_Matrix(cem.block_list(*pr.x).M,b);
            // TODO: transforms
        }
        cem.Add_Block_Matrix_Entry(Value(b),Value(b),*pr.x);

        // regular connections
        for(int c=0;c<bl.connections.m;c++)
        {
            if(bl.connections(c).con_id>=0)
            {
                if(bl.connections(c).master)
                {
                    BLOCK_ID b2=bl.connections(c).id;
                    int c2=bl.connections(c).con_id;
                    int h=Compute_Connection_Hash(b,c,b2,c2);
                    auto pr=lookup_matrix_by_hash.Insert(h,{});
                    if(pr.y)
                    {
                        *pr.x=cem.Create_Matrix_Block(false);
                        Fill_Connection_Matrix(cem.block_list(*pr.x).M,b,c,b2,c2);
                        // TODO: transforms
                    }
                    cem.Add_Block_Matrix_Entry(Value(b),Value(b2),*pr.x);
                }
            }
        }
    }

    // regular block, con_id, irregular block, {regular-dof,irregular-dof,regular-owns-id}
    HASHTABLE<std::tuple<CANONICAL_BLOCK_ID,int,CANONICAL_BLOCK_ID,ARRAY<TRIPLE<int,int,bool> > >,int> irregular_matrix_blocks;

    for(auto& ic:irregular_connections)
    {
        HASHTABLE<BLOCK_ID,ARRAY<TRIPLE<int,int,bool> > > h;
        CANONICAL_BLOCK_ID id=blocks(ic.regular).block;
        CROSS_SECTION& cs=canonical_blocks(id).cross_sections(ic.con_id);
        int o=cs.owned.Size();
        int u=cs.used.Size();
        for(int i=0;i<o;i++)
            h.Get_Or_Insert(ic.edge_on(i).x).Append({cs.owned.min_corner+i,ic.edge_on(i).y,true});
        if(cs.master_index>=0)
            h.Get_Or_Insert(ic.edge_on(o).x).Append({cs.master_index,ic.edge_on(o).y,true});
        for(int i=0;i<u;i++)
            h.Get_Or_Insert(ic.edge_on(o+1+i).x).Append({cs.used.min_corner+i,ic.edge_on(o+1+i).y,false});
        for(auto& i:h)
        {
            i.data.Sort();
            auto pr=irregular_matrix_blocks.Insert(std::make_tuple(id,ic.con_id,blocks(i.key).block,i.data),{});
            if(pr.y)
            {
                *pr.x=cem.Create_Matrix_Block(false);
                Fill_Irregular_Connection_Matrix(cem.block_list(*pr.x).M,ic.regular,ic.con_id,i.key,i.data);
                // TODO
            }
            // TODO: transforms
            cem.Add_Block_Matrix_Entry(Value(ic.regular),Value(i.key),*pr.x);
        }
    }

    // TODO: rhs

    cem.End_Fill_Blocks();
}
//#####################################################################
// Function Compute_Block_Hash
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Block_Hash(BLOCK_ID b)
{
    auto& bl=blocks(b);
    int h=Hash_Reduce(bl.block);
    for(int cc=0;cc<bl.connections.m;cc++)
    {
        auto& c=bl.connections(cc);
        if(c.con_id>=0)
        {
            int t=Hash_Reduce(blocks(c.id).block);
            h=int_hash(h,int_hash(t,Hash_Reduce(c.con_id),Hash_Reduce(c.master)));
        }
        else
        {
            auto& ic=irregular_connections(~c.con_id);
            CROSS_SECTION& cs=canonical_blocks(blocks(c.id).block).cross_sections(cc);
            ARRAY<PAIR<int,CANONICAL_BLOCK_ID> > a;
            // Only visit owned and, if necessary, master.
            for(int i=0,n=cs.used.Size()+c.master;i<n;i++)
                a.Append({ic.edge_on(i).y,blocks(ic.edge_on(i).x).block});
            a.Sort();
            h=int_hash(h,Hash(a));
        }
    }
    if(bl.edge_on.m)
    {
        ARRAY<TRIPLE<int,CANONICAL_BLOCK_ID,int> > edge_ons;
        for(int i:bl.edge_on)
        {
            auto& ic=irregular_connections(i);
            for(auto p:ic.edge_on)
                if(p.x==b)
                    edge_ons.Append({p.y,blocks(ic.regular).block,ic.con_id});
        }
        edge_ons.Sort();
        h=int_hash(h,Hash(edge_ons));
    }
    return h;
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
// Function Fill_Irregular_Connection_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Irregular_Connection_Matrix(MATRIX_MXN<T>& mat,BLOCK_ID b0,int con_id0,BLOCK_ID b1,
    const ARRAY<TRIPLE<int,int,bool> >& t)
{
    // TODO:
}
template class COMPONENT_LAYOUT_FEM<VECTOR<double,2> >;
}
