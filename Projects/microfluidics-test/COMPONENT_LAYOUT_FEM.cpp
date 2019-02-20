//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <fstream>
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
// Function Compute_Xform
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Set_Cornector(VERTEX_DATA& vd,BLOCK_ID id,int con_id)
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
// Function Compute_Xform
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Emit_Component_Blocks(CANONICAL_COMPONENT* cc,const XFORM& xf,ARRAY<VERTEX_DATA>& vd)
{
    int offset=Value(blocks.m);
    for(BLOCK_ID b(0);b<cc->blocks.m;b++)
    {
        auto& bb=cc->blocks(b);
        for(int c=0;c<bb.connections.m;c++)
        {
            auto& cc=bb.connections(c);
            if(cc.id>=BLOCK_ID()) cc.id+=offset;
            else Set_Cornector(vd(~Value(cc.id)),b,c);
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
        T y=cst.width/n*i-cst.width/2;
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
// Function Make_Canonical_Joint
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Make_Canonical_Joint(const JOINT_KEY& key) -> PAIR<CANONICAL_COMPONENT*,ARRAY<T> >
{
    auto it=canonical_joints.insert({key,{}});
    if(!it.second) return it.first->second;

    // TODO

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
    // BLOCK::master_mask
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
template<class T> CANONICAL_BLOCK_ID COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Merge_Canonical_Blocks(CANONICAL_BLOCK_ID id0,int con_id0,XFORM xf0,
    CANONICAL_BLOCK_ID id1,int con_id1,XFORM xf1)
{
    auto pr=merge_canonical_blocks.Insert(std::make_tuple(id0,con_id0,id1,con_id1),{});
    if(pr.y) return *pr.x;

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

    *pr.x=canonical_blocks.Add_End();
    CANONICAL_BLOCK& cb=canonical_blocks(*pr.x);
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

    return *pr.x;
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

    bl.block=Merge_Canonical_Blocks(bl.block,con_id,bl.xform,bl2.block,con_id2,bl2.xform);
    
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
            irregular_connections(~d).regular;
            irregular_connections(~d).con_id=bl.connections.m;
        }
        else
        {
            blocks(b).connections(d).id=id;
            blocks(b).connections(d).con_id=bl.connections.m;
        }
        bl.connections.Append(bl2.connections(c));
    }

    bl2={CANONICAL_BLOCK_ID(-1)};
}
template class COMPONENT_LAYOUT_FEM<VECTOR<float,2> >;
template class COMPONENT_LAYOUT_FEM<VECTOR<double,2> >;
}
