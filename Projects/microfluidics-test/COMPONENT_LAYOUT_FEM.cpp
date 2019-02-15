//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <fstream>
#include "COMPONENT_LAYOUT_FEM.h"
namespace PhysBAM{

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
                    ARRAY<VERTEX_DATA> vd;
                    vd(0)=connection_points.Get(name2);
                    vd(1)=connection_points.Get(name3);
                    connection_points.Delete(name2);
                    connection_points.Delete(name3);
                    PIPE_KEY key;
                    key.type=cross_section_hash.Get(name);
                    TV dir=vd(1).X-vd(0).X;
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
                    TV u=(B-A).Normalized();
                    TV C=A+u*t0;
                    TV D=C+u*t1;
                    key.length=t1;
                    auto cc=Make_Canonical_Pipe_Change(key);
                    XFORM xf={Compute_Xform((B-A).Normalized()),A};
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
        connections.Append({{vd.con.regular,id},{vd.con.con_id,con_id}});
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
    for(BLOCK_ID i(0);i<cc->blocks.m;i++)
        blocks.Append({cc->blocks(i).block,Compose_Xform(xf,cc->blocks(i).xform)});

    for(auto a:cc->connections)
    {
        if(a.id[0]>=BLOCK_ID()) a.id[0]+=offset;
        if(a.id[1]>=BLOCK_ID()) a.id[1]+=offset;
        if(a.id[0]<BLOCK_ID()) Set_Cornector(vd(~Value(a.id[0])),a.id[1],a.con_id[1]);
        else if(a.id[1]<BLOCK_ID()) Set_Cornector(vd(~Value(a.id[1])),a.id[0],a.con_id[0]);
        else connections.Append(a);
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

    XFORM_ID id;
    if(XFORM_ID* pid=xform_comp_table.Get_Pointer({a.id,b.id}))
        id=*pid;
    else
    {
        id=xforms.Append(xforms(a.id)*xforms(b.id));
        xform_comp_table.Set({a.id,b.id},id);
    }
    return {id,a.b+xforms(a.id)*b.b};
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
        CANONICAL_BLOCK_ID target_id=Make_Canonical_Pipe_Block({key.type,target_length});
        while(length>target_length*(T)1.5)
        {
            cc->blocks.Append({target_id,{XFORM_ID(),TV(offset,0)},1});
            length-=target_length;
            offset+=target_length;
        }
    }
    CANONICAL_BLOCK_ID extra_id=Make_Canonical_Pipe_Block({key.type,length});
    cc->blocks.Append({extra_id,{XFORM_ID(),TV(offset,0)},1});
    for(BLOCK_ID i(0);i<cc->blocks.m;i++)
        cc->connections.Append({{i-1,i},{1,0}});
    cc->connections(0).id[0]=BLOCK_ID(~0);
    cc->connections.Last().id[1]=BLOCK_ID(~1);
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
        cc->blocks.Append({id,{XFORM_ID(),TV(i*wid,0)},1});
    }
    for(BLOCK_ID i(0);i<cc->blocks.m;i++)
        cc->connections.Append({{i-1,i},{1,0}});
    cc->connections(0).id[0]=BLOCK_ID(~0);
    cc->connections.Last().id[1]=BLOCK_ID(~1);
    return cc;
}
template<class F> // func(a,b) means val(a)<val(b)
void Cross_Section_Topology(ARRAY<VECTOR<int,3> >& E,F func,int n0,int n1)
{
    int a=0,b=n0,c=n0+n1;
    while(a<n0 || b<c)
    {
        if(b>=c || (a<b && func(a+1,b+1)))
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
    // BLOCK::master_mask
    // COMPONENT::master_mask











    

    

















    
}
template class COMPONENT_LAYOUT_FEM<VECTOR<float,2> >;
template class COMPONENT_LAYOUT_FEM<VECTOR<double,2> >;
}
