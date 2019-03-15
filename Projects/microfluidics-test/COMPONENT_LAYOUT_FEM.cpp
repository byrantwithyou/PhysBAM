//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Data_Structures/TUPLE.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <fstream>
#include "CACHED_ELIMINATION_MATRIX.h"
#include "COMPONENT_BC.h"
#include "COMPONENT_CHANGE.h"
#include "COMPONENT_JOINT.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "COMPONENT_PIPE.h"
#include <tuple>
namespace PhysBAM{

double comp_tol=1e-10;

template<class T>
bool Canonical_Direction(VECTOR<T,2> u)
{
    if(u.x>comp_tol) return true;
    if(u.x<-comp_tol) return false;
    if(u.x>comp_tol) return true;
    return false;
}

// l <target-length>
// c <cross-section-name> <num-elements> <width>
// v <vertex-name> <vertex-location-2d>
// j <cross-section-name> <num-pipes> <origin-vertex> [<vertex-name> <connection-name>]*
// p <cross-section-name> <connection-name> <connection-name>
// g <cross-section-name> <cross-section-name> <vertex-name> <vertex-name> <distance> <length> <connection-name> <connection-name>
// u <cross-section-name> <origin-vertex> <vertex-name> <connection-name> <flow-rate>
// t <cross-section-name> <origin-vertex> <vertex-name> <connection-name> <traction-2d>
// T <cross-section-name> <origin-vertex> <vertex-name> <connection-name>
// U <analytic-velocity>
// P <analytic-pressure>
//#####################################################################
// Function Destructor
//#####################################################################
template<class T> COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
~COMPONENT_LAYOUT_FEM()
{
    delete analytic_velocity;
    delete analytic_pressure;
    for(auto& p:canonical_block_matrices) delete p.key;
}
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
    HASHTABLE<std::string,PAIR<int,T> > cross_section_hash;
    HASHTABLE<std::string,TV> vertices;
    HASHTABLE<std::string,VERTEX_DATA> connection_points;
    COMPONENT_PIPE<T> comp_pipe;
    COMPONENT_CHANGE<T> comp_change;
    COMPONENT_BC<T> comp_bc;
    COMPONENT_JOINT<T> comp_joint;
    comp_pipe.target_length=target_length;
    comp_change.target_length=target_length;
    comp_joint.target_length=target_length;
    
    while(getline(fin,line))
    {
        std::stringstream ss(line);
        ss>>c;
        if(isspace(c) || c=='#') continue;
        switch(c)
        {
            case 'l':
                ss>>target_length;
                comp_pipe.target_length=target_length;
                comp_change.target_length=target_length;
                comp_joint.target_length=target_length;
                break;

            case 'c':
                ss>>name>>i0>>t0;
                cross_section_hash.Set(name,{i0,t0});
                break;

            case 'v':
                ss>>name>>v0;
                vertices.Set(name,{v0});
                break;

            case 'j':
                {
                    ss>>name2>>i0>>name3;
                    auto cs=cross_section_hash.Get(name2);
                    TV O=vertices.Get(name3);
                    ARRAY<PAIR<TV,std::string> > verts;
                    for(int i=0;i<i0;i++)
                    {
                        ss>>name>>name2;
                        verts.Append({(vertices.Get(name)-O).Normalized(),name2});
                    }
                    TV first=verts(0).x;
                    auto comp=[&first](const auto& a,const auto& b)
                    {
                        T xa=TV::Oriented_Angle_Between(first,a.x);
                        if(xa<0) xa+=2*pi;
                        T xb=TV::Oriented_Angle_Between(first,b.x);
                        if(xb<0) xb+=2*pi;
                        return xa<xb;
                    };
                    std::sort(verts.begin()+1,verts.end(),comp);
                    ARRAY<T> angles;
                    for(int i=1;i<i0;i++)
                    {
                        T a=TV::Oriented_Angle_Between(verts(i-1).x,verts(i).x);
                        if(a<0) a+=2*pi;
                        angles.Append(a);
                    }
                    auto cj=comp_joint.Make_Component(cs.x,cs.y,angles);
                    XFORM<TV> xf={Compute_Xform(verts(0).x),O};
                    ARRAY<VERTEX_DATA> vd(i0);
                    Emit_Component_Blocks(cj.x,xf,vd);
                    for(int i=0;i<i0;i++)
                    {
                        vd(i).X=O+verts(i).x*cj.y(i);
                        connection_points.Set(verts(i).y,vd(i));
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
                    auto cs=cross_section_hash.Get(name);
                    TV dir=vd(1).X-vd(0).X;
                    if(!Canonical_Direction(dir))
                    {
                        dir=-dir;
                        std::swap(vd(0),vd(1));
                    }
                    T len=dir.Normalize();
                    XFORM<TV> xf={Compute_Xform(dir),vd(0).X};
                    auto cc=comp_pipe.Make_Component(cs.x,cs.y,len);
                    Emit_Component_Blocks(cc,xf,vd);
                }
                break;

            case 'g':
                {
                    ss>>name>>name2;
                    auto cs0=cross_section_hash.Get(name);
                    auto cs1=cross_section_hash.Get(name2);
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
                    auto cc=comp_change.Make_Component(cs0.x,cs0.y,cs1.x,cs1.y,t1);
                    XFORM<TV> xf={Compute_Xform(dir),A};
                    ss>>name>>name2;

                    ARRAY<VERTEX_DATA> vd(i0);
                    Emit_Component_Blocks(cc,xf,vd);
                    vd(0).X=C;
                    vd(1).X=D;
                    connection_points.Set(name,vd(0));
                    connection_points.Set(name2,vd(1));
                }
                break;
            case 'u':
            case 't':
                {
                    T len=target_length;
                    ss>>name>>name2>>name3>>name4;
                    auto cs=cross_section_hash.Get(name);
                    if(c=='u') ss>>t0;
                    else ss>>v0;
                    TV A=vertices.Get(name2);
                    TV B=vertices.Get(name3);
                    TV dir=(B-A).Normalized();
                    TV C=A+dir*len;
                    auto pr=comp_bc.Make_Block(cs.x,cs.y,len,c=='u');
                    CANONICAL_BLOCK<T>* cb=pr.x;

                    BLOCK_ID b=blocks.Append({pr.x,{Compute_Xform(dir),A},{BLOCK_ID(-1)}});
                    VERTEX_DATA vd={C,{b,CON_ID()}};
                    connection_points.Set(name4,vd);

                    BOUNDARY_CONDITION bc={b,pr.y,pr.z,{},{},-dir};
                    bc.data_v.Resize(bc.bc_v.Size());
                    bc.data_e.Resize(bc.bc_e.Size());

                    if(c=='u')
                    {
                        T y0=cb->X(bc.bc_v.min_corner).y;
                        T y1=cb->X(bc.bc_v.max_corner).y-1;
                        T a=6*t0/cube(y1-y0);
                        for(int i:bc.bc_v)
                        {
                            T y=cb->X(i).y;
                            bc.data_v(i)=a*(y-y0)*(y-y1)*dir;
                        }
                        for(int i:bc.bc_e)
                        {
                            T y=cb->X.Subset(cb->S(i)).Sum().y/2;
                            bc.data_e(i)=a*(y-y0)*(y-y1)*dir;
                        }
                        bc_v.Append(bc);
                    }
                    else
                    {
                        bc.data_v.Fill(v0);
                        bc.data_e.Fill(v0);
                        bc_t.Append(bc);
                    }
                }
                break;
            case 'U':
                delete analytic_velocity;
                analytic_velocity=new ANALYTIC_VECTOR_PROGRAM<TV>(ss.str().c_str()+2);
                break;
                break;
            case 'P':
                delete analytic_pressure;
                analytic_pressure=new ANALYTIC_SCALAR_PROGRAM<TV>(ss.str().c_str()+2);
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
Set_Connector(VERTEX_DATA& vd,BLOCK_ID id,CON_ID con_id)
{
    if(vd.con.is_regular)
    {
        if(vd.con.id>=BLOCK_ID())
        {
            blocks(vd.con.id).connections(vd.con.con_id)={id,con_id};
            blocks(id).connections(con_id)={vd.con.id,vd.con.con_id};
        }
        else
        {
            vd.con.id=id;
            vd.con.con_id=con_id;
        }
    }
    else
    {
        IRREGULAR_CONNECTION& con=irregular_connections(vd.con.irreg_id);
        con.regular=id;
        con.con_id=con_id;
        blocks(id).connections(con_id).Set_Irreg(vd.con.irreg_id);
    }
}
//#####################################################################
// Function Emit_Component_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Emit_Component_Blocks(const CANONICAL_COMPONENT<T>* cc,const XFORM<TV>& xf,ARRAY<VERTEX_DATA>& vd)
{
    int offset=Value(blocks.m),offset_edge_on=Value(irregular_connections.m);
    auto mp_bl=[=](CC_BLOCK_ID b){return BLOCK_ID(Value(b)+offset);};
    auto mp_ic=[=](CC_IRREG_ID i){return IRREG_ID(Value(i)+offset_edge_on);};

    for(CC_BLOCK_ID b(0);b<cc->blocks.m;b++)
    {
        const auto& cbl=cc->blocks(b);
        auto& bl=blocks(blocks.Add_End());
        bl.block=cbl.block;
        bl.xform=xf*cbl.xform;
        bl.flags=cbl.flags;
        bl.connections.Resize(cbl.connections.m);
        for(CON_ID c(0);c<cbl.connections.m;c++)
        {
            const auto& cn=cbl.connections(c);
            if(cn.is_regular)
            {
                if(cn.id>=CC_BLOCK_ID())
                    bl.connections(c)={mp_bl(cn.id),cn.con_id};
                else
                    Set_Connector(vd(~Value(cn.id)),mp_bl(b),c);
            }
            else
                bl.connections(c)={mp_ic(cn.irreg_id)};
        }
        for(auto e:cbl.edge_on)
            bl.edge_on.Append(mp_ic(e));
    }

    for(const auto& a:cc->irregular_connections)
    {
        IRREG_ID index=irregular_connections.Add_End();
        IRREGULAR_CONNECTION& ic=irregular_connections(index);
        for(auto& b:a.edge_on_v) ic.edge_on_v.Append({mp_bl(b.x),b.y});
        for(auto& b:a.edge_on_e) ic.edge_on_e.Append({mp_bl(b.x),b.y});

        if(a.regular>=CC_BLOCK_ID())
        {
            ic.regular=mp_bl(a.regular);
            ic.con_id=a.con_id;
        }
        else
            vd(~Value(a.regular)).con.Set_Irreg(index);
    }
}
//#####################################################################
// Function Compute_Xform
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Xform(const TV& dir) -> MATRIX<T,2>
{
    return MATRIX<T,TV::m>(dir,dir.Orthogonal_Vector());
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute()
{
    Merge_Blocks();

    for(auto& bl:blocks)
    {
        auto pr=canonical_block_matrices.Insert(bl.block,{});
        if(pr.y) Fill_Canonical_Block_Matrix(*pr.x,bl.block);
    }
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Update_Masters()
{
    for(BLOCK_ID b0(0);b0<blocks.m;b0++)
    {
        for(CON_ID c0(0);c0<blocks(b0).connections.m;c0++)
        {
            auto& n0=blocks(b0).connections(c0);
            BLOCK_ID b1=n0.id;
            CON_ID c1=n0.con_id;
            if(!n0.is_regular)
            {
                n0.master=true;
                continue;
            }
            if(b0>b1) continue; // smaller index does the logic
            auto& n1=blocks(b1).connections(c1);
            bool m=true;
            if(c0==CON_ID())
            {
                if(c1!=CON_ID()) m=true;
                else
                {
                    if(blocks(b0).connections.m==CON_ID(2)) m=true;
                    else if(blocks(b1).connections.m==CON_ID(2)) m=false;
                }
            }
            else
            {
                if(!c1) m=false;
                else if(blocks(b0).connections.m==CON_ID(2)) m=true;
                else if(blocks(b1).connections.m==CON_ID(2)) m=false;
            }
            n0.master=m;
            n1.master=!m;
        }
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
        int m=n/2+(n&m0);
        for(int i=0;i<m;i++) func(c+i,a+b*i,true);
        for(int i=m;i<n;i++) func(c+i,a+b*i,false);
    }
    else
    {
        int m=n/2+(n&(1-m0));
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
    CANONICAL_BLOCK<T>* cb=blocks(b).block;
    int mask=0;
    for(CON_ID c(0);c<blocks(b).connections.m;c++)
        if(blocks(b).connections(c).master)
            mask|=1<<Value(c);
    auto pr=separates_dofs.Insert({cb,mask},0);
    if(!pr.y) return *pr.x;

    int valid=0;
    ARRAY<CON_ID> owner(cb->X.m,use_init,CON_ID(-1));
    for(CON_ID i(0);i<cb->cross_sections.m;i++)
    {
        const auto& cs=cb->cross_sections(i);
        int master=(mask>>Value(i))&1;
        INTERVAL<int> iv=cs.v;
        int m=(iv.min_corner+iv.max_corner)/2;
        if(cs.own_first) iv.min_corner=m+(iv.Size()&master);
        else iv.max_corner=m+(iv.Size()&(1-master));
        for(int j:iv) owner(j)=i;
    }
    for(auto e:cb->E)
    {
        VECTOR<CON_ID,3> f(owner.Subset(e));
        CON_ID o(-1);
        for(auto i:f)
        {
            if(i==CON_ID(-1) || i==o) continue;
            if(o==CON_ID(-1))
            {
                o=i;
                continue;
            }
            valid|=1<<Value(o);
            valid|=1<<Value(i);
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
Merge_Canonical_Blocks(CANONICAL_BLOCK<T>* cb0,CON_ID con_id0,XFORM<TV> xf0,
    CANONICAL_BLOCK<T>* cb1,CON_ID con_id1,XFORM<TV> xf1) -> PAIR<CANONICAL_BLOCK<T>*,ARRAY<int> >*
{
    auto pr=merge_canonical_blocks.Insert(std::make_tuple(cb0,con_id0,cb1,con_id1),{});
    if(pr.y) return pr.x;

    const CROSS_SECTION& cs0=cb0->cross_sections(con_id0);
    const CROSS_SECTION& cs1=cb1->cross_sections(con_id1);

    ARRAY<int> index_v_map,index_e_map;
    int num_v_dofs=Merge_Dofs(index_v_map,cb0->X.m,cb1->X.m,cs0.v,cs1.v,
        cs0.own_first,cs1.own_first);
    int num_e_dofs=Merge_Dofs(index_e_map,cb0->X.m,cb1->X.m,cs0.e,cs1.e,
        cs0.own_first,cs1.own_first);

    pr.x->x=new CANONICAL_BLOCK<T>;
    CANONICAL_BLOCK<T>* cb=pr.x->x;
    for(CON_ID i(0);i<cb0->cross_sections.m;i++)
        if(i!=con_id0)
            cb->cross_sections.Append(cb0->cross_sections(i));
    for(CON_ID i(0);i<cb1->cross_sections.m;i++)
        if(i!=con_id1)
            cb->cross_sections.Append(
                Map_Cross_Section(cb1->cross_sections(i),index_v_map,index_e_map));

    XFORM<TV> M01i=xf0*xf1.Inverse();

    cb->X=cb0->X;
    cb->X.Resize(num_v_dofs);
    for(int i=0;i<cb1->X.m;i++)
    {
        TV X=M01i*cb1->X(i);
        int j=index_v_map(i);
        if(j>=cb0->X.m) cb->X(j)=X;
        else assert((cb->X(j)-X).Magnitude()<(T)1e-6);
    }

    cb->S=cb0->S;
    cb->S.Resize(num_e_dofs);
    cb->bc_v=cb0->bc_v;
    cb->bc_e=cb0->bc_e;
    for(int i=0;i<cb1->E.m;i++)
    {
        int j=index_e_map(i);
        IV s(index_v_map.Subset(cb1->E(i)));
        s.Sort();
        if(j>=cb0->S.m) cb->S(j)=s;
        else assert(cb->S(j)==s);
    }

    cb->bc_v.Append_Elements(index_v_map.Subset(cb1->bc_v));
    cb->bc_e.Append_Elements(index_e_map.Subset(cb1->bc_e));

    cb->E=cb0->E;
    ARRAY<IV3> E=cb1->E;
    E.Flattened()=index_v_map.Subset(E.Flattened());
    cb->E.Append_Elements(E);

    return pr.x;
}
//#####################################################################
// Function Merge_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Merge_Blocks(BLOCK_ID id,CON_ID con_id)
{
    PHYSBAM_ASSERT(con_id>=CON_ID());
    PHYSBAM_ASSERT(id>=BLOCK_ID());
    BLOCK& bl=blocks(id);
    BLOCK_ID id2=bl.connections(con_id).id;
    CON_ID con_id2=bl.connections(con_id).con_id;
    BLOCK& bl2=blocks(id2);

    auto pr=Merge_Canonical_Blocks(bl.block,con_id,bl.xform,bl2.block,con_id2,bl2.xform);
    bl.block=pr->x;

    for(CON_ID c=con_id;c<bl.connections.m;c++)
    {
        const auto& cn=bl.connections(c);
        if(cn.is_regular) blocks(cn.id).connections(cn.con_id).con_id=c-1;
        else irregular_connections(cn.irreg_id).con_id=c-1;
    }
    bl.connections.Pop();

    for(CON_ID c(0);c<bl2.connections.m;c++)
    {
        if(c==con_id2) continue;
        auto& cn=bl2.connections(c);
        if(cn.is_regular)
        {
            blocks(cn.id).connections(cn.con_id).id=id;
            blocks(cn.id).connections(cn.con_id).con_id=bl.connections.m;
        }
        else
        {
            irregular_connections(cn.irreg_id).regular=id;
            irregular_connections(cn.irreg_id).con_id=bl.connections.m;
        }
        bl.connections.Append(cn);
    }

    for(IRREG_ID i:bl2.edge_on)
    {
        auto& ic=irregular_connections(i);
        for(auto& j:ic.edge_on_v)
            if(j.x==id2)
            {
                j.x=id;
                j.y=pr->y(j.y);
            }
    }

    bl2={0};
}
//#####################################################################
// Function Merge_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Merge_Blocks()
{
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        if(!blocks(b).block) continue;
        if(int mask=Separates_Dofs(b))
        {
            PHYSBAM_ASSERT(!(blocks(b).flags&1));
            CON_ID besti(-1);
            int best=INT_MAX;
            for(CON_ID i(0);i<blocks(b).connections.m;i++)
            {
                if(!(mask&(1<<Value(i)))) continue;
                if(!blocks(b).connections(i).is_regular) continue;
                BLOCK_ID d=blocks(b).connections(i).id;
                if(blocks(d).flags&1) continue;
                int c=Approx_Dof_Count(d);
                if(c<best) // keep blocks small
                {
                    best=c;
                    besti=i;
                }
            }
            PHYSBAM_ASSERT(besti>=CON_ID());
            Merge_Blocks(b,besti);
            b--; // repeat the check on this block
        }
    }
    // TODO: renumber the BLOCK_ID's
}
//#####################################################################
// Function Approx_Dof_Count
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Approx_Dof_Count(BLOCK_ID b)
{
    const auto& bl=blocks(b);
    const auto* cb=bl.block;
    int num=cb->X.m;
    for(const auto& cs:cb->cross_sections)
        num-=cs.v.Size()/2;
    return num;
}


int fem_visc_table[12][12]=
{
    {3,1,0,0,0,-4,3,1,0,0,0,-4},
    {1,3,0,0,0,-4,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,-1,0,4,-4,0},
    {0,0,0,8,-8,0,0,4,0,4,-4,-4},
    {0,0,0,-8,8,0,-4,0,0,-4,4,4},
    {-4,-4,0,0,0,8,0,-4,0,-4,4,4},
    {3,0,1,0,-4,0,3,0,1,0,-4,0},
    {1,0,-1,4,0,-4,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,3,0,-4,0},
    {0,0,4,4,-4,-4,0,0,0,8,0,-8},
    {0,0,-4,-4,4,4,-4,0,-4,0,8,0},
    {-4,0,0,-4,4,4,0,0,0,-8,0,8}
};

int fem_pres_table[3][12]=
{
    {-1,0,0,1,-1,1,-1,0,0,1,1,-1},
    {0,1,0,1,-1,-1,0,0,0,2,0,-2},
    {0,0,0,2,-2,0,0,0,1,1,-1,-1}
};

//#####################################################################
// Function Fill_Canonical_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Canonical_Block_Matrix(BLOCK_MATRIX<T>& mat,const CANONICAL_BLOCK<T>* cb)
{
    mat.nv_r=mat.nv_c=cb->X.m;
    mat.ne_r=mat.ne_c=cb->S.m;
    mat.np_r=mat.np_c=cb->X.m;
    mat.Resize();

    HASHTABLE<IV,int> edge_lookup;
    for(int i=0;i<cb->S.m;i++)
        edge_lookup.Set(cb->S(i).Sorted(),i);

    for(IV3 v:cb->E)
    {
        IV3 dof[2]={v};
        for(int i=0;i<3;i++)
            dof[1](i)=edge_lookup.Get(v.Remove_Index(i).Sorted());
        MATRIX<T,2> F(cb->X(v.y)-cb->X(v.x),cb->X(v.z)-cb->X(v.x)),G=F.Inverse();
        T scale=mu*F.Determinant()/6;
        T p_scale=F.Determinant()/6;

        for(int r=0;r<2;r++)
            for(int s=0;s<2;s++)
                for(int i=0;i<3;i++)
                    for(int j=0;j<3;j++)
                    {
                        MATRIX<T,2> M;
                        for(int a=0;a<2;a++)
                            for(int b=0;b<2;b++)
                                M(a,b)=fem_visc_table[6*a+3*r+i][6*b+3*s+j];
                        M=scale*G.Transpose_Times(M*G);
                        mat.Add_uu(dof[r](i),r,dof[s](j),s,M+M.Trace());
                    }

        for(int s=0;s<2;s++)
            for(int i=0;i<3;i++)
                for(int j=0;j<3;j++)
                {
                    TV u;
                    for(int b=0;b<2;b++)
                        u(b)=fem_pres_table[i][6*b+3*s+j];
                    u=p_scale*G.Transpose_Times(u);
                    mat.Add_pu(v(i),dof[s](j),s,u);
                    mat.Add_up(dof[s](j),s,v(i),u);
                }
    }
}


// entry * J / 360
int fem_u_dot_v_table[6][6]=
{
    {6,-1,-1,-4,0,0},
    {-1,6,-1,0,-4,0},
    {-1,-1,6,0,0,-4},
    {-4,0,0,32,16,16},
    {0,-4,0,16,32,16},
    {0,0,-4,16,16,32}
};

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Times_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& w,const BLOCK_VECTOR<T>& u) const
{
    const auto& bl=blocks(b);
    const auto* cb=bl.block;
    MATRIX<T,2> M=bl.xform.M;

    HASHTABLE<IV,int> edge_lookup;
    for(int i=0;i<cb->S.m;i++)
        edge_lookup.Set(cb->S(i).Sorted(),i);

    for(IV3 v:cb->E)
    {
        IV3 dof[2]={v};
        for(int i=0;i<3;i++)
            dof[1](i)=edge_lookup.Get(v.Remove_Index(i).Sorted());
        MATRIX<T,2> F(cb->X(v.y)-cb->X(v.x),cb->X(v.z)-cb->X(v.x));
        T scale=(M*F).Determinant()/360;

        VECTOR<TV,6> r,s;
        for(int a=0;a<2;a++) for(int i=0;i<3;i++) r(i+3*a)=u.Get_u(dof[a](i),a);
        for(int i=0;i<6;i++)
            for(int j=0;j<6;j++)
                s(i)+=r(j)*fem_u_dot_v_table[i][j];
        s*=scale;
        for(int a=0;a<2;a++) for(int i=0;i<3;i++) w.Add_u(dof[a](i),a,s(i+3*a));
    }
}

int fem_p_u_table[3][6]=
{
    {2, -1, -1, 4, 8, 8},
    {-1, 2, -1, 8, 4, 8},
    {-1, -1, 2, 8, 8, 4}
};

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Times_P_U(BLOCK_ID b,BLOCK_VECTOR<T>& w,const ARRAY<T>& div_v,const ARRAY<T>& div_e) const
{
    const auto& bl=blocks(b);
    const auto* cb=bl.block;
    MATRIX<T,2> M=bl.xform.M;

    HASHTABLE<IV,int> edge_lookup;
    for(int i=0;i<cb->S.m;i++)
        edge_lookup.Set(cb->S(i).Sorted(),i);

    for(IV3 v:cb->E)
    {
        IV3 dof[2]={v};
        for(int i=0;i<3;i++)
            dof[1](i)=edge_lookup.Get(v.Remove_Index(i).Sorted());
        MATRIX<T,2> F(cb->X(v.y)-cb->X(v.x),cb->X(v.z)-cb->X(v.x));
        T scale=(M*F).Determinant()/120;

        VECTOR<T,6> r;
        VECTOR<T,3> s;
        for(int i=0;i<3;i++) r(i)=div_v(dof[0](i));
        for(int i=0;i<3;i++) r(i+3)=div_e(dof[1](i));
        for(int i=0;i<3;i++)
            for(int j=0;j<6;j++)
                s(i)+=r(j)*fem_p_u_table[i][j];
        s*=scale;
        for(int i=0;i<3;i++) w.Add_p(v(i),s(i));
    }
}

int fem_line_int_u_dot_v_table[3][3] = {{4, -1, 2}, {-1, 4, 2}, {2, 2, 16}};

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Times_Line_Integral_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& w,const BLOCK_VECTOR<T>& u) const
{
    const auto& bl=blocks(b);
    const auto* cb=bl.block;
    MATRIX<T,2> M=bl.xform.M;
    for(int e:cb->bc_e)
    {
        VECTOR<TV,3> r(u.Get_v(cb->S(e).x),u.Get_v(cb->S(e).y),u.Get_e(e)),s;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                s(i)+=r(j)*fem_line_int_u_dot_v_table[i][j];

        VECTOR<TV,2> X(cb->X.Subset(cb->S(e)));
        T scale=(M*(X.x-X.y)).Magnitude()/30;
        w.Add_v(cb->S(e).x,s.x*scale);
        w.Add_v(cb->S(e).y,s.y*scale);
        w.Add_e(e,s.z*scale);
    }
}

//#####################################################################
// Function Copy_Matrix_Data
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Copy_Matrix_Data(BLOCK_MATRIX<T>& A,BLOCK_ID b,
    const DOF_PAIRS& dpa,const DOF_PAIRS& dpb,BLOCK_ID ar,BLOCK_ID ac) const
{
    const BLOCK_MATRIX<T>& B=canonical_block_matrices.Get(blocks(b).block);
    MATRIX<T,2> G=blocks(b).xform.M.Inverse();
    MATRIX<T,2> Ma=G*blocks(ar).xform.M;
    MATRIX<T,2> Mb=G*blocks(ac).xform.M;
    T sa=1/sqrt(Ma.Determinant());
    T sb=1/sqrt(Mb.Determinant());
    Ma*=sa;
    Mb*=sb;

    for(auto p:dpa.v)
    {
        for(auto r:dpb.v) A.Add_vv(p.x,r.x,Ma.Transpose_Times(B.Get_vv(p.y,r.y)*Mb));
        for(auto r:dpb.e) A.Add_ve(p.x,r.x,Ma.Transpose_Times(B.Get_ve(p.y,r.y)*Mb));
        for(auto r:dpb.p) A.Add_vp(p.x,r.x,Ma.Transpose_Times(B.Get_vp(p.y,r.y)*sb));
    }
    for(auto p:dpa.e)
    {
        for(auto r:dpb.v) A.Add_ev(p.x,r.x,Ma.Transpose_Times(B.Get_ev(p.y,r.y)*Mb));
        for(auto r:dpb.e) A.Add_ee(p.x,r.x,Ma.Transpose_Times(B.Get_ee(p.y,r.y)*Mb));
        for(auto r:dpb.p) A.Add_ep(p.x,r.x,Ma.Transpose_Times(B.Get_ep(p.y,r.y)*sb));
    }
    for(auto p:dpa.p)
    {
        for(auto r:dpb.v) A.Add_pv(p.x,r.x,Mb.Transpose_Times(B.Get_pv(p.y,r.y)*sa));
        for(auto r:dpb.e) A.Add_pe(p.x,r.x,Mb.Transpose_Times(B.Get_pe(p.y,r.y)*sa));
    }
}
//#####################################################################
// Function Copy_Vector_Data
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Copy_Vector_Data(const BLOCK_VECTOR<T>& B,BLOCK_ID b,const DOF_PAIRS& dp,BLOCK_ID a)
{
    BLOCK_VECTOR<T>& A=rhs_block_list(b);
    if(!A.V.m) Init_Block_Vector(A,b);

    MATRIX<T,2> M=blocks(b).xform.M;
    T s=1/sqrt(M.Determinant());
    M*=s;

    for(auto p:dp.v) A.Add_v(p.x,M.Transpose_Times(B.Get_v(p.y)));
    for(auto p:dp.e) A.Add_e(p.x,M.Transpose_Times(B.Get_e(p.y)));
    for(auto p:dp.p) A.Add_p(p.x,s*B.Get_p(p.y));
}
//#####################################################################
// Function Fill_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Block_Matrix(REFERENCE_BLOCK_DATA& rd)
{
    BLOCK_ID b=rd.b;
    const auto& bl=blocks(b);
    Init_Block_Matrix(rd.M,b,b);

    Copy_Matrix_Data(rd.M,b,rd.pairs,rd.pairs,b,b);

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
        {
            const auto& p=reference_block_data(blocks(c.id).ref_id).regular_pairs(c.con_id);
            Copy_Matrix_Data(rd.M,c.id,p,p,b,b);
        }
        else
        {
            const auto& ic=irregular_connections(c.irreg_id);
            const auto& irbd=reference_irregular_data(ic.ref_id);
            for(const auto& h:irbd.pairs)
                Copy_Matrix_Data(rd.M,h.b,h.irreg_pairs[0][1],h.irreg_pairs[0][1],b,b);
        }
    }

    for(IRREG_ID e:bl.edge_on)
    {
        const auto& ic=irregular_connections(e);
        const auto& irbd=reference_irregular_data(ic.ref_id);
        for(const auto& h:irbd.pairs)
            if(h.b==b)
                Copy_Matrix_Data(rd.M,h.b,h.irreg_pairs[1][0],h.irreg_pairs[1][0],ic.regular,ic.regular);
    }
}
//#####################################################################
// Function Fill_Connection_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Connection_Matrix(REFERENCE_CONNECTION_DATA& cd)
{
    Init_Block_Matrix(cd.M,cd.b[0],cd.b[1]);
    PHYSBAM_ASSERT(blocks(cd.b[0]).connections(cd.con_id[0]).master);
    auto& rd0=reference_block_data(blocks(cd.b[0]).ref_id);
    auto& rd1=reference_block_data(blocks(cd.b[1]).ref_id);
    Copy_Matrix_Data(cd.M,cd.b[0],rd0.pairs,rd0.regular_pairs(cd.con_id[0]),cd.b[0],cd.b[1]);
    Copy_Matrix_Data(cd.M,cd.b[1],rd1.regular_pairs(cd.con_id[1]),rd1.pairs,cd.b[0],cd.b[1]);
}
//#####################################################################
// Function Fill_Irregular_Connection_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Fill_Irregular_Connection_Matrix(REFERENCE_IRREGULAR_DATA& ri)
{
    BLOCK_ID bb=irregular_connections(ri.ic_id).regular;
    for(auto&z:ri.pairs)
    {
        Init_Block_Matrix(z.M,bb,z.b);
        BLOCK_ID b[2]={bb,z.b};
        for(int i=0;i<2;i++)
            Copy_Matrix_Data(z.M,b[i],z.irreg_pairs[0][i],z.irreg_pairs[1][i],bb,z.b);
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
    Compute_Reference_Regular_Connections();
    Compute_Reference_Irregular_Connections();

    for(auto& rd:reference_block_data)
    {
        Compute_Dof_Remapping(rd);
        Compute_Dof_Pairs(rd);
    }

    for(auto& rc:reference_connection_data)
        Compute_Dof_Pairs(rc);

    for(auto& ri:reference_irregular_data)
        Compute_Dof_Pairs(ri);

    for(auto& rd:reference_block_data)
    {
        Fill_Block_Matrix(rd);
        rd.mat_id=cem.Create_Matrix_Block(true);
        cem.block_list(rd.mat_id).M=rd.M.M;
    }

    for(auto& rc:reference_connection_data)
    {
        Fill_Connection_Matrix(rc);
        rc.mat_id=cem.Create_Matrix_Block(false);
        cem.block_list(rc.mat_id).M=rc.M.M;
    }

    for(auto& ri:reference_irregular_data)
    {
        Fill_Irregular_Connection_Matrix(ri);
        for(auto& d:ri.pairs)
        {
            d.mat_id=cem.Create_Matrix_Block(false);
            cem.block_list(d.mat_id).M=d.M.M;
        }
    }

    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        int id=reference_block_data(blocks(b).ref_id).mat_id;
        cem.Add_Block_Matrix_Entry(Value(b),Value(b),id);
        nonzero_blocks.Append({b,b,id});
    }

    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        const auto& bl=blocks(b);
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular && c.master)
            {
                auto reg_id=regular_connection_hash.Get(std::make_tuple(bl.ref_id,cc,blocks(c.id).ref_id,c.con_id));
                int id=reference_connection_data(reg_id).mat_id;
                cem.Add_Block_Matrix_Entry(Value(b),Value(c.id),id);
                nonzero_blocks.Append({b,c.id,id});
            }
        }
    }

    for(auto& ic:irregular_connections)
    {
        auto &ri=reference_irregular_data(ic.ref_id);
        for(const auto& p:ri.pairs)
        {
            cem.Add_Block_Matrix_Entry(Value(ic.regular),Value(p.b),p.mat_id);
            nonzero_blocks.Append({ic.regular,p.b,p.mat_id});
        }
    }

    rhs_block_list.Resize(blocks.m);

    if(analytic_velocity && analytic_pressure)
    {
        for(BLOCK_ID b(0);b<blocks.m;b++)
        {
            BLOCK& bl=blocks(b);
            CANONICAL_BLOCK<T>* cb=bl.block;

            BLOCK_VECTOR<T> w,u;
            Init_Block_Vector(w,cb);
            Init_Block_Vector(u,cb);
            for(auto i:cb->bc_v)
                u.Add_v(i,analytic_velocity->v(cb->X(i),0));
            for(auto i:cb->bc_e)
                u.Add_e(i,analytic_velocity->v(cb->X.Subset(cb->S(i)).Average(),0));
            w.V=canonical_block_matrices.Get(bl.block).M*u.V;

            u.V.Fill(0);
            ARRAY<T> div_v(cb->X.m),div_e(cb->S.m);
            for(int i=0;i<cb->X.m;i++)
            {
                TV Z=cb->X(i);
                u.Add_v(i,Force(Z));
                div_v(i)=analytic_velocity->dX(Z,0).Trace();
            }
            for(int i=0;i<cb->S.m;i++)
            {
                TV Z=cb->X.Subset(cb->S(i)).Average();
                u.Add_e(i,Force(Z));
                div_e(i)=analytic_velocity->dX(Z,0).Trace();
            }
            Times_U_Dot_V(b,w,u);
            Times_P_U(b,w,div_v,div_e);
            Apply_To_RHS(b,w);
        }

        for(const auto& bc:bc_t)
        {
            BLOCK& bl=blocks(bc.b);
            CANONICAL_BLOCK<T>* cb=bl.block;

            BLOCK_VECTOR<T> w,u;
            Init_Block_Vector(w,cb);
            Init_Block_Vector(u,cb);
            for(auto i:cb->bc_v)
            {
                TV Z=cb->X(i);
                u.Add_v(i,Traction(bc.normal,Z));
            }
            for(auto i:cb->bc_e)
            {
                TV Z=cb->X.Subset(cb->S(i)).Average();
                u.Add_e(i,Traction(bc.normal,Z));
            }
            Times_Line_Integral_U_Dot_V(bc.b,w,u);
            Apply_To_RHS(bc.b,w);
        }
    }
    else
    {
        for(const auto& bc:bc_v)
        {
            BLOCK& bl=blocks(bc.b);
            CANONICAL_BLOCK<T>* cb=bl.block;

            BLOCK_VECTOR<T> w,u;
            Init_Block_Vector(w,cb);
            Init_Block_Vector(u,cb);
            for(int i=0;i<bc.bc_v.Size();i++) u.Add_v(bc.bc_v.min_corner+i,bc.data_v(i));
            for(int i=0;i<bc.bc_e.Size();i++) u.Add_e(bc.bc_e.min_corner+i,bc.data_e(i));
            w.V=canonical_block_matrices.Get(bl.block).M*u.V;
            Apply_To_RHS(bc.b,w);
        }

        for(const auto& bc:bc_t)
        {
            BLOCK& bl=blocks(bc.b);
            CANONICAL_BLOCK<T>* cb=bl.block;

            BLOCK_VECTOR<T> w,u;
            Init_Block_Vector(w,cb);
            Init_Block_Vector(u,cb);
            for(int i=0;i<bc.bc_v.Size();i++) u.Add_v(bc.bc_v.min_corner+i,bc.data_v(i));
            for(int i=0;i<bc.bc_e.Size();i++) u.Add_e(bc.bc_e.min_corner+i,bc.data_e(i));
            Times_Line_Integral_U_Dot_V(bc.b,w,u);
            Apply_To_RHS(bc.b,w);
        }
    }

    for(BLOCK_ID i(0);i<rhs_block_list.m;i++)
    {
        int j=-1;
        if(rhs_block_list(i).V.m)
            j=cem.vector_list.Append(rhs_block_list(i).V);
        cem.rhs(Value(i))=j;
    }

    // TODO: rhs
    // multiply rhs(b) by F.Transpose(), where
    // A = blocks(b).xform.M;
    // F = A/sqrt(A.Determinant())
    //
    // Multiply solution sol(b) by F when done.

    cem.End_Fill_Blocks();
    cem.valid_row.Resize(Value(blocks.m),use_init,true);

    LOG::printf("CM\n");
    for(auto&p:canonical_block_matrices) LOG::printf("%P %P\n",std::make_tuple(p.data.nv_r,p.data.ne_r,p.data.np_r,p.data.nv_c,p.data.ne_c,p.data.np_c),p.data.M);

    for(int i=0;i<cem.block_list.m;i++)
        LOG::printf("%i %P\n",i,cem.block_list(i).M);
}
//#####################################################################
// Function Compute_Block_Hash
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Reference_Blocks()
{
    typedef ARRAY<std::tuple<CANONICAL_BLOCK<T>*,CON_ID,bool> > REG_CON;
    typedef ARRAY<std::tuple<int,CANONICAL_BLOCK<T>*,int> > IRREG_CON_V;
    typedef ARRAY<std::tuple<int,CANONICAL_BLOCK<T>*,int> > IRREG_CON_E;
    typedef std::tuple<CANONICAL_BLOCK<T>*,REG_CON,IRREG_CON_V,IRREG_CON_E> KEY;

    HASHTABLE<KEY,REFERENCE_BLOCK_ID> h;

    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        REG_CON reg_con;
        IRREG_CON_V irreg_con_v;
        IRREG_CON_E irreg_con_e;

        auto& bl=blocks(b);
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular)
                reg_con.Append(std::make_tuple(blocks(c.id).block,c.con_id,c.master));
            else
            {
                reg_con.Append(std::make_tuple((CANONICAL_BLOCK<T>*)0,CON_ID(-1),true));

                const auto& ic=irregular_connections(c.irreg_id);
                Visit_Irregular_Cross_Section_Dofs(ic,
                    blocks(ic.regular).block->cross_sections(ic.con_id),
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
        {
            *pr.x=reference_block_data.Add_End();
            auto& rd=reference_block_data(*pr.x);
            rd.b=b;
            rd.regular_pairs.Resize(bl.connections.m);
        }
        bl.ref_id=*pr.x;
    }
}
//#####################################################################
// Function Compute_Reference_Regular_Connections
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Reference_Regular_Connections()
{
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        const auto& bl=blocks(b);
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular && c.master)
            {
                auto pr=regular_connection_hash.Insert(std::make_tuple(bl.ref_id,cc,blocks(c.id).ref_id,c.con_id),{});
                if(pr.y) *pr.x=reference_connection_data.Append({{b,c.id},{cc,c.con_id}});
            }
        }
    }
}
//#####################################################################
// Function Compute_Reference_Irregular_Connections
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Reference_Irregular_Connections()
{
    typedef std::tuple<CANONICAL_BLOCK<T>*,CON_ID,ARRAY<PAIR<CANONICAL_BLOCK<T>*,int> >,ARRAY<PAIR<CANONICAL_BLOCK<T>*,int> > > KEY;
    HASHTABLE<KEY,IRREG_ID> ref;

    for(IRREG_ID i(0);i<irregular_connections.m;i++)
    {
        auto& ic=irregular_connections(i);
        ARRAY<PAIR<CANONICAL_BLOCK<T>*,int> > av,ae;
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
Compute_Connection_Hash(BLOCK_ID b0,CON_ID con_id0,BLOCK_ID b1,CON_ID con_id1)
{
    auto id0=blocks(b0).block,id1=blocks(b1).block;
    return Hash(std::make_tuple(id0,con_id0,id1,con_id1));
}
//#####################################################################
// Function Compute_Dof_Remapping
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Dof_Remapping(REFERENCE_BLOCK_DATA& rd)
{
    const auto& bl=blocks(rd.b);
    const auto* cb=bl.block;
    rd.dof_map_v.Resize(cb->X.m,use_init,1);
    rd.dof_map_e.Resize(cb->S.m,use_init,1);

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
        {
            const auto* cb2=blocks(c.id).block;
            Visit_Regular_Cross_Section_Dofs(cb->cross_sections(cc),
                cb2->cross_sections(c.con_id),c.master,
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
            const auto& ic=irregular_connections(c.irreg_id);
            Visit_Irregular_Cross_Section_Dofs(ic,
                blocks(ic.regular).block->cross_sections(ic.con_id),
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

    int iv=0,ie=0,ip=0;
    rd.dof_map_p=rd.dof_map_v;
    rd.dof_map_v.Subset(cb->bc_v).Fill(0);
    rd.dof_map_e.Subset(cb->bc_e).Fill(0);
    for(int i=0;i<rd.dof_map_v.m;i++)
        rd.dof_map_v(i)=rd.dof_map_v(i)?iv++:-1;
    for(int i=0;i<rd.dof_map_e.m;i++)
        rd.dof_map_e(i)=rd.dof_map_e(i)?ie++:-1;
    for(int i=0;i<rd.dof_map_p.m;i++)
        rd.dof_map_p(i)=rd.dof_map_p(i)?ip++:-1;
    rd.num_dofs_v=iv;
    rd.num_dofs_e=ie;
    rd.num_dofs_p=ip;
}
//#####################################################################
// Function Init_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Init_Block_Matrix(BLOCK_MATRIX<T>& M,BLOCK_ID a,BLOCK_ID b) const
{
    const auto& c=reference_block_data(blocks(a).ref_id);
    const auto& d=reference_block_data(blocks(b).ref_id);
    M.nv_r=c.num_dofs_v;
    M.ne_r=c.num_dofs_e;
    M.np_r=c.num_dofs_p;
    M.nv_c=d.num_dofs_v;
    M.ne_c=d.num_dofs_e;
    M.np_c=d.num_dofs_p;
    M.Resize();
}
//#####################################################################
// Function Init_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Init_Block_Vector(BLOCK_VECTOR<T>& V,BLOCK_ID b) const
{
    const auto& c=reference_block_data(blocks(b).ref_id);
    V.nv=c.num_dofs_v;
    V.ne=c.num_dofs_e;
    V.np=c.num_dofs_p;
    V.Resize();
}
//#####################################################################
// Function Init_Block_Matrix
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Init_Block_Vector(BLOCK_VECTOR<T>& V,const CANONICAL_BLOCK<T>* cb) const
{
    V.nv=cb->X.m;
    V.ne=cb->S.m;
    V.np=cb->X.m;
    V.Resize();
}
//#####################################################################
// Function Apply_To_RHS
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Apply_To_RHS(BLOCK_ID b,const BLOCK_VECTOR<T>& w)
{
    const auto& rd=reference_block_data(blocks(b).ref_id);
    const auto& bl=blocks(b);

    Copy_Vector_Data(w,b,rd.pairs,b);

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
            Copy_Vector_Data(w,b,rd.regular_pairs(cc),c.id);
        else
        {
            const auto& ic=irregular_connections(c.irreg_id);
            const auto& irbd=reference_irregular_data(ic.ref_id);
            for(const auto& h:irbd.pairs)
                Copy_Vector_Data(w,b,h.irreg_pairs[1][0],h.b);
        }
    }

    for(IRREG_ID e:bl.edge_on)
    {
        const auto& ic=irregular_connections(e);
        const auto& irbd=reference_irregular_data(ic.ref_id);
        for(const auto& h:irbd.pairs)
            if(h.b==b)
                Copy_Vector_Data(w,b,h.irreg_pairs[0][1],ic.regular);
    }
}
//#####################################################################
// Function Compute_Bounding_Box
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Bounding_Box() const -> RANGE<TV>
{
    RANGE<TV> box=RANGE<TV>::Empty_Box();
    for(const auto&bl:blocks)
        for(auto i:bl.block->bc_v)
            box.Enlarge_To_Include_Point(bl.xform*bl.block->X(i));
    return box;
}
//#####################################################################
// Function Eliminate_Strip
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Eliminate_Strip(CACHED_ELIMINATION_MATRIX<T>& cem,const ARRAY<BLOCK_ID>& a)
{
    for(int sep=1;sep-1<a.m;sep*=2)
        for(int i=sep-1;i<a.m;i+=sep*2)
            if(cem.valid_row(Value(a(i))))
                cem.Eliminate_Row(Value(a(i)));
}
//#####################################################################
// Function Eliminate_Strip
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Eliminate_Simple(CACHED_ELIMINATION_MATRIX<T>& cem,BLOCK_ID first,CON_ID con_id_source)
{
    ARRAY<BLOCK_ID> list;
    BLOCK_ID b=first;
    CON_ID con_id=con_id_source;
    while(1)
    {
        const auto& bl=blocks(b);
        if(bl.connections.m>CON_ID(2)) break;
        if(bl.flags&1) break;
        list.Append(b);
        if(bl.connections.m==CON_ID(1)) break;
        CON_ID o(1-Value(con_id));
        b=bl.connections(o).id;
        con_id=bl.connections(o).con_id;
    }
    Eliminate_Strip(cem,list);
}
//#####################################################################
// Function Eliminate_Rows
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Eliminate_Irregular_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem)
{
    // Get rid of irregular connections first.
    for(const auto& ic:irregular_connections)
    {
        BLOCK_ID prev(-1);
        ARRAY<BLOCK_ID> found;
        for(const auto& p:ic.edge_on_v)
        {
            if(p.x!=prev)
            {
                found.Append(p.x);
                prev=p.x;
            }
        }
        Eliminate_Strip(cem,found);
    }
}
//#####################################################################
// Function Eliminate_Irregular_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Eliminate_Non_Seperators(CACHED_ELIMINATION_MATRIX<T>& cem)
{
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        if(!cem.valid_row(Value(b))) continue;
        const auto& bl=blocks(b);
        if(bl.connections.m>CON_ID(2) || bl.flags&1)
        {
            for(const auto& c:bl.connections)
                Eliminate_Simple(cem,c.id,c.con_id);
            if(bl.connections.m>CON_ID(2))
                cem.Eliminate_Row(Value(b));
        }
    }
    for(BLOCK_ID b(0);b<blocks.m;b++)
        assert(!cem.valid_row(Value(b)) || blocks(b).flags&1);
}
//#####################################################################
// Function Compute_Dof_Pairs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Dof_Pairs(REFERENCE_BLOCK_DATA& rd)
{
    // Copy from self
    for(int i=0;i<rd.dof_map_v.m;i++)
        if(rd.dof_map_v(i)>=0)
            rd.pairs.v.Append({rd.dof_map_v(i),i});
    for(int i=0;i<rd.dof_map_e.m;i++)
        if(rd.dof_map_e(i)>=0)
            rd.pairs.e.Append({rd.dof_map_e(i),i});
    for(int i=0;i<rd.dof_map_p.m;i++)
        if(rd.dof_map_p(i)>=0)
            rd.pairs.p.Append({rd.dof_map_p(i),i});
}
//#####################################################################
// Function Compute_Dof_Pairs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Dof_Pairs(REFERENCE_CONNECTION_DATA& rc)
{
    const auto* cb0=blocks(rc.b[0]).block;
    const auto* cb1=blocks(rc.b[1]).block;
    REFERENCE_BLOCK_DATA* rd[2];
    rd[0]=&reference_block_data(blocks(rc.b[0]).ref_id);
    rd[1]=&reference_block_data(blocks(rc.b[1]).ref_id);
    DOF_PAIRS* pairs[2];
    pairs[0]=&rd[0]->regular_pairs(rc.con_id[0]);
    pairs[1]=&rd[1]->regular_pairs(rc.con_id[1]);
    Visit_Regular_Cross_Section_Dofs(cb0->cross_sections(rc.con_id[0]),
        cb1->cross_sections(rc.con_id[1]),true, // 0 is master
        [&](int a,int b,bool o)
        {
            int f[2]={a,b};
            int i=!o,r=rd[i]->dof_map_v(f[i]),s=rd[i]->dof_map_p(f[i]);
            if(r>=0) pairs[o]->v.Append({r,f[o]});
            if(s>=0) pairs[o]->p.Append({s,f[o]});
        },
        [&](int a,int b,bool o)
        {
            int f[2]={a,b};
            int i=!o,r=rd[i]->dof_map_e(f[i]);
            if(r>=0) pairs[o]->e.Append({r,f[o]});
        });
}
//#####################################################################
// Function Compute_Dof_Pairs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Dof_Pairs(REFERENCE_IRREGULAR_DATA& ri)
{
    auto& ic=irregular_connections(ri.ic_id);
    HASHTABLE<BLOCK_ID,int> prs; // index into d.pairs

    // Copy from irregular neighbors
    REFERENCE_BLOCK_DATA* rd[2];
    rd[0]=&reference_block_data(blocks(ic.regular).ref_id);
    Visit_Irregular_Cross_Section_Dofs(ic,
        blocks(ic.regular).block->cross_sections(ic.con_id),
        [&](int a,BLOCK_ID b,int c,bool o)
        {
            auto y=prs.Insert(b,0);
            if(y.y) *y.x=ri.pairs.Add_End();
            auto& z=ri.pairs(*y.x);
            rd[1]=&reference_block_data(blocks(b).ref_id);
            int f[2]={a,c};
            int i=!o,r=rd[i]->dof_map_v(f[i]),s=rd[i]->dof_map_p(f[i]);
            for(int j=0;j<2;j++)
            {
                if(r>=0) z.irreg_pairs[i][j].v.Append({r,f[j]});
                if(s>=0) z.irreg_pairs[i][j].p.Append({s,f[j]});
            }
        },
        [&](int a,BLOCK_ID b,int c,bool o)
        {
            auto y=prs.Insert(b,0);
            if(y.y) *y.x=ri.pairs.Add_End();
            auto& z=ri.pairs(*y.x);
            rd[1]=&reference_block_data(blocks(b).ref_id);
            int f[2]={a,c};
            int i=!o,r=rd[i]->dof_map_e(f[i]);
            for(int j=0;j<2;j++)
                if(r>=0) z.irreg_pairs[i][j].e.Append({r,f[j]});
        });
}
//#####################################################################
// Function Compute_Dof_Pairs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Compute_Dof_Pairs()
{
}
//#####################################################################
// Function Visualize_Block_State
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Visualize_Block_State(BLOCK_ID b) const
{
    auto& bl=blocks(b);
    auto* cb=bl.block;
    auto Z=[=](int i){return bl.xform*cb->X(i);};
    for(auto t:cb->E)
    {
        VECTOR<TV,3> P;
        for(int i=0;i<3;i++) P(i)=bl.xform*cb->X(t(i));
        for(auto p:P) Add_Debug_Object(VECTOR<TV,2>(p,P.Average()),VECTOR<T,3>(.5,.5,.5));
    }
    HASHTABLE<int> he,hv;
    he.Set_All(cb->bc_e);
    hv.Set_All(cb->bc_v);
    for(int i=0;i<cb->S.m;i++)
    {
        TV X=Z(cb->S(i).x);
        TV Y=Z(cb->S(i).y);
        Add_Debug_Object(VECTOR<TV,2>(X,Y),he.Contains(i)?VECTOR<T,3>(0,1,1):VECTOR<T,3>(0,0,1));
    }
    for(int i=0;i<cb->X.m;i++)
    {
        Add_Debug_Particle(Z(i),hv.Contains(i)?VECTOR<T,3>(1,1,0):VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>("display_size",.1);
    }

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
        {
            const auto* cb2=blocks(c.id).block;
            Visit_Regular_Cross_Section_Dofs(cb->cross_sections(cc),
                cb2->cross_sections(c.con_id),c.master,
                [&](int a,int b,bool o)
                {
                    Add_Debug_Particle(Z(a),o?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                },
                [&](int a,int b,bool o)
                {
                    Add_Debug_Particle((Z(cb->S(a).x)+Z(cb->S(a).y))/2,o?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                });
        }
        else
        {
            const auto& ic=irregular_connections(c.irreg_id);
            Visit_Irregular_Cross_Section_Dofs(ic,
                blocks(ic.regular).block->cross_sections(ic.con_id),
                [&](int a,BLOCK_ID b,int c,bool o)
                {
                    Add_Debug_Particle(Z(a),o?VECTOR<T,3>(1,0,1):VECTOR<T,3>(0,0,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                },
                [&](int a,BLOCK_ID b,int c,bool o)
                {
                    Add_Debug_Particle((Z(cb->S(a).x)+Z(cb->S(a).y))/2,o?VECTOR<T,3>(1,0,1):VECTOR<T,3>(0,0,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                });
        }
    }

    for(IRREG_ID i:bl.edge_on)
    {
        auto& ic=irregular_connections(i);
        Visit_Irregular_Cross_Section_Dofs(ic,
            blocks(ic.regular).block->cross_sections(ic.con_id),
            [&](int a,BLOCK_ID b2,int c,bool o)
            {
                if(b==b2)
                {
                    Add_Debug_Particle(Z(c),o?VECTOR<T,3>(.5,.5,.5):VECTOR<T,3>(1,1,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.05);
                }
            },
            [&](int a,BLOCK_ID b2,int c,bool o)
            {
                if(b==b2)
                {
                    Add_Debug_Particle((Z(cb->S(c).x)+Z(cb->S(c).y))/2,o?VECTOR<T,3>(.5,.5,.5):VECTOR<T,3>(1,1,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.05);
                }
            });
    }
}
//#####################################################################
// Function Visualize_Solution
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Visualize_Solution(BLOCK_ID b) const
{
    const auto& bl=blocks(b);
    const auto* cb=bl.block;
    auto Z=[=](int i){return bl.xform*cb->X(i);};
    const auto& rd=reference_block_data(blocks(b).ref_id);
    const auto& U=rhs_block_list(b);
    for(int i=0;i<cb->X.m;i++)
        if(rd.dof_map_v(i)>=0)
        {
            Add_Debug_Particle(Z(i),VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV>("V",U.Get_v(rd.dof_map_v(i)));
        }

    for(int i=0;i<cb->S.m;i++)
        if(rd.dof_map_e(i)>=0)
        {
            Add_Debug_Particle((Z(cb->S(i).x)+Z(cb->S(i).y))/2,VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV>("V",U.Get_e(rd.dof_map_e(i)));
        }

    for(int i=0;i<cb->X.m;i++)
        if(rd.dof_map_p(i)>=0)
        {
            Add_Debug_Particle(Z(i),VECTOR<T,3>(0,1,0));
            Debug_Particle_Set_Attribute<TV>("display_size",U.Get_p(rd.dof_map_p(i)));
        }
}
//#####################################################################
// Function Transform_Solution
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Transform_Solution(const CACHED_ELIMINATION_MATRIX<T>& cem)
{
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        int j=cem.rhs(Value(b));
        auto& U=rhs_block_list(b);
        U.V=cem.vector_list(j);
        const auto& A = blocks(b).xform.M;
        T s=1/sqrt(A.Determinant());
        auto F = A*s;

        for(int i=0;i<U.nv;i++) U.Set_v(i,F*U.Get_v(i));
        for(int i=0;i<U.ne;i++) U.Set_e(i,F*U.Get_e(i));
        for(int i=0;i<U.np;i++) U.Set_p(i,s*U.Get_p(i));
    }
}
//#####################################################################
// Function Dump_World_Space_System
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Dump_World_Space_System() const
{
    // SYSTEM_MATRIX_HELPER<T> h;

    // for(auto t:nonzero_blocks)
    // {
    //     const BLOCK_MATRIX<T>& B=matrix_block_list(t.z);
    //     BLOCK_MATRIX<T> M(B);
    //     M.M.x.Fill(0);
    //     Transform_To_World_Space(M,B,t.x,t.y);
        
    // }
// struct SYSTEM_MATRIX_HELPER
// {
//     ARRAY<TRIPLE<int,int,T> > data;
//     int start;
//     bool compacted;
}
//#####################################################################
// Function Transform_To_World_Space
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<VECTOR<T,2> >::
Transform_To_World_Space(BLOCK_MATRIX<T>& M,const BLOCK_MATRIX<T>& B,BLOCK_ID a,BLOCK_ID b) const
{
    MATRIX<T,2> Ma=blocks(a).xform.M;
    MATRIX<T,2> Mb=blocks(b).xform.M;
    T sa=1/sqrt(Ma.Determinant());
    T sb=1/sqrt(Mb.Determinant());
    Ma*=sa;
    Mb*=sb;

    for(int p=0;p<B.nv_r;p++)
    {
        for(int r=0;r<B.nv_c;r++) M.Add_vv(p,r,Ma.Transpose_Times(B.Get_vv(p,r)*Mb));
        for(int r=0;r<B.ne_c;r++) M.Add_ve(p,r,Ma.Transpose_Times(B.Get_ve(p,r)*Mb));
        for(int r=0;r<B.np_c;r++) M.Add_vp(p,r,Ma.Transpose_Times(B.Get_vp(p,r)*sb));
    }
    for(int p=0;p<B.ne_r;p++)
    {
        for(int r=0;r<B.nv_c;r++) M.Add_ev(p,r,Ma.Transpose_Times(B.Get_ev(p,r)*Mb));
        for(int r=0;r<B.ne_c;r++) M.Add_ee(p,r,Ma.Transpose_Times(B.Get_ee(p,r)*Mb));
        for(int r=0;r<B.np_c;r++) M.Add_ep(p,r,Ma.Transpose_Times(B.Get_ep(p,r)*sb));
    }
    for(int p=0;p<B.np_r;p++)
    {
        for(int r=0;r<B.nv_c;r++) M.Add_pv(p,r,Mb.Transpose_Times(B.Get_pv(p,r)*sa));
        for(int r=0;r<B.ne_c;r++) M.Add_pe(p,r,Mb.Transpose_Times(B.Get_pe(p,r)*sa));
    }
}
template class COMPONENT_LAYOUT_FEM<VECTOR<double,2> >;
}
