//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <fstream>
#include "LAYOUT_BUILDER_FEM.h"
namespace PhysBAM{

template<class T>
bool Canonical_Direction(VECTOR<T,2> u)
{
    double comp_tol=1e-10;
    if(u.x>comp_tol) return true;
    if(u.x<-comp_tol) return false;
    if(u.x>comp_tol) return true;
    return false;
}

//#####################################################################
// Function Constructor
//#####################################################################
template<class T> LAYOUT_BUILDER_FEM<T>::
LAYOUT_BUILDER_FEM(COMPONENT_LAYOUT_FEM<T>& cl):
    cl(cl)
{
}
//#####################################################################
// Function From_File
//#####################################################################
template<class T> void LAYOUT_BUILDER_FEM<T>::
From_File(const std::string& file)
{
    std::ifstream fin(file);
    std::string line;

    std::string name,name2,name3,name4;
    int i0;
    T t0,t1;
    char c;
    TV v0;
    HASHTABLE<std::string,CS_ID> cs_names;
    HASHTABLE<std::string,VERT_ID> vert_names;
    HASHTABLE<std::string,CONNECTOR_ID> con_names;
    
    while(getline(fin,line))
    {
        std::stringstream ss(line);
        ss>>c;
        if(isspace(c) || c=='#') continue;
        switch(c)
        {
            case 'l':
                ss>>t0;
                Set_Target_Length(t0);
                break;

            case 'c':
                {
                    ss>>name>>i0>>t0;
                    CS_ID cs=Cross_Section(i0,t0);
                    cs_names.Set(name,cs);
                }
                break;

            case 'v':
                {
                    ss>>name>>v0;
                    VERT_ID vid=Vertex(v0);
                    vert_names.Set(name,vid);
                }
                break;

            case 'j':
                {
                    ss>>name2>>i0>>name3;
                    CS_ID cs=cs_names.Get(name2);
                    VERT_ID o=vert_names.Get(name3);
                    ARRAY<VERT_ID> arms(i0);
                    ARRAY<std::string> c(i0);
                    for(int i=0;i<i0;i++)
                    {
                        ss>>name>>c(i);
                        arms(i)=vert_names.Get(name);
                    }
                    ARRAY<CONNECTOR_ID> r=Joint(cs,i0,o,arms);
                    for(int i=0;i<i0;i++)
                        con_names.Set(c(i),r(i));
                }
                break;

            case 'p':
                {
                    ss>>name>>name2>>name3;
                    Pipe(cs_names.Get(name),con_names.Get(name2),con_names.Get(name3));
                }
                break;

            case 'g':
                {
                    ss>>name>>name2;
                    CS_ID cs0=cs_names.Get(name);
                    CS_ID cs1=cs_names.Get(name2);
                    ss>>name>>name2>>t0>>t1;
                    VERT_ID v0=vert_names.Get(name);
                    VERT_ID v1=vert_names.Get(name2);
                    ss>>name>>name2;
                    VECTOR<CONNECTOR_ID,2> r=Pipe(cs0,cs1,v0,v1,t0,t1);
                    con_names.Set(name,r(0));
                    con_names.Set(name2,r(1));
                }
                break;
            case 'u':
                {
                    ss>>name>>name2>>name3>>name4>>t0;
                    auto c=Set_BC(cs_names.Get(name),vert_names.Get(name2),vert_names.Get(name3),t0);
                    con_names.Set(name4,c.x);
                }
                break;
            case 't':
                {
                    ss>>name>>name2>>name3>>name4>>v0;
                    auto c=Set_BC(cs_names.Get(name),vert_names.Get(name2),vert_names.Get(name3),v0);
                    con_names.Set(name4,c.x);
                }
                break;
            default:
                LOG::printf("PARSE FAIL: %c %s\n",c,ss.str());
        }
    }
}
//#####################################################################
// Function To_String
//#####################################################################
template<class T> std::string LAYOUT_BUILDER_FEM<T>::
To_String() const
{
    std::ostringstream os;
    for(int i=0;i<commands.m;i++)
    {
        const auto& c=commands(i);
        switch(c.x)
        {
            case TAG::SET_LEN:
                os<<"l "<<cl.target_length;
                break;
            case TAG::DECL_CS:
                {
                    CS_ID cs(c.y);
                    os<<"c C"<<c.y<<" "<<cross_sections(cs).x<<" "<<cross_sections(cs).y;
                }
                break;
            case TAG::DECL_VERT:
                os<<"v "<<"V"<<c.y<<" "<<verts(VERT_ID(c.y))(0)<<" "<<verts(VERT_ID(c.y))(1);
                break;
            case TAG::SET_BC_U:
                {
                    PHYSBAM_ASSERT(commands(i+1).x==TAG::CS);
                    PHYSBAM_ASSERT(commands(i+2).x==TAG::VERT);
                    PHYSBAM_ASSERT(commands(i+3).x==TAG::VERT);
                    PHYSBAM_ASSERT(commands(i+4).x==TAG::CONNECTOR);
                    CS_ID cs(commands(i+1).y);
                    VERT_ID v0(commands(i+2).y),v1(commands(i+3).y);
                    os<<"u "<<"C"<<cs<<" V"<<v0<<" V"<<v1<<" K"<<commands(i+4).y<<" "<<cl.bc_v(c.y).flow_rate;
                    i+=4;
                }
                break;
            case TAG::SET_BC_T:
                {
                    PHYSBAM_ASSERT(commands(i+1).x==TAG::CS);
                    PHYSBAM_ASSERT(commands(i+2).x==TAG::VERT);
                    PHYSBAM_ASSERT(commands(i+3).x==TAG::VERT);
                    PHYSBAM_ASSERT(commands(i+4).x==TAG::CONNECTOR);
                    CS_ID cs(commands(i+1).y);
                    VERT_ID v0(commands(i+2).y),v1(commands(i+3).y);
                    os<<"t C"<<cs<<" V"<<v0<<" V"<<v1<<" K"<<commands(i+4).y;
                    os<<" "<<cl.bc_t(c.y).traction(0)<<" "<<cl.bc_t(c.y).traction(1);
                    i+=4;
                }
                break;
            case TAG::DECL_JT:
                {
                    PHYSBAM_ASSERT(commands(i+1).x==TAG::CS);
                    PHYSBAM_ASSERT(commands(i+2).x==TAG::VERT);
                    int n=commands(i).y;
                    os<<"j C"<<commands(i+1).y<<" "<<n<<" V"<<commands(i+2).y;
                    for(int j=0;j<n;j++)
                    {
                        PHYSBAM_ASSERT(commands(i+3+j*2).x==TAG::VERT);
                        PHYSBAM_ASSERT(commands(i+3+j*2+1).x==TAG::CONNECTOR);
                        os<<" V"<<commands(i+3+j*2).y<<" K"<<commands(i+3+j*2+1).y;
                    }
                    i+=2+2*n;
                }
                break;
            case TAG::DECL_PIPE:
                {
                    if(commands(i).y==-1)
                    {
                        PHYSBAM_ASSERT(commands(i+1).x==TAG::CS);
                        PHYSBAM_ASSERT(commands(i+2).x==TAG::CONNECTOR);
                        PHYSBAM_ASSERT(commands(i+3).x==TAG::CONNECTOR);
                        os<<"p C"<<commands(i+1).y<<" K"<<commands(i+2).y<<" K"<<commands(i+3).y;
                        i+=3;
                    }
                    else
                    {
                        PHYSBAM_ASSERT(commands(i+1).x==TAG::CS);
                        PHYSBAM_ASSERT(commands(i+2).x==TAG::CS);
                        PHYSBAM_ASSERT(commands(i+3).x==TAG::VERT);
                        PHYSBAM_ASSERT(commands(i+4).x==TAG::VERT);
                        PHYSBAM_ASSERT(commands(i+5).x==TAG::CONNECTOR);
                        PHYSBAM_ASSERT(commands(i+6).x==TAG::CONNECTOR);
                        int g=commands(i).y;
                        os<<"g C"<<commands(i+1).y<<" C"<<commands(i+2).y;
                        os<<" V"<<commands(i+3).y<<" V"<<commands(i+4).y;
                        os<<" "<<var_size_pipes(g).x<<" "<<var_size_pipes(g).y;
                        os<<" K"<<commands(i+5).y<<" K"<<commands(i+6).y;
                        i+=6;
                    }
                }
                break;
            default:
                LOG::printf("Failed to serialize at %d\n",i);
        }
        if(i!=commands.m-1) os<<"\n";
    }
    return os.str();
}
//#####################################################################
// Function Set_Target_Length
//#####################################################################
template<class T> void LAYOUT_BUILDER_FEM<T>::
Set_Target_Length(T l)
{
    cl.target_length=l*cl.unit_m;
    comp_pipe.target_length=cl.target_length;
    comp_change.target_length=cl.target_length;
    comp_bc.target_length=cl.target_length;
    comp_joint.target_length=cl.target_length;
    commands.Append({TAG::SET_LEN,-1});
}
//#####################################################################
// Function Cross_Section
//#####################################################################
template<class T> auto LAYOUT_BUILDER_FEM<T>::
Cross_Section(int d,T w) -> CS_ID
{
    w*=cl.unit_m;
    CS_ID c=cross_sections.Append({d,w});
    commands.Append({TAG::DECL_CS,Value(c)});
    return c;
}
//#####################################################################
// Function Vertex
//#####################################################################
template<class T> auto LAYOUT_BUILDER_FEM<T>::
Vertex(const TV& X) -> VERT_ID
{
    TV Y=X*cl.unit_m;
    VERT_ID v=verts.Append(Y);
    commands.Append({TAG::DECL_VERT,Value(v)});
    return v;
}
//#####################################################################
// Function Set_BC
//#####################################################################
template<class T> auto LAYOUT_BUILDER_FEM<T>::
Set_BC(CS_ID cs,VERT_ID from,VERT_ID to,T flow_rate) -> PAIR<CONNECTOR_ID,BC_U_ID>
{
    T len=cl.target_length;
    auto crs=cross_sections(cs);
    flow_rate*=cl.unit_m/cl.unit_s;
    TV A=verts(from),B=verts(to);
    TV dir=(B-A).Normalized();
    TV C=A+dir*len;
    auto pr=comp_bc.Make_Block(crs.x,crs.y,len,true);
    BLOCK_ID b=cl.blocks.Append({pr.x,{Compute_Xform(dir),A},{BLOCK_ID(-1)}});
    VERTEX_DATA vd={C,{b,CON_ID()}};
    connectors.Set(last_cid,vd);
    BOUNDARY_CONDITION<T> bc={b,pr.y,pr.z,flow_rate,{},-dir};
    int bc_idx=cl.bc_v.Append(bc);
    ARRAY<PAIR<TAG,int> > com=
    {
        {TAG::SET_BC_U,bc_idx},{TAG::CS,Value(cs)},
        {TAG::VERT,Value(from)},{TAG::VERT,Value(to)},{TAG::CONNECTOR,Value(last_cid)}
    };
    commands.Append_Elements(com);
    return {last_cid++,BC_U_ID(bc_idx)};
}
//#####################################################################
// Function Set_BC
//#####################################################################
template<class T> auto LAYOUT_BUILDER_FEM<T>::
Set_BC(CS_ID cs,VERT_ID from,VERT_ID to,const TV& traction) -> PAIR<CONNECTOR_ID,BC_T_ID>
{
    T len=cl.target_length;
    auto crs=cross_sections(cs);
    TV t=traction*cl.unit_kg/sqr(cl.unit_s);
    TV A=verts(from),B=verts(to);
    TV dir=(B-A).Normalized();
    TV C=A+dir*len;
    auto pr=comp_bc.Make_Block(crs.x,crs.y,len,false);
    BLOCK_ID b=cl.blocks.Append({pr.x,{Compute_Xform(dir),A},{BLOCK_ID(-1)}});
    VERTEX_DATA vd={C,{b,CON_ID()}};
    connectors.Set(last_cid,vd);
    BOUNDARY_CONDITION<T> bc={b,pr.y,pr.z,0,t,-dir};
    int bc_idx=cl.bc_t.Append(bc);
    ARRAY<PAIR<TAG,int> > com=
    {
        {TAG::SET_BC_T,bc_idx},{TAG::CS,Value(cs)},
        {TAG::VERT,Value(from)},{TAG::VERT,Value(to)},{TAG::CONNECTOR,Value(last_cid)}
    };
    commands.Append_Elements(com);
    return {last_cid++,BC_T_ID(bc_idx)};
}
//#####################################################################
// Function Joint
//#####################################################################
template<class T> auto LAYOUT_BUILDER_FEM<T>::
Joint(CS_ID cs,int n,VERT_ID o,const ARRAY<VERT_ID>& arms) -> ARRAY<CONNECTOR_ID>
{
    auto crs=cross_sections(cs);
    TV O=verts(o);
    ARRAY<TV> vertx;
    ARRAY<PAIR<TV,CONNECTOR_ID> > vc;
    for(VERT_ID p:arms)
    {
        CONNECTOR_ID c=last_cid++;
        vc.Append({(verts(p)-O).Normalized(),c});
    }
    TV first=vc(0).x;
    auto comp=[&first](const auto& a,const auto& b)
    {
        T xa=TV::Oriented_Angle_Between(first,a.x);
        if(xa<0) xa+=2*pi;
        T xb=TV::Oriented_Angle_Between(first,b.x);
        if(xb<0) xb+=2*pi;
        return xa<xb;
    };
    std::sort(vc.begin()+1,vc.end(),comp);
    ARRAY<T> angles;
    for(int i=1;i<n;i++)
    {
        T a=TV::Oriented_Angle_Between(vc(i-1).x,vc(i).x);
        if(a<0) a+=2*pi;
        angles.Append(a);
    }
    auto cj=comp_joint.Make_Component(crs.x,crs.y,angles);
    XFORM<TV> xf={Compute_Xform(vc(0).x),O};
    ARRAY<VERTEX_DATA> vd(n);
    Emit_Component_Blocks(cj.x,xf,vd);
    ARRAY<CONNECTOR_ID> r;
    for(int i=0;i<n;i++)
    {
        vd(i).X=O+vc(i).x*cj.y(i);
        connectors.Set(vc(i).y,vd(i));
        r.Append(vc(i).y);
    }
    ARRAY<PAIR<TAG,int> > com={{TAG::DECL_JT,n},{TAG::CS,Value(cs)},{TAG::VERT,Value(o)}};
    commands.Append_Elements(com);
    for(int i=0;i<n;i++)
    {
        commands.Append({TAG::VERT,Value(arms(i))});
        commands.Append({TAG::CONNECTOR,Value(vc(i).y)});
    }
    return r;
}
//#####################################################################
// Function Pipe
//#####################################################################
template<class T> void LAYOUT_BUILDER_FEM<T>::
Pipe(CS_ID cs,CONNECTOR_ID a,CONNECTOR_ID b)
{
    ARRAY<VERTEX_DATA> vd(2);
    vd(0)=connectors.Get(a);
    vd(1)=connectors.Get(b);
    connectors.Delete(a);
    connectors.Delete(b);
    auto crs=cross_sections(cs);
    TV dir=vd(1).X-vd(0).X;
    if(!Canonical_Direction(dir))
    {
        dir=-dir;
        std::swap(vd(0),vd(1));
    }
    T len=dir.Normalize();
    XFORM<TV> xf={Compute_Xform(dir),vd(0).X};
    auto cc=comp_pipe.Make_Component(crs.x,crs.y,len);
    Emit_Component_Blocks(cc,xf,vd);
    ARRAY<PAIR<TAG,int> > com={{TAG::DECL_PIPE,-1},{TAG::CS,Value(cs)},{TAG::CONNECTOR,Value(a)},{TAG::CONNECTOR,Value(b)}};
    commands.Append_Elements(com);
}
//#####################################################################
// Function Pipe
//#####################################################################
template<class T> auto LAYOUT_BUILDER_FEM<T>::
Pipe(CS_ID cs0,CS_ID cs1,VERT_ID v0,VERT_ID v1,T offset,T length) -> VECTOR<CONNECTOR_ID,2>
{
    auto crs0=cross_sections(cs0);
    auto crs1=cross_sections(cs1);
    offset*=cl.unit_m;
    length*=cl.unit_m;
    TV A=verts(v0),B=verts(v1);
    TV dir=(B-A).Normalized();
    if(!Canonical_Direction(dir))
    {
        dir=-dir;
        std::swap(A,B);
    }
    TV C=A+dir*offset;
    TV D=C+dir*length;
    auto cc=comp_change.Make_Component(crs0.x,crs0.y,crs1.x,crs1.y,length);
    XFORM<TV> xf={Compute_Xform(dir),C};

    ARRAY<VERTEX_DATA> vd(2);
    Emit_Component_Blocks(cc,xf,vd);
    vd(0).X=C;
    vd(1).X=D;
    VECTOR<CONNECTOR_ID,2> r(last_cid,last_cid+1);
    last_cid+=2;
    connectors.Set(r(0),vd(0));
    connectors.Set(r(1),vd(1));
    int g=var_size_pipes.Append({offset,length});
    ARRAY<PAIR<TAG,int> > com=
    {
        {TAG::DECL_PIPE,g},
        {TAG::CS,Value(cs0)},{TAG::CS,Value(cs1)},
        {TAG::VERT,Value(v0)},{TAG::VERT,Value(v1)},
        {TAG::CONNECTOR,Value(r(0))},{TAG::CONNECTOR,Value(r(1))}
    };
    commands.Append_Elements(com);
    return r;
}
//#####################################################################
// Function Set_Connector
//#####################################################################
template<class T> void LAYOUT_BUILDER_FEM<T>::
Set_Connector(VERTEX_DATA& vd,BLOCK_ID id,CON_ID con_id)
{
    if(vd.con.is_regular)
    {
        if(vd.con.id>=BLOCK_ID())
        {
            cl.blocks(vd.con.id).connections(vd.con.con_id)={id,con_id};
            cl.blocks(id).connections(con_id)={vd.con.id,vd.con.con_id};
        }
        else
        {
            vd.con.id=id;
            vd.con.con_id=con_id;
        }
    }
    else
    {
        IRREGULAR_CONNECTION& con=cl.irregular_connections(vd.con.irreg_id);
        con.regular=id;
        con.con_id=con_id;
        cl.blocks(id).connections(con_id).Set_Irreg(vd.con.irreg_id);
    }
}
//#####################################################################
// Function Emit_Component_Blocks
//#####################################################################
template<class T> void LAYOUT_BUILDER_FEM<T>::
Emit_Component_Blocks(const CANONICAL_COMPONENT<T>* cc,const XFORM<TV>& xf,ARRAY<VERTEX_DATA>& vd)
{
    int offset=Value(cl.blocks.m),offset_edge_on=Value(cl.irregular_connections.m);
    auto mp_bl=[=](CC_BLOCK_ID b){return BLOCK_ID(Value(b)+offset);};
    auto mp_ic=[=](CC_IRREG_ID i){return IRREG_ID(Value(i)+offset_edge_on);};

    for(CC_BLOCK_ID b(0);b<cc->blocks.m;b++)
    {
        const auto& cbl=cc->blocks(b);
        auto& bl=cl.blocks(cl.blocks.Add_End());
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
            bl.edge_on.Append({mp_ic(e.x),e.y});
    }

    for(const auto& a:cc->irregular_connections)
    {
        IRREG_ID index=cl.irregular_connections.Add_End();
        IRREGULAR_CONNECTION& ic=cl.irregular_connections(index);
        for(auto& b:a.edge_on) ic.edge_on.Append({mp_bl(b.b),b.e,b.v0,b.v1});

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
template<class T> auto LAYOUT_BUILDER_FEM<T>::
Compute_Xform(const TV& dir) -> MATRIX<T,2>
{
    return MATRIX<T,TV::m>(dir,dir.Orthogonal_Vector());
}
template class LAYOUT_BUILDER_FEM<double>;
}
