//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "FEM_TABLE.h"
#include "FLAT_SYSTEM_FEM.h"
#include "FLUID_LAYOUT_FEM.h"

namespace PhysBAM{

template<class T>
struct general_fem_coefficients
{
    MATRIX<T,2> viscosity[6][6];
    VECTOR<T,2> pressure[3][6];
};

template<class T, class TV>
void Generate_Hash_FEM_Entries(ARRAY<T,LOCAL_V_CODE_ID>& viscosity_hash,
    ARRAY<T,LOCAL_P_CODE_ID>& pressure_hash,
    const TV& u,const TV& v,T mu)
{
    T S[4] = {u.x, v.x, u.y, v.y}, R[10], scale=mu/(6*u.Cross(v).x);
    for(int i=0,k=0;i<4;i++)
        for(int j=i;j<4;j++)
            R[k++]=scale*S[i]*S[j];

    LOCAL_V_CODE_ID num_hashes_visc=Number_Unique_Visc_Codes();
    viscosity_hash.Resize(num_hashes_visc,no_init);
    for(LOCAL_V_CODE_ID i(0);i<num_hashes_visc;i++)
    {
        T v=0;
        for(int j=0;j<10;j++)
            v+=Unique_Entries_Visc(i,j)*R[j];
        viscosity_hash(i)=v;
    }

    LOCAL_P_CODE_ID num_hashes_pres=Number_Unique_Pres_Codes();
    pressure_hash.Resize(num_hashes_pres,no_init);
    for(LOCAL_P_CODE_ID i(0);i<num_hashes_pres;i++)
    {
        T v=0;
        for(int j=0;j<4;j++)
            v+=Unique_Entries_Pres(i,j)*S[j];
        pressure_hash(i)=v;
    }
}

template<class T>
void Generate_FEM_Entries(general_fem_coefficients<T>& out,
    const ARRAY<T,LOCAL_V_CODE_ID>& viscosity_hash,
    const ARRAY<T,LOCAL_P_CODE_ID>& pressure_hash)
{
    for(int i=0;i<6;i++)
        for(int j=0;j<6;j++)
            for(int k=0;k<2;k++)
                for(int l=0;l<2;l++)
                    out.viscosity[i][j](k,l)=viscosity_hash(Main_Table_Visc(i,k,j,l));

    for(int i=0;i<3;i++)
        for(int j=0;j<6;j++)
            for(int k=0;k<2;k++)
                out.pressure[i][j](k)=pressure_hash(Main_Table_Pres(i,j,k));
}

static int bc_table[3][3]={{4,-1,2},{-1,4,2},{2,2,16}};
template<class TV, class T>
VECTOR<TV,3> Times_BC_NdotN(VECTOR<TV,3>& bc, T edge_length)
{
    VECTOR<TV,3> out;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            out(i)+=bc_table[i][j]*bc(j);
    return edge_length/30*out;
}

template<class T,class TV>
VECTOR<TV,6> Times_force_NdotN(VECTOR<TV,6>& f, T tri_area)
{
    VECTOR<TV,6> out;
    for(int i=0;i<6;i++)
        for(int j=0;j<6;j++)
            out(i)+=force_table[i][j]*f(j);
    return tri_area/360*out;
}

static const auto num_v_codes=Number_Unique_Visc_Codes();
static const auto num_p_codes=Number_Unique_Pres_Codes();

inline CODE_ID Code(PIPE_ID p,LOCAL_V_CODE_ID c)
{
    return CODE_ID(Value(p)*(Value(num_v_codes)+Value(num_p_codes))+Value(num_p_codes)+Value(c));
};
inline CODE_ID Code(PIPE_ID p,LOCAL_P_CODE_ID c)
{
    return CODE_ID(Value(p)*(Value(num_v_codes)+Value(num_p_codes))+Value(c));
};
inline LOCAL_P_CODE_ID Neg(LOCAL_P_CODE_ID c)
{
    return LOCAL_P_CODE_ID(1^Value(c));
}
inline LOCAL_P_CODE_ID Cond_Neg(LOCAL_P_CODE_ID c,int neg)
{
    return LOCAL_P_CODE_ID(neg^Value(c));
}

int Missing_Vertex(VECTOR<PARTICLE_ID,3> t, VECTOR<PARTICLE_ID,2> e)
{
    for(int i=0;i<3;i++)
        if(!e.Contains(t(i)))
            return i;
    PHYSBAM_FATAL_ERROR();
}

template<class T,class TV>
void Generate_Discretization(ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values,const FLUID_LAYOUT_FEM<TV>& fl,
    const PARSE_DATA_FEM<TV>& pd,T mu)
{

    PIPE_ID unknown_pipe(-2);
    PIPE_ID irregular_pipe(-1);

    
    ARRAY<general_fem_coefficients<T>,TRIANGLE_ID> irregular_data(fl.Number_Triangles());
    ARRAY<general_fem_coefficients<T>,PIPE_ID> pipe_data(fl.pipes.m);
    
    for(PIPE_ID i(0);i<fl.pipes.m;i++)
    {
        ARRAY<T,LOCAL_V_CODE_ID> viscosity_hash;
        ARRAY<T,LOCAL_P_CODE_ID> pressure_hash;
        Generate_Hash_FEM_Entries(viscosity_hash,pressure_hash,fl.pipes(i).u,fl.pipes(i).v,mu);
        Generate_FEM_Entries(pipe_data(i),viscosity_hash,pressure_hash);
        for(auto t:pressure_hash) code_values.Append(t);
        for(auto t:viscosity_hash) code_values.Append(t);
    }
    for(TRIANGLE_ID i(0);i<irregular_data.m;i++)
    {
        if(fl.blocks(fl.elem_data(i).block_id).pipe_id>=PIPE_ID()) continue;
        ARRAY<T,LOCAL_V_CODE_ID> viscosity_hash;
        ARRAY<T,LOCAL_P_CODE_ID> pressure_hash;
        auto tri=fl.Triangle(i);
        VECTOR<TV,3> X(fl.X(tri.x),fl.X(tri.y),fl.X(tri.z));
        Generate_Hash_FEM_Entries(viscosity_hash,pressure_hash,X(1)-X(0),X(2)-X(0),mu);
        Generate_FEM_Entries(irregular_data(i),viscosity_hash,pressure_hash);
    }

    for(PARTICLE_ID p(0);p<fl.Number_Particles();p++)
    {
        PIPE_ID regular=unknown_pipe;
        auto elems=fl.Incident_Triangles(p);
        for(auto t:elems)
        {
            PIPE_ID pipe=fl.blocks(fl.elem_data(t).block_id).pipe_id;
            if(pipe<PIPE_ID()) regular=irregular_pipe;
            else if(regular==unknown_pipe) regular=pipe;
            else if(regular!=pipe) regular=irregular_pipe;
        }

        DOF_ID dof_p=fl.pressure_dofs(p);
        DOF_ID dof_u[2]={fl.vel_node_dofs(p),fl.vel_node_dofs(p)+1};

        if(regular<PIPE_ID())
        {
            MATRIX<T,TV::m> visc;
            TV pres;
            for(auto t:elems)
            {
                PIPE_ID pipe=fl.blocks(fl.elem_data(t).block_id).pipe_id;
                int k=fl.Triangle(t).Find(p);
                const general_fem_coefficients<T>& fem=pipe<PIPE_ID()?irregular_data(t):pipe_data(pipe);
                visc+=fem.viscosity[k][k];
                int sign=(pipe>=PIPE_ID() && (t-fl.pipes(regular).first_element)%2)?-1:1;
                pres+=fem.pressure[k][k]*sign;
            }
            CODE_ID codes_u[3];
            CODE_ID codes_p[2];
            for(int i=0;i<2;i++)
                for(int j=i;j<2;j++)
                    codes_u[i+j]=code_values.Append(visc(i,j));
            for(int i=0;i<2;i++)
                codes_p[i]=code_values.Append(-pres(i));

            for(int i=0;i<2;i++)
                for(int j=0;j<2;j++)
                    coded_entries.Append({dof_u[i],dof_u[j],codes_u[i+j]});

            for(int i=0;i<2;i++){
                coded_entries.Append({dof_p,dof_u[i],codes_p[i]});
                coded_entries.Append({dof_u[i],dof_p,codes_p[i]});}
        }
        else
        {
            LOCAL_V_CODE_ID visc[2][2]={};
            LOCAL_P_CODE_ID pres[2]={};
            for(int i=0;i<2;i++)
                for(int j=0;j<2;j++)
                    visc[i][j]=Vertex_Table_Visc(i,j,elems.m);

            if(elems.m<6)
            {
                TRIANGLE_ID fe=fl.pipes(regular).first_element;
                int x=0,z=0;
                for(auto t:elems)
                {
                    int k=fl.Triangle(t).Find(p);
                    int n=(t-fe)&1;
                    if(k==0) z+=n+1;
                    if(k==1) x+=n+1;
                }
                for(int i=0;i<2;i++)
                    pres[i]=Vertex_Table_Pres(i,x,z);
            }

            for(int i=0;i<2;i++)
                for(int j=0;j<2;j++)
                    coded_entries.Append({dof_u[i],dof_u[j],Code(regular,visc[i][j])});

            for(int i=0;i<2;i++){
                CODE_ID code=Code(regular,Neg(pres[i]));
                coded_entries.Append({dof_p,dof_u[i],code});
                coded_entries.Append({dof_u[i],dof_p,code});}
        }
    }

    for(EDGE_ID e(0);e<fl.Number_Edges();e++)
    {
        PIPE_ID regular=unknown_pipe;
        auto tris=fl.Edge_Triangles(e);
        for(auto t:tris)
        {
            PIPE_ID pipe=fl.blocks(fl.elem_data(t).block_id).pipe_id;
            if(pipe<PIPE_ID()) regular=irregular_pipe;
            else if(regular==unknown_pipe) regular=pipe;
            else if(regular!=pipe) regular=irregular_pipe;
        }


        TRIANGLE_ID t0=tris(0);
        auto tri0=fl.Triangle(t0);
        int c0=Missing_Vertex(tri0,fl.Edge(e));
        PARTICLE_ID v0[2]={tri0((c0+1)%3),tri0((c0+2)%3)};
        DOF_ID dof_p[2];
        DOF_ID dof_u[3][2];
        for(int i=0;i<2;i++)
        {
            dof_p[i]=fl.pressure_dofs(v0[i]);
            dof_u[i][0]=fl.vel_node_dofs(v0[i]);
            dof_u[i][1]=fl.vel_node_dofs(v0[i])+1;
        }
        dof_u[2][0]=fl.vel_edge_dofs(e);
        dof_u[2][1]=fl.vel_edge_dofs(e)+1;

        if(regular<PIPE_ID())
        {
            MATRIX<T,TV::m> visc[3][3];
            TV pres[2][3];
            for(int i=0;i<tris.m;i++)
            {
                TRIANGLE_ID t=tris(i);
                PIPE_ID pipe=fl.blocks(fl.elem_data(t).block_id).pipe_id;
                const general_fem_coefficients<T>& fem=pipe<PIPE_ID()?irregular_data(t):pipe_data(pipe);

                auto tri=fl.Triangle(t);
                int c=Missing_Vertex(tri,fl.Edge(e));
                int di[3]={(c+1+i)%3,(c+2-i)%3,c+3};
                int sign=(pipe>=PIPE_ID() && (t-fl.pipes(regular).first_element)%2)?-1:1;
                for(int m=0;m<3;m++)
                    for(int n=0;n<3;n++)
                    {
                        if(m<2 && m==n) continue;
                        visc[m][n]+=fem.viscosity[di[m]][di[n]];
                    }
                for(int m=0;m<2;m++)
                    for(int n=0;n<3;n++)
                        pres[m][n]+=fem.pressure[di[m]][di[n]]*sign;
            }

            for(int m=0;m<3;m++)
                for(int n=0;n<3;n++)
                    for(int a=0;a<2;a++)
                        for(int b=0;b<2;b++)
                        {
                            if(m<2 && m==n) continue;
                            DOF_ID d0=dof_u[m][a];
                            DOF_ID d1=dof_u[n][b];
                            if(d1<d0) continue;
                            CODE_ID code=code_values.Append(visc[m][n](a,b));
                            coded_entries.Append({d0,d1,code});
                            coded_entries.Append({d1,d0,code});
                        }

            for(int m=0;m<2;m++)
                for(int n=0;n<3;n++)
                    for(int b=0;b<2;b++)
                    {
                        DOF_ID d0=dof_p[m];
                        DOF_ID d1=dof_u[n][b];
                        CODE_ID code=code_values.Append(-pres[m][n](b));
                        coded_entries.Append({d0,d1,code});
                        coded_entries.Append({d1,d0,code});
                    }

        }
        else
        {
            TRIANGLE_ID t0=tris(0);
            auto tri0=fl.Triangle(t0);
            int c0=Missing_Vertex(tri0,fl.Edge(e));
            int di[3]={(c0+1)%3,(c0+2)%3,c0+3};
            int sign=(regular>=PIPE_ID() && Value(t0-fl.pipes(regular).first_element)%2);
            LOCAL_V_CODE_ID visc[3][3][2][2];
            LOCAL_P_CODE_ID pres[2][3][2];
            memset(visc,-1,sizeof(visc));
            memset(pres,-1,sizeof(pres));

            if(tris.m==1)
            {
                for(int m=0;m<3;m++)
                    for(int n=0;n<3;n++)
                        for(int a=0;a<2;a++)
                            for(int b=0;b<2;b++)
                            {
                                if(m<2 && m==n) continue;
                                visc[m][n][a][b]=Main_Table_Visc(di[m],a,di[n],b);
                            }

                for(int m=0;m<2;m++)
                    for(int n=0;n<3;n++)
                        for(int b=0;b<2;b++)
                            pres[m][n][b]=Cond_Neg(Main_Table_Pres(di[m],di[n],b),sign);
            }
            else
            {
                assert(tris.m==2);
                TRIANGLE_ID t1=tris(1);
                auto tri1=fl.Triangle(t1);
                int c1=Missing_Vertex(tri1,fl.Edge(e));

                for(int m=0;m<3;m++)
                    for(int n=0;n<3;n++)
                        for(int a=0;a<2;a++)
                            for(int b=0;b<2;b++)
                            {
                                if(m<2 && m==n) continue;
                                visc[m][n][a][b]=Edge_Table_Visc(c0,c1,m,a,n,b);
                            }

                for(int m=0;m<2;m++)
                    for(int n=0;n<3;n++)
                        for(int b=0;b<2;b++)
                            pres[m][n][b]=Cond_Neg(Edge_Table_Pres(c0,c1,m,n,b),sign);
            }

            for(int m=0;m<3;m++)
                for(int n=0;n<3;n++)
                    for(int a=0;a<2;a++)
                        for(int b=0;b<2;b++)
                        {
                            if(m<2 && m==n) continue;
                            DOF_ID d0=dof_u[m][a];
                            DOF_ID d1=dof_u[n][b];
                            if(d1<d0) continue;
                            CODE_ID code=Code(regular,visc[m][n][a][b]);
                            coded_entries.Append({d0,d1,code});
                            coded_entries.Append({d1,d0,code});
                        }

            for(int m=0;m<2;m++)
                for(int n=0;n<3;n++)
                    for(int b=0;b<2;b++)
                    {
                        DOF_ID d0=dof_p[m];
                        DOF_ID d1=dof_u[n][b];
                        CODE_ID code=Code(regular,Neg(pres[m][n][b]));
                        coded_entries.Append({d0,d1,code});
                        coded_entries.Append({d1,d0,code});
                    }
        }
    }

    for(TRIANGLE_ID t(0);t<fl.Number_Triangles();t++)
    {
        PIPE_ID regular=fl.blocks(fl.elem_data(t).block_id).pipe_id;
        auto tri=fl.Triangle(t);

        DOF_ID dof_p[3];
        DOF_ID dof_u[6][2];
        for(int i=0;i<3;i++)
        {
            dof_p[i]=fl.pressure_dofs(tri(i));
            dof_u[i][0]=fl.vel_node_dofs(tri(i));
            dof_u[i][1]=fl.vel_node_dofs(tri(i))+1;
        }
        for(auto e:fl.Triangle_Edges(t))
        {
            int c=Missing_Vertex(tri,fl.Edge(e));
            dof_u[c+3][0]=fl.vel_edge_dofs(e);
            dof_u[c+3][1]=fl.vel_edge_dofs(e)+1;
        }

        if(regular<PIPE_ID())
        {
            MATRIX<T,TV::m> visc_ee[3]; // i+3 -> (i+1)%3+3
            MATRIX<T,TV::m> visc_ve[3]; // i -> i+3
            TV pres_ve[3]; // i -> i+3

            const general_fem_coefficients<T>& fem=irregular_data(t);

            for(int m=0;m<3;m++)
            {
                visc_ee[m]=fem.viscosity[m+3][(m+1)%3+3];
                visc_ve[m]=fem.viscosity[m][m+3];
            }
            for(int m=0;m<2;m++)
                pres_ve[m]=fem.pressure[m][m+3];

            for(int m=0;m<3;m++)
                for(int a=0;a<2;a++)
                    for(int b=0;b<2;b++)
                    {
                        DOF_ID d0=dof_u[m+3][a];
                        DOF_ID d1=dof_u[(m+1)%3+3][b];
                        CODE_ID code=code_values.Append(visc_ee[m](a,b));
                        coded_entries.Append({d0,d1,code});
                        coded_entries.Append({d1,d0,code});
                        d0=dof_u[m][a];
                        d1=dof_u[m+3][b];
                        code=code_values.Append(visc_ve[m](a,b));
                        coded_entries.Append({d0,d1,code});
                        coded_entries.Append({d1,d0,code});
                    }

            for(int m=0;m<3;m++)
                for(int b=0;b<2;b++)
                {
                    DOF_ID d0=dof_p[m];
                    DOF_ID d1=dof_u[m+3][b];
                    CODE_ID code=code_values.Append(-pres_ve[m](b));
                    coded_entries.Append({d0,d1,code});
                    coded_entries.Append({d1,d0,code});
                }
        }
        else
        {
            int sign=(t-fl.pipes(regular).first_element)%2;

            for(int m=0;m<3;m++)
                for(int a=0;a<2;a++)
                    for(int b=0;b<2;b++)
                    {
                        LOCAL_V_CODE_ID visc_ee=Main_Table_Visc(m+3,a,(m+1)%3+3,b);
                        LOCAL_V_CODE_ID visc_ve=Main_Table_Visc(m,a,m+3,b);
                        DOF_ID d0=dof_u[m+3][a];
                        DOF_ID d1=dof_u[(m+1)%3+3][b];
                        CODE_ID code=Code(regular,visc_ee);
                        coded_entries.Append({d0,d1,code});
                        coded_entries.Append({d1,d0,code});
                        d0=dof_u[m][a];
                        d1=dof_u[m+3][b];
                        code=Code(regular,visc_ve);
                        coded_entries.Append({d0,d1,code});
                        coded_entries.Append({d1,d0,code});
                    }

            for(int m=0;m<2;m++)
                for(int b=0;b<2;b++)
                {
                    LOCAL_P_CODE_ID pres_ve=Cond_Neg(Main_Table_Pres(m,m+3,b),sign);
                    DOF_ID d0=dof_p[m];
                    DOF_ID d1=dof_u[m+3][b];
                    CODE_ID code=Code(regular,Neg(pres_ve));
                    coded_entries.Append({d0,d1,code});
                    coded_entries.Append({d1,d0,code});
                }
        }
    }
}
template void Generate_Discretization<double,VECTOR<double,2> >(ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID>,int>&,ARRAY<double,CODE_ID>&,FLUID_LAYOUT_FEM<VECTOR<double,2> > const&,PARSE_DATA_FEM<VECTOR<double,2> > const&,double);

}
