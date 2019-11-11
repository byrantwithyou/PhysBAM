//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
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
    assert(scale>=0);
    for(int i=0,k=0;i<4;i++)
        for(int j=i;j<4;j++)
            R[k++]=scale*S[i]*S[j];
    for(auto& s:S) s/=6;

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
        {
            auto M=Main_Table_Visc(i,j);
            for(int k=0;k<2;k++)
                for(int l=0;l<2;l++)
                    out.viscosity[i][j](k,l)=viscosity_hash(M(k,l));
        }

    for(int i=0;i<3;i++)
        for(int j=0;j<6;j++)
        {
            auto V=Main_Table_Pres(i,j);
            for(int k=0;k<2;k++)
                out.pressure[i][j](k)=pressure_hash(V(k));
        }
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
inline VECTOR<LOCAL_P_CODE_ID,2> Neg(VECTOR<LOCAL_P_CODE_ID,2> c)
{
    for(auto&i:c) if(i) i=LOCAL_P_CODE_ID(Value(i)^1);
    return c;
}
inline VECTOR<LOCAL_P_CODE_ID,2> Cond_Neg(VECTOR<LOCAL_P_CODE_ID,2> c,int neg)
{
    for(auto&i:c) if(i) i=LOCAL_P_CODE_ID(Value(i)^neg);
    return c;
}

int Missing_Vertex(VECTOR<PARTICLE_ID,3> t, VECTOR<PARTICLE_ID,2> e)
{
    for(int i=0;i<3;i++)
        if(!e.Contains(t(i)))
            return i;
    PHYSBAM_FATAL_ERROR();
}

template<class T,class TV>
void Apply_Analytic_BC(const PARSE_DATA_FEM<TV,TV>& pd,BC_ID bc_id,
    ARRAY<T,DOF_ID>& rhs,const MATRIX<T,TV::m>& visc, DOF_ID dof_u[2], const TV& X)
{
    TV V=pd.Velocity(X,bc_id),MV=visc*V;
    for(int i=0;i<2;i++) rhs(dof_u[i])-=MV(i);
}

template<class T,class TV>
void Apply_Analytic_BC(const PARSE_DATA_FEM<TV,TV>& pd,BC_ID bc_id,
    ARRAY<T,DOF_ID>& rhs, const TV& pres, DOF_ID dof_p, const TV& X)
{
    rhs(dof_p)+=pres.Dot(pd.Velocity(X,bc_id));
}

template<class T,class TV>
void Apply_Analytic_BC(const PARSE_DATA_FEM<TV,TV>& pd,BC_ID bc_id,
    ARRAY<T,DOF_ID>& rhs,MATRIX<LOCAL_V_CODE_ID,TV::m> visc, 
    DOF_ID dof_u[2], const TV& X, PIPE_ID pipe,
    const ARRAY<T,CODE_ID>& code_values)
{
    MATRIX<T,TV::m> visc_mat;
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            visc_mat(i,j)=code_values(Code(pipe,visc(i,j)));
    Apply_Analytic_BC(pd,bc_id,rhs,visc_mat,dof_u,X);
}

template<class T,class TV>
void Apply_Analytic_BC(const PARSE_DATA_FEM<TV,TV>& pd,BC_ID bc_id,
    ARRAY<T,DOF_ID>& rhs,VECTOR<LOCAL_P_CODE_ID,TV::m> pres,
    DOF_ID dof_p, const TV& X, PIPE_ID pipe,
    const ARRAY<T,CODE_ID>& code_values)
{
    TV pres_vec;
    for(int i=0;i<2;i++)
        pres_vec(i)=-code_values(Code(pipe,pres[i]));
    Apply_Analytic_BC(pd,bc_id,rhs,pres_vec,dof_p,X);
}

template<class T,class TV>
void Add_To_Matrix(const PARSE_DATA_FEM<TV,TV>& pd,
    ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    const ARRAY<T,CODE_ID>& code_values, MATRIX<LOCAL_V_CODE_ID,TV::m> visc, 
    DOF_ID dof_u[2], DOF_ID dof_v[2], BC_ID bc_u, BC_ID bc_v,
    ARRAY<T,DOF_ID>& rhs, const TV& Xu, const TV& Xv, PIPE_ID pipe, bool diag)
{
    if(dof_u[0]>=DOF_ID() && dof_v[0]>=DOF_ID())
    {
        for(int a=0;a<2;a++)
            for(int b=0;b<2;b++)
            {
                DOF_ID d0=dof_u[a];
                DOF_ID d1=dof_v[b];
                if(visc(a,b))
                {
                    CODE_ID code=Code(pipe,visc(a,b));
                    coded_entries.Append({d0,d1,code});
                    if(!diag) coded_entries.Append({d1,d0,code});
                }
            }
    }
    else if(dof_u[0]>=DOF_ID())
    {
        Apply_Analytic_BC(pd,bc_v,rhs,visc,dof_u,Xv,pipe,code_values);
    }
    else if(dof_v[0]>=DOF_ID() && !diag)
    {
        Apply_Analytic_BC(pd,bc_u,rhs,visc.Transposed(),dof_v,Xu,pipe,code_values);
    }
}

template<class T,class TV>
void Add_To_Matrix(const PARSE_DATA_FEM<TV,TV>& pd,
    ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values, const MATRIX<T,TV::m>& visc,
    DOF_ID dof_u[2], DOF_ID dof_v[2], BC_ID bc_u, BC_ID bc_v,
    ARRAY<T,DOF_ID>& rhs, const TV& Xu, const TV& Xv, bool diag)
{
    if(dof_u[0]>=DOF_ID() && dof_v[0]>=DOF_ID())
    {
        for(int a=0;a<2;a++)
            for(int b=(diag?a:0);b<2;b++)
            {
                DOF_ID d0=dof_u[a];
                DOF_ID d1=dof_v[b];
                CODE_ID code=code_values.Append(visc(a,b));
                coded_entries.Append({d0,d1,code});
                if(!diag || a!=b) coded_entries.Append({d1,d0,code});
            }
    }
    else if(dof_u[0]>=DOF_ID())
    {
        Apply_Analytic_BC(pd,bc_v,rhs,visc,dof_u,Xv);
    }
    else if(dof_v[0]>=DOF_ID() && !diag)
    {
        Apply_Analytic_BC(pd,bc_u,rhs,visc.Transposed(),dof_v,Xu);
    }
}

template<class T,class TV>
void Add_To_Matrix(const PARSE_DATA_FEM<TV,TV>& pd,
    ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    const ARRAY<T,CODE_ID>& code_values, VECTOR<LOCAL_P_CODE_ID,TV::m> pres,
    DOF_ID dof_u[2], DOF_ID dof_p, BC_ID bc_id,
    ARRAY<T,DOF_ID>& rhs, const TV& Xu, PIPE_ID pipe)
{
    pres=Neg(pres);
    if(dof_u[0]>=DOF_ID())
    {
        for(int b=0;b<2;b++)
        {
            DOF_ID d0=dof_p;
            DOF_ID d1=dof_u[b];
            if(pres(b))
            {
                CODE_ID code=Code(pipe,pres(b));
                coded_entries.Append({d0,d1,code});
                coded_entries.Append({d1,d0,code});
            }
        }
    }
    else
    {
        Apply_Analytic_BC(pd,bc_id,rhs,pres,dof_p,Xu,pipe,code_values);
    }
}

template<class T,class TV>
void Add_To_Matrix(const PARSE_DATA_FEM<TV,TV>& pd,
    ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values, const TV& pres,
    DOF_ID dof_u[2], DOF_ID dof_p, BC_ID bc_id,
    ARRAY<T,DOF_ID>& rhs, const TV& Xu)
{
    if(dof_u[0]>=DOF_ID())
    {
        for(int b=0;b<2;b++)
        {
            DOF_ID d0=dof_p;
            DOF_ID d1=dof_u[b];
            CODE_ID code=code_values.Append(-pres(b));
            coded_entries.Append({d0,d1,code});
            coded_entries.Append({d1,d0,code});
        }
    }
    else
    {
        Apply_Analytic_BC(pd,bc_id,rhs,pres,dof_p,Xu);
    }
}

template<class T,class TV>
void Generate_Discretization(ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values,const FLUID_LAYOUT_FEM<TV>& fl,
    const PARSE_DATA_FEM<TV,TV>& pd,T mu,ARRAY<T,DOF_ID>& rhs)
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
        BC_ID bc_id(-1);
        fl.particle_bc_map.Get(p,bc_id);
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
        TV X=fl.X(p);
        
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
                int sign=(pipe>=PIPE_ID() && (t-fl.pipes(pipe).first_element)%2)?-1:1;
                pres+=fem.pressure[k][k]*sign;
            }

            Add_To_Matrix(pd,coded_entries,code_values,visc,dof_u,dof_u,bc_id,bc_id,rhs,X,X,true);
            Add_To_Matrix(pd,coded_entries,code_values,pres,dof_u,dof_p,bc_id,rhs,X);
        }
        else
        {
            Add_To_Matrix(pd,coded_entries,code_values,Vertex_Table_Visc(elems.m),
                dof_u,dof_u,bc_id,bc_id,rhs,X,X,regular,true);

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

                Add_To_Matrix(pd,coded_entries,code_values,Vertex_Table_Pres(x,z),
                    dof_u,dof_p,bc_id,rhs,X,regular);
            }
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
        BC_ID bc_id[3];
        for(int i=0;i<2;i++)
        {
            dof_p[i]=fl.pressure_dofs(v0[i]);
            dof_u[i][0]=fl.vel_node_dofs(v0[i]);
            dof_u[i][1]=fl.vel_node_dofs(v0[i])+1;
            bc_id[i]=BC_ID(-1);
            fl.particle_bc_map.Get(v0[i],bc_id[i]);
        }
        dof_u[2][0]=fl.vel_edge_dofs(e);
        dof_u[2][1]=fl.vel_edge_dofs(e)+1;
        bc_id[2]=BC_ID(-1);
        fl.bc_map.Get(fl.Edge(e).Sorted(),bc_id[2]);

        assert(dof_p[0]>=DOF_ID());
        assert(dof_p[1]>=DOF_ID());
        TV X[3]={fl.X(v0[0]),fl.X(v0[1])};
        X[2]=(T).5*(X[0]+X[1]);

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
                int sign=(pipe>=PIPE_ID() && (t-fl.pipes(pipe).first_element)%2)?-1:1;
                for(int m=0;m<3;m++)
                    for(int n=0;n<3;n++)
                    {
                        if(m<2 && m==n) continue;
                        visc[m][n]+=fem.viscosity[di[m]][di[n]];
                    }
                for(int m=0;m<2;m++)
                    for(int n=0;n<3;n++){
                        if(m==n) continue;
                        pres[m][n]+=fem.pressure[di[m]][di[n]]*sign;}
            }

            for(int m=0;m<3;m++)
                for(int n=m;n<3;n++)
                    if(m>=2 || m!=n)
                        Add_To_Matrix(pd,coded_entries,code_values,visc[m][n],dof_u[m],dof_u[n],bc_id[m],bc_id[n],rhs,X[m],X[n],m==n);

            for(int m=0;m<2;m++)
                for(int n=0;n<3;n++)
                    if(m!=n)
                        Add_To_Matrix(pd,coded_entries,code_values,pres[m][n],dof_u[n],dof_p[m],bc_id[n],rhs,X[n]);
        }
        else
        {
            TRIANGLE_ID t0=tris(0);
            auto tri0=fl.Triangle(t0);
            int c0=Missing_Vertex(tri0,fl.Edge(e));
            int di[3]={(c0+1)%3,(c0+2)%3,c0+3};
            int sign=(regular>=PIPE_ID() && Value(t0-fl.pipes(regular).first_element)%2);
            MATRIX<LOCAL_V_CODE_ID,TV::m> visc[3][3];
            VECTOR<LOCAL_P_CODE_ID,TV::m> pres[2][3];
            memset(visc,-1,sizeof(visc));
            memset(pres,-1,sizeof(pres));

            if(tris.m==1)
            {
                for(int m=0;m<3;m++)
                    for(int n=0;n<3;n++)
                    {
                        if(m<2 && m==n) continue;
                        visc[m][n]=Main_Table_Visc(di[m],di[n]);
                    }

                for(int m=0;m<2;m++)
                    for(int n=0;n<3;n++)
                    {
                        if(m==n) continue;
                        pres[m][n]=Cond_Neg(Main_Table_Pres(di[m],di[n]),sign);
                    }
            }
            else
            {
                assert(tris.m==2);
                TRIANGLE_ID t1=tris(1);
                auto tri1=fl.Triangle(t1);
                int c1=Missing_Vertex(tri1,fl.Edge(e));

                for(int m=0;m<3;m++)
                    for(int n=0;n<3;n++)
                    {
                        if(m<2 && m==n) continue;
                        visc[m][n]=Edge_Table_Visc(c0,c1,m,n);
                    }
                
                for(int m=0;m<2;m++)
                    for(int n=0;n<3;n++)
                    {
                        if(m==n) continue;
                        pres[m][n]=Cond_Neg(Edge_Table_Pres(c0,c1,m,n),sign);
                    }
            }

            for(int m=0;m<3;m++)
                for(int n=m;n<3;n++)
                    if(m>=2 || m!=n)
                        Add_To_Matrix(pd,coded_entries,code_values,visc[m][n],dof_u[m],dof_u[n],bc_id[m],bc_id[n],rhs,X[m],X[n],regular,m==n);

            for(int m=0;m<2;m++)
                for(int n=0;n<3;n++)
                    if(m!=n)
                        Add_To_Matrix(pd,coded_entries,code_values,pres[m][n],dof_u[n],dof_p[m],bc_id[n],rhs,X[n],regular);
        }

        if(tris.m==1)
        {
            if(bc_id[2]!=BC_ID(-1) && pd.bc(bc_id[2]).type==traction)
            {
                VECTOR<TV,3> bc;
                TV YA=fl.X(tri0((c0+1)%3));
                TV YB=fl.X(tri0((c0+2)%3));
                TV N=(YB-YA).Rotate_Clockwise_90();
                for(int i=0;i<3;i++)
                    bc(i)=pd.Traction(X[i],N,mu,bc_id[2]);
                VECTOR<TV,3> B=Times_BC_NdotN(bc);
                for(int i=0;i<3;i++)
                    for(int j=0;j<2;j++)
                        if(dof_u[i][j]>=DOF_ID())
                            rhs(dof_u[i][j])+=B(i)(j);
            }
        }
    }

    for(TRIANGLE_ID t(0);t<fl.Number_Triangles();t++)
    {
        PIPE_ID regular=fl.blocks(fl.elem_data(t).block_id).pipe_id;
        auto tri=fl.Triangle(t);

        DOF_ID dof_p[3];
        DOF_ID dof_u[6][2];
        TV X[6];
        BC_ID bc_id[6];
        for(int i=0;i<3;i++)
        {
            dof_p[i]=fl.pressure_dofs(tri(i));
            dof_u[i][0]=fl.vel_node_dofs(tri(i));
            dof_u[i][1]=fl.vel_node_dofs(tri(i))+1;
            X[i]=fl.X(tri(i));
            bc_id[i]=BC_ID(-1);
            fl.particle_bc_map.Get(tri(i),bc_id[i]);
        }
        for(auto e:fl.Triangle_Edges(t))
        {
            int c=Missing_Vertex(tri,fl.Edge(e));
            dof_u[c+3][0]=fl.vel_edge_dofs(e);
            dof_u[c+3][1]=fl.vel_edge_dofs(e)+1;
            bc_id[c+3]=BC_ID(-1);
            fl.bc_map.Get(fl.Edge(e).Sorted(),bc_id[c+3]);
        }
        for(int i=0;i<3;i++)
            X[i+3]=(T).5*(X[(i+1)%3]+X[(i+2)%3]);

        if(regular<PIPE_ID())
        {
            const general_fem_coefficients<T>& fem=irregular_data(t);
            for(int m=0;m<3;m++)
            {
                Add_To_Matrix(pd,coded_entries,code_values,fem.viscosity[m+3][(m+1)%3+3],
                    dof_u[m+3],dof_u[(m+1)%3+3],bc_id[m+3],bc_id[(m+1)%3+3],rhs,X[m+3],X[(m+1)%3+3],false);
                Add_To_Matrix(pd,coded_entries,code_values,fem.viscosity[m][m+3],
                    dof_u[m],dof_u[m+3],bc_id[m],bc_id[m+3],rhs,X[m],X[m+3],false);
            }

            for(int m=0;m<3;m++)
                Add_To_Matrix(pd,coded_entries,code_values,fem.pressure[m][m+3],
                    dof_u[m+3],dof_p[m],bc_id[m+3],rhs,X[m+3]);
        }
        else
        {
            int sign=(t-fl.pipes(regular).first_element)%2;

            for(int m=0;m<3;m++)
            {
                Add_To_Matrix(pd,coded_entries,code_values,Main_Table_Visc(m+3,(m+1)%3+3),
                    dof_u[m+3],dof_u[(m+1)%3+3],bc_id[m+3],bc_id[(m+1)%3+3],rhs,X[m+3],X[(m+1)%3+3],regular,false);
                Add_To_Matrix(pd,coded_entries,code_values,Main_Table_Visc(m,m+3),
                    dof_u[m],dof_u[m+3],bc_id[m],bc_id[m+3],rhs,X[m],X[m+3],regular,false);
                Add_To_Matrix(pd,coded_entries,code_values,Cond_Neg(Main_Table_Pres(m,m+3),sign),
                    dof_u[m+3],dof_p[m],bc_id[m+3],rhs,X[m+3],regular);
            }
        }

        VECTOR<TV,6> f;
        for(int i=0;i<6;i++)
            f(i)=pd.Force(X[i],mu);
        T area=fl.Area(t);
        VECTOR<TV,6> F=Times_force_NdotN(f, area);
        for(int i=0;i<6;i++)
            for(int j=0;j<2;j++)
                if(dof_u[i][j]>=DOF_ID())
                    rhs(dof_u[i][j])+=F(i)(j);
        if(pd.analytic_velocity && pd.analytic_pressure)
        {
            VECTOR<T,6> div;
            for(int i=0;i<6;i++)
                div(i)=pd.Divergence(X[i]);
            VECTOR<T,3> D=Times_div_PdotN(div, area);
            for(int i=0;i<3;i++)
                rhs(dof_p[i])-=D[i];
        }
    }
}

template<class T,class TV>
void Solve_And_Display_Solution(const FLUID_LAYOUT_FEM<TV>& fl,const PARSE_DATA_FEM<TV,TV>& pd,
    const SYSTEM_MATRIX_HELPER<T>& MH,const ARRAY<T,DOF_ID>& rhs_vector,
    ARRAY<T,DOF_ID>* sol_out)
{
    SPARSE_MATRIX_FLAT_MXN<T> M;
    MH.Set_Matrix(Value(fl.num_dofs),Value(fl.num_dofs),M);
    
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > KRY_VEC;
    typedef MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<T>,T,KRY_VEC> KRY_MAT;
    KRY_MAT sys(M);
    KRY_VEC rhs,sol;
    rhs.v=reinterpret_cast<const ARRAY<T>&>(rhs_vector);
    sol.v.Resize(Value(fl.num_dofs));

    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    sys.Test_System(sol);
    OCTAVE_OUTPUT<T>("M.txt").Write("M",sys,rhs);
    OCTAVE_OUTPUT<T>("b.txt").Write("b",rhs);

    MINRES<T> mr;
    bool converged=mr.Solve(sys,sol,rhs,av,1e-8,0,100000);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    OCTAVE_OUTPUT<T>("x.txt").Write("x",sol);
    if(sol_out) reinterpret_cast<ARRAY<T>&>(*sol_out)=sol.v;
    
    for(PARTICLE_ID i(0);i<fl.Number_Particles();i++){
        int dof=Value(fl.vel_node_dofs(i));
        TV v=dof<0?pd.Velocity(fl.X(i),fl.particle_bc_map.Get(i)):TV(sol.v(dof),sol.v(dof+1));
        Add_Debug_Particle(fl.X(i),VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>("V",v);}
    for(EDGE_ID i(0);i<fl.Number_Edges();i++){
        auto X=.5*(fl.X(fl.Edge(i)(0))+fl.X(fl.Edge(i)(1)));
        int dof=Value(fl.vel_edge_dofs(i));
        TV v=dof<0?pd.Velocity(X,fl.bc_map.Get(fl.Edge(i).Sorted())):TV(sol.v(dof),sol.v(dof+1));
        Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>("V",v);}
    fl.Dump_Mesh();
    Flush_Frame("minres solve");
}

template void Generate_Discretization<double,VECTOR<double,2> >(
    ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID>,int>&,ARRAY<double,CODE_ID>&,
    FLUID_LAYOUT_FEM<VECTOR<double,2> > const&,PARSE_DATA_FEM<VECTOR<double,2>,VECTOR<double,2> > const&,double,
    ARRAY<double,DOF_ID>&);
template void Solve_And_Display_Solution<double,VECTOR<double,2> >(FLUID_LAYOUT_FEM<VECTOR<double,2> > const&,
    PARSE_DATA_FEM<VECTOR<double,2>,VECTOR<double,2> > const&,
    SYSTEM_MATRIX_HELPER<double> const&,ARRAY<double,DOF_ID> const&,ARRAY<double,DOF_ID>*);
}
