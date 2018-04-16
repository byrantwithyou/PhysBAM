//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef MATRIX<T,3> MV;
typedef VECTOR<int,TV::m> TV_INT;


struct KV:public KRYLOV_VECTOR_BASE<T>
{
    MV M;
    TV V;

    KV(){}
    virtual ~KV(){}

    KRYLOV_VECTOR_BASE& operator= (const KRYLOV_VECTOR_BASE& bv)
        {
            M=((const KV&)bv).M;
            V=((const KV&)bv).V;
            return *this;
        }

    virtual KRYLOV_VECTOR_BASE& operator+=(const KRYLOV_VECTOR_BASE& bv)
    {
        M+=((const KV&)bv).M;
        V+=((const KV&)bv).V;
        return *this;
    }

    virtual KRYLOV_VECTOR_BASE& operator-=(const KRYLOV_VECTOR_BASE& bv)
    {
        M-=((const KV&)bv).M;
        V-=((const KV&)bv).V;
        return *this;
    }

    virtual KRYLOV_VECTOR_BASE& operator*=(const T a)
    {
        M*=a;
        V*=a;
        return *this;
    }

    virtual void Copy(const T c,const KRYLOV_VECTOR_BASE& bv)
    {
        M=c*((const KV&)bv).M;
        V=c*((const KV&)bv).V;
    }

    virtual void Copy(const T c1,const KRYLOV_VECTOR_BASE& bv1,const KRYLOV_VECTOR_BASE& bv2)
    {
        M=c1*((const KV&)bv1).M+((const KV&)bv2).M;
        V=c1*((const KV&)bv1).V+((const KV&)bv2).V;
    }

    virtual int Raw_Size() const {return 12;}
    virtual T& Raw_Get(int i) {if(i<9) return M.x[i];return V(i-9);}
    const T& Raw_Get(int i) const {if(i<9) return M.x[i];return V(i-9);}
    virtual KRYLOV_VECTOR_BASE* Clone_Default() const {return new KV;}
    virtual void Resize(const KRYLOV_VECTOR_BASE& v) {}
};

struct KM:public KRYLOV_SYSTEM_BASE<T>
{
    const ARRAY<TV>& X;
    IMPLICIT_OBJECT<TV>& io;
    T k;
    MV M;
    TV V;

    KM(const ARRAY<TV>& X,IMPLICIT_OBJECT<TV>& io,T k)
        :KRYLOV_SYSTEM_BASE<T>(false,true),X(X),io(io),k(k)
    {}

    virtual ~KM(){}

    virtual void Multiply(const KRYLOV_VECTOR_BASE<T>& bx,KRYLOV_VECTOR_BASE<T>& result,bool transpose=false) const
    {
        const KV& u=(const KV&)bx;
        KV& r=(KV&)result;
        MV dM,dM_;
        TV dV,dV_;
        MV M_=u.M;
        TV V_=u.V;

        for(int i=0;i<X.m;i++)
        {
            TV Y=M*X(i)+V;
            TV Y_=M_*X(i)+V_;
            T p=io(Y);
            TV N=io.Normal(Y);
            T p_=N.Dot(Y_);
            SYMMETRIC_MATRIX<T,3> H=io.Hessian(Y);
            TV N_=H*Y_;
            dV+=p*N;
            dV_+=p_*N+p*N_;
            dM+=Outer_Product(p*N,X(i));
            dM_+=Outer_Product(p_*N+p*N_,X(i));
        }
        dV/=2*X.m;
        dV_/=2*X.m;
        dM/=2*X.m;
        dM_/=2*X.m;
        MV A=M.Transpose_Times(M);
        MV A_=M_.Transpose_Times(M)+M.Transpose_Times(M_);
        MV B=A-A.Trace()/3;
        MV B_=A_-A_.Trace()/3;
        MV C=B.Twice_Symmetric_Part()-2./3*B.Trace();
        MV C_=B_.Twice_Symmetric_Part()-2./3*B_.Trace();
        dM+=k*M*C;
        dM_+=k*M_*C+k*M*C_;
        r.M=dM_;
        r.V=dV_;
    }

    virtual double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bx,const KRYLOV_VECTOR_BASE<T>& by) const
    {
        const KV& x=(const KV&)bx;
        const KV& y=(const KV&)by;
        return x.V.Dot(y.V)+x.M.Inner_Product(x.M,y.M);
    }

    virtual T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
    {return sqrt(Inner_Product(x,x));}

    virtual void Project(KRYLOV_VECTOR_BASE<T>& x) const{}
    virtual void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const{}
    virtual void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const{}
};

struct OBJECTIVE: public NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>
{
    ARRAY<TV>& X;
    T k;
    SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV> io;
    OBJECTIVE(ARRAY<TV>& X,int k,GRID<TV>& grid,ARRAY<T,TV_INT>& phi)
        :X(X),k(k),io(grid,phi)
    {}

    virtual ~OBJECTIVE(){}

    void Compute(const KRYLOV_VECTOR_BASE<T>& xx,KRYLOV_SYSTEM_BASE<T>* hh,KRYLOV_VECTOR_BASE<T>* gg,T* e) const override
    {
        const KV& x=dynamic_cast<const KV&>(xx);
        KV* g=dynamic_cast<KV*>(gg);
        KM* h=dynamic_cast<KM*>(hh);
        MV dM;
        TV dV;

        T en=0;
        for(int i=0;i<X.m;i++)
        {
            TV Y=x.M*X(i)+x.V;
            T p=io(Y);
            TV N=io.Normal(Y);
            SYMMETRIC_MATRIX<T,3> H=io.Hessian(Y);
            en+=sqr(p)/2;
            dV+=p*N;
            dM+=Outer_Product(p*N,X(i));
        }
        en/=2*X.m;
        dV/=2*X.m;
        dM/=2*X.m;
        MV A=x.M.Transpose_Times(x.M);
        // A_{ij} = M_{ki} M_{kj}
        // A_{ij,rs} = delta_{is} M_{rj} + M_{ri} delta_{js}
        MV B=A-A.Trace()/3;
        // B_{ij,rs} = delta_{is} M_{rj} + M_{ri} delta_{js} - 2 M_{rs} delta_{ij}/3
        // B_{ij} B_{ij,rs} = B_{sj} M_{rj} + B_{is} M_{ri} - 2 B_{ii} M_{rs}/3
        // B_{ij} B_{ij,rs} = M_{rk} (B_{sk} + B_{ks}) - 2 B_{ii} M_{rs}/3
        T norm=B.Frobenius_Norm_Squared();
        en+=k*norm/2;
        MV C=B.Twice_Symmetric_Part()-2./3*B.Trace();
        dM+=k*x.M*C;
        if(e) *e=en;
        if(g){g->M=dM;g->V=dV;}
        if(h){h->M=x.M;h->V=x.V;}
    }
};

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    std::string filename;
    std::string levelset_filename;
    int inside_vertex=-1;
    int outside_vertex=-1;
    int front_vertex=-1;
    bool compute_levelset=false;
    T levelset_radius=2;
    T regularity_parameter=1;
    int levelset_res=10;
    bool is_left=false;
    parse_args.Add("-i",&filename,"file","file containing bone geometry");
    parse_args.Add("-ls",&levelset_filename,"file","file containing template level set");
    parse_args.Add("-in",&inside_vertex,"vertex","index of vertex on inside top round part");
    parse_args.Add("-out",&outside_vertex,"vertex","index of vertex on outside top round part");
    parse_args.Add("-front",&front_vertex,"vertex","index of vertex on front flat part");
    parse_args.Add("-compute_levelset",&compute_levelset,"compute level set");
    parse_args.Add("-levelset_radius",&levelset_radius,"radius","scale level set radius compared with fit");
    parse_args.Add("-levelset_res",&levelset_res,"resolution","level set resolution");
    parse_args.Add("-left",&is_left,"is left leg");
    parse_args.Add_Not("-right",&is_left,"is left leg");
    parse_args.Add("-reg",&regularity_parameter,"param","regularization parameter");
    parse_args.Parse();
    LOG::Initialize_Logging(false,false,1<<30,true);

    PHYSBAM_ASSERT(inside_vertex>=0);
    PHYSBAM_ASSERT(outside_vertex>=0);
    PHYSBAM_ASSERT(front_vertex>=0);

    TRIANGULATED_SURFACE<T> tsa;
    tsa.Read_Obj(filename);
    TRIANGULATED_SURFACE<T> &ts=*tsa.Create_Compact_Copy();
    ts.Fill_Holes();
    ts.mesh.Make_Orientations_Consistent();

    TV X_i=tsa.particles.X(inside_vertex);
    TV X_o=tsa.particles.X(outside_vertex);
    TV X_f=tsa.particles.X(front_vertex);

    TV C=(X_i+X_o+X_f)/3;
    TV A_i=X_f-C;
    T a_i_len=A_i.Normalize();
    TV A_k=(X_o-X_i).Cross(A_i).Normalized();
    TV A_j=A_k.Cross(A_i);
    if(is_left) A_k=-A_k;
    MV R(A_i,A_j,A_k);
    R/=a_i_len;
    for(int i=0;i<ts.particles.X.m;i++)
        ts.particles.X(i)=R.Transpose_Times(ts.particles.X(i)-C);

    ts.Write_Obj("orient.obj");

    if(compute_levelset)
    {
        GRID<TV> ls_grid(TV_INT()+levelset_res,RANGE<TV>::Centered_Box()*levelset_radius,true);
        ARRAY<T,TV_INT> phi(ls_grid.Domain_Indices());
        LEVELSET<TV> levelset(ls_grid,phi,0);
        SIGNED_DISTANCE::Calculate(ts,ls_grid,phi,true);
        Write_To_File<double>("levelset.phi",levelset);
        return 0;
    }
    else
    {
        GRID<TV> ls_grid;
        ARRAY<T,TV_INT> phi;
        LEVELSET<TV> levelset(ls_grid,phi,0);
        Read_From_File("levelset.phi",levelset);

        ARRAY<TV> X;
        for(int i=0;i<ts.particles.X.m;i++)
            if(ts.particles.X(i).Magnitude_Squared()<=levelset_radius*levelset_radius)
                X.Append(ts.particles.X(i));

        OBJECTIVE obj(X,regularity_parameter,ls_grid,phi);

        NEWTONS_METHOD<T> nm;
        nm.tolerance=(T)1e-4;
        nm.max_iterations=100;
        nm.krylov_tolerance=(T)1e-10;
        nm.angle_tolerance=(T)1e-4;
        nm.debug=true;
        nm.require_one_iteration=true;

        KV x;
        x.M+=1;
        SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV> lio(ls_grid,phi);
        KM sys(X,lio,regularity_parameter);

        ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
        obj.Test(x,sys);
        bool conv=nm.Newtons_Method(obj,sys,x,av);
        obj.Test(x,sys);
        LOG::printf("converged: %i\n",conv);

        VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.),ls_grid,"output");
        Dump_Levelset(ls_grid,lio,TV(1,0,0));
        Flush_Frame<TV>("A");
        Dump_Levelset(ls_grid,lio,TV(1,0,0));
        Dump_Surface(ts,TV(0,1,0));
        Flush_Frame<TV>("B");
        Dump_Levelset(ls_grid,lio,TV(1,0,0));
        for(int i=0;i<X.m;i++)
            Add_Debug_Particle(X(i),TV(0,1,0));
        Flush_Frame<TV>("C");
        
    }

    LOG::Finish_Logging();
    return 0;
}

void Compute_Levelset(TRIANGULATED_SURFACE<T>& ts)
{

}

