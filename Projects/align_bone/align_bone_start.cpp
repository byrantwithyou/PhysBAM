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
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
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
            
    KM(const ARRAY<TV>& X,IMPLICIT_OBJECT<TV>& io,
        T k,const bool use_preconditioner,
        const bool preconditioner_commutes_with_projection)
        :KRYLOV_SYSTEM_BASE<T>(false,true),X(X),io(io),k(k)
    {}

    virtual ~KM(){}

    virtual void Multiply(const KRYLOV_VECTOR_BASE<T>& bx,KRYLOV_VECTOR_BASE<T>& result,bool transpose=false) const
    {
        const KV& x=(const KV&)bx;
        KV& r=(KV&)result;
        MV M,dM;
        TV V,dV;
        for(int i=0;i<X.m;i++)
        {
            TV Y=M*X(i)+V;
            TV dY=x.M*X(i)+x.V;
            T p=io(Y);
            TV N=io.Normal(Y);
            T dp=N.Dot(dY);
            SYMMETRIC_MATRIX<T,3> H=io.Hessian(Y);
            TV dN=H*dY;
            TV dU=dp*N+p*dN;
            dV+=dU;
            dM+=Outer_Product(dU,X(i));
        }
        auto A=M.Normal_Equations_Matrix();
        auto dA=dM.Transpose_Times(M).Twice_Symmetric_Part();
        A=A-A.Trace()/3;
        dA=dA-dA.Trace()/3;
        T norm=A.Frobenius_Norm_Squared();
        T dnorm=A.Inner_Product(A,dA)*2;
        auto B=A.Twice_Symmetric_Part()-2./3*A.Trace();
        auto dB=dA.Twice_Symmetric_Part()-2./3*dA.Trace();
        auto C=M*B;
        auto dC=dM*B+M*dB;
        dM+=k*(dnorm*C+norm*dC);
        r.M=dM;
        r.V=dV;
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
        MV A=x.M.Transpose_Times(x.M);
        // A_{ij} = M_{ki} M_{kj}
        // A_{ij,rs} = delta_{is} M_{rj} + M_{ri} delta_{js}
        A=A-A.Trace()/3;
        // A_{ij,rs} = delta_{is} M_{rj} + M_{ri} delta_{js} - M_{rs}  delta_{ij}/3 - M_{rs} delta_{ij}/3
        // A_{ij,rs} = delta_{is} M_{rj} + M_{ri} delta_{js} - 2 M_{rs}  delta_{ij}/3
        // A_{ij} A_{ij,rs} = A_{si} M_{ri} + M_{ri} A_{is} - 2 A_{ii} M_{rs}/3
        T norm=A.Frobenius_Norm_Squared();
        en+=k*sqr(norm)/4;
        dM+=k*norm*x.M*(A.Twice_Symmetric_Part()-2./3*A.Trace());
        if(e) *e=en;
        if(g){g->M=dM;g->V=dV;}
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

    PHYSBAM_ASSERT(inside_vertex>=0);
    PHYSBAM_ASSERT(outside_vertex>=0);
    PHYSBAM_ASSERT(front_vertex>=0);

    TRIANGULATED_SURFACE<T> ts;
    ts.Read_Obj(filename);
    ts.Fill_Holes();
    ts.mesh.Make_Orientations_Consistent();

    TV X_i=ts.particles.X(inside_vertex);
    TV X_o=ts.particles.X(outside_vertex);
    TV X_f=ts.particles.X(front_vertex);

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

        
//        NONLINEAR_FUNCTION
        
/*
  
  Choose M, S;

  1/num_pts * sum[i] phi(M * X(i) + S)^2/2 + k * Frobenius_Norm_Squared(M^T*M - Trace(M^T*M)/3)^4

  


 */        


        
    }
    
    return 0;
}

void Compute_Levelset(TRIANGULATED_SURFACE<T>& ts)
{
    
}

