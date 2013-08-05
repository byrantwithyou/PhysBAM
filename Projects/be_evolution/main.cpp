//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <climits>

using namespace PhysBAM;

template<class TV>
class MINIMIZATION_OBJECTIVE:public NONLINEAR_FUNCTION<typename TV::SCALAR(KRYLOV_VECTOR_BASE<typename TV::SCALAR>&)>, public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    T dt,time,a,b,c;
    ARRAY<TWIST<TV> > fake_rigid_x;
    GENERALIZED_VELOCITY<TV> x1,v1,&x0,&v0,&d;
    GENERALIZED_VELOCITY<TV> &tmp0,&tmp1;
    BACKWARD_EULER_SYSTEM<TV> system;

    MINIMIZATION_OBJECTIVE(SOLID_BODY_COLLECTION<TV>& solid_body_collection,T dt,T time)
        :KRYLOV_SYSTEM_BASE<T>(false,false),solid_body_collection(solid_body_collection),dt(dt),time(time),
        fake_rigid_x(solid_body_collection.rigid_body_collection.rigid_body_particles.twist.m),
        x1(solid_body_collection.deformable_body_collection.particles.X,fake_rigid_x,solid_body_collection),
        v1(solid_body_collection.deformable_body_collection.particles.V,solid_body_collection.rigid_body_collection.rigid_body_particles.twist,solid_body_collection),
        x0(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),v0(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),
        d(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),tmp0(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),
        tmp1(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),system(0,solid_body_collection,dt,time,time,0,0,0,true,true)
    {
        x0=x1;
        v0=v1;
        a=1/(dt*dt);
        b=-a;
        c=-1/dt;
        d.Copy(b,x0);
        d.Copy(c,v0,d);
    }

    virtual ~MINIMIZATION_OBJECTIVE()
    {
        delete &x0;
        delete &v0;
        delete &d;
        delete &tmp0;
        delete &tmp1;
    }

    virtual void Compute(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const
    {
        PHYSBAM_ASSERT(!h || h==this);

        tmp0.Copy(a,x,d);

        const_cast<MINIMIZATION_OBJECTIVE<TV>*>(this)->v1.Copy(-1,x0,x);
        const_cast<MINIMIZATION_OBJECTIVE<TV>*>(this)->v1*=1/dt;
        const_cast<MINIMIZATION_OBJECTIVE<TV>*>(this)->x1.Copy(1,x);
        solid_body_collection.Update_Position_Based_State(time,true);

        if(e){
            T ke=0,pe=0;
            solid_body_collection.Compute_Energy(time,ke,pe);
            *e=system.Inner_Product(tmp0,tmp0)/(2*a)+pe;}

        if(g){
            tmp1*=0;
            solid_body_collection.Add_Velocity_Independent_Forces(tmp1.V.array,tmp1.rigid_V.array,time);
            system.projection_data.mass.Multiply(tmp0,debug_cast<GENERALIZED_VELOCITY<TV>&>(*g),false);
            *g-=tmp1;}
    }

    void Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
    {
        const GENERALIZED_VELOCITY<TV>& V=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV);
        GENERALIZED_VELOCITY<TV>& F=debug_cast<GENERALIZED_VELOCITY<TV>&>(BF);
        solid_body_collection.Implicit_Velocity_Independent_Forces(V.V.array,V.rigid_V.array,F.V.array,F.rigid_V.array,-1,time);
        system.projection_data.mass.Multiply(V,tmp0,false);
        F.Copy(a,tmp0,F);
    }

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
    {
        const GENERALIZED_VELOCITY<TV>& V1=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV1),&V2=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV2);
        return V1.V.Dot_Double_Precision(V2.V)+V1.rigid_V.Dot_Double_Precision(V2.rigid_V);
    }

    typename TV::SCALAR Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
    {
        return sqrt(Inner_Product(BR,BR));
    }

    void Project(KRYLOV_VECTOR_BASE<T>& V) const {}
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const {}
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const {}
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const {}
};

template<class TV>
class SIMULATION
{
    typedef typename TV::SCALAR T;
public:
    SOLID_BODY_COLLECTION<TV> solid_body_collection;
    T time;
    NEWTONS_METHOD<T> nm;

    SIMULATION()
    {
        nm.max_iterations=100000;
        nm.max_krylov_iterations=2000;
        nm.krylov_tolerance=1;
        nm.fail_on_krylov_not_converged=false;
        nm.tolerance=1e-5;
        nm.angle_tolerance=1e-2;
    }

    ~SIMULATION() {}

    void Advance_One_Time_Step_Position(const T dt)
    {
        LOG::SCOPE scope("Advance_One_Time_Step_Position");
        MINIMIZATION_OBJECTIVE<TV> obj(solid_body_collection,dt,time);
        KRYLOV_VECTOR_BASE<T>* x0 = obj.x0.Clone_Default();
        *x0=obj.x0;
        obj.Test(*x0,obj);

        bool converged=nm.Newtons_Method(obj,obj,*x0);
        PHYSBAM_ASSERT(converged);
        delete x0;
    }
};

typedef double T;
typedef float RW;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;
RW rw=RW();
STREAM_TYPE stream_type(rw); // gcc 3.3.2 workaround
std::string output_directory="output";

VIEWER_OUTPUT<TV>* pvo=0;
SIMULATION<TV>* psim=0;

void Flush_State(const char* str)
{
    static int frame=0;
    printf("Flush %i, '%s'\n",frame,str);
    pvo->Flush_Frame(STRING_UTILITIES::string_sprintf(str,frame).c_str());
    psim->solid_body_collection.Write(stream_type,output_directory,frame,-1,frame==0,true,true,true,false);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
    frame++;
}
extern void (*NM_Flush_State)(const char*);

int main(int argc,char* argv[])
{
    NM_Flush_State=&Flush_State;
    SIMULATION<TV> simulation;
    psim=&simulation;
    int res=6,seed=-1,steps=10;
    bool enforce_definiteness=false,do_pt=false;
    T dt=.1;

    PARSE_ARGS parse_args(argc,argv);
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;
    parse_args.Add_Not("-mr",&simulation.nm.use_cg,"use minres instead of cg");
    parse_args.Add("-o",&output_directory,"dir","output directory");
    parse_args.Add("-resolution",&res,"res","resolution");
    parse_args.Add("-enf_def",&enforce_definiteness,"enforce definiteness in system");
    parse_args.Add("-kry_it",&simulation.nm.max_krylov_iterations,"iter","maximum iterations for Krylov solver");
    parse_args.Add("-kry_tol",&simulation.nm.krylov_tolerance,"tol","tolerance for Krylov solver");
    parse_args.Add("-newton_it",&simulation.nm.max_iterations,"iter","maximum iterations for Newton");
    parse_args.Add("-newton_tol",&simulation.nm.tolerance,"tol","tolerance for Newton");
    parse_args.Add("-kry_fail",&simulation.nm.fail_on_krylov_not_converged,"terminate if Krylov solver fails to converge");
    parse_args.Add("-angle_tol",&simulation.nm.angle_tolerance,"tol","gradient descent tolerance");
    parse_args.Add("-seed",&seed,"fixed seed","set random seed");
    parse_args.Add("-dt",&dt,"step","time step size");
    parse_args.Add("-pt",&do_pt,"point test");
    parse_args.Add("-steps",&steps,"steps","number of time steps");
    parse_args.Add_Not("-gss",&simulation.nm.use_wolfe_search,"use golden section search instead of wolfe conditions line search");
    parse_args.Parse();
    if(!simulation.nm.use_wolfe_search) simulation.nm.use_golden_section_search=true;

    LOG::cout<<std::setprecision(16);

    GRID<TV> grid(TV_INT(),RANGE<TV>::Unit_Box());
    VIEWER_OUTPUT<TV> vo(stream_type,grid,output_directory);
    pvo=&vo;
    SOLIDS_STANDARD_TESTS<TV> tests(stream_type,getenv("PHYSBAM_DATA_DIRECTORY"),simulation.solid_body_collection);
    GRID<TV> cube_grid(TV_INT()+res,RANGE<TV>::Centered_Box());
    TETRAHEDRALIZED_VOLUME<T>& tv=tests.Create_Mattress(cube_grid);
    simulation.solid_body_collection.Add_Force(Create_Finite_Volume(tv,new COROTATED_FIXED<T,TV::m>(1e6,0.3,0)));
    for(int i=0;i<simulation.solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
        simulation.solid_body_collection.deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;
    if(enforce_definiteness) simulation.solid_body_collection.Enforce_Definiteness(true);

    RANDOM_NUMBERS<T> random;
    if(seed!=-1) random.Set_Seed(seed);
    if(do_pt) simulation.solid_body_collection.deformable_body_collection.particles.X.Fill(TV());
    else random.Fill_Uniform(simulation.solid_body_collection.deformable_body_collection.particles.X,-1,1);

    simulation.solid_body_collection.Update_Simulated_Particles();
    Flush_State("frame %d");
    for(int frame=1;frame<=steps;frame++)
    {
        simulation.Advance_One_Time_Step_Position(dt);
        Flush_State("frame %d");
    }

    LOG::Finish_Logging();
    return 0;
}
//#####################################################################
