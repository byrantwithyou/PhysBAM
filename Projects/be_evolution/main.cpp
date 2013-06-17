//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
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

    SIMULATION()
    {
    }

    ~SIMULATION() {}

    void Advance_One_Time_Step_Position(const T dt)
    {
        MINIMIZATION_OBJECTIVE<TV> obj(solid_body_collection,dt,time);
        KRYLOV_VECTOR_BASE<T>* x0 = obj.x0.Clone_Default();
        *x0=obj.x0;
        obj.Test(*x0,obj);

        NEWTONS_METHOD<T> nm;
        nm.max_iterations=100000;
        nm.max_krylov_iterations=2000;
        nm.krylov_tolerance=1e-3;
        nm.fail_on_krylov_not_converged=false;
        bool converged=nm.Newtons_Method(obj,obj,*x0);
        PHYSBAM_ASSERT(converged);
        delete x0;
    }
};

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    RW rw=RW();STREAM_TYPE stream_type(rw); // gcc 3.3.2 workaround
    std::string output_directory="output";
    LOG::cout<<std::setprecision(16);

    GRID<TV> grid(TV_INT(),RANGE<TV>::Unit_Box());
    VIEWER_OUTPUT<TV> vo(stream_type,grid,output_directory);
    SIMULATION<TV> simulation;
    SOLIDS_STANDARD_TESTS<TV> tests(stream_type,getenv("PHYSBAM_DATA_DIRECTORY"),simulation.solid_body_collection);
    GRID<TV> cube_grid(TV_INT()+6,RANGE<TV>::Centered_Box());
    TETRAHEDRALIZED_VOLUME<T>& tv=tests.Create_Mattress(cube_grid);
    simulation.solid_body_collection.Add_Force(Create_Finite_Volume(tv,new COROTATED_FIXED<T,3>(1e6,0.3,0)));
    for(int i=0;i<simulation.solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
        simulation.solid_body_collection.deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;


    RANDOM_NUMBERS<T> random;
    random.Fill_Uniform(simulation.solid_body_collection.deformable_body_collection.particles.X,-1,1);


    simulation.solid_body_collection.Update_Simulated_Particles();

    T dt=.1;
    vo.Flush_Frame(STRING_UTILITIES::string_sprintf("frame %d",0).c_str());
    simulation.solid_body_collection.Write(stream_type,output_directory,0,-1,true,true,true,true,false);
    for(int frame=1;frame<10;frame++)
    {
        simulation.Advance_One_Time_Step_Position(dt);
        vo.Flush_Frame(STRING_UTILITIES::string_sprintf("frame %d",frame).c_str());
        simulation.solid_body_collection.Write(stream_type,output_directory,frame,-1,true,true,true,true,false);
    }

    return 0;
}
//#####################################################################
