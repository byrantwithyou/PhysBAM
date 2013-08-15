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
#include <Tools/Read_Write/TYPED_STREAM.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <climits>
#include "MINIMIZATION_OBJECTIVE.h"
#include "SIMULATION.h"
#include <boost/function.hpp>
using namespace PhysBAM;

extern boost::function<void(const char*)> NM_Flush_State;

template<class TV>
class VIEWER_WRAPPER
{
public:
    VIEWER_OUTPUT<TV> vo;
    SIMULATION<TV>& sim;
    int frame;

    VIEWER_WRAPPER(STREAM_TYPE stream_type,const GRID<TV>& grid,const std::string& output_directory,SIMULATION<TV>& sim)
        :vo(stream_type,grid,output_directory),sim(sim),frame(0)
    {
    }

    void Flush_State(const char* str)
    {
        printf("Flush %i, '%s'\n",frame,str);
        vo.Flush_Frame(STRING_UTILITIES::string_sprintf(str,frame).c_str());
        sim.solid_body_collection.Write(vo.stream_type,vo.output_directory,frame,-1,frame==0,true,true,true,false);
        FILE_UTILITIES::Write_To_Text_File(vo.output_directory+"/common/last_frame",frame,"\n");
        frame++;
    }
};

template<class T,class TV,class TV_INT> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT&
Insert_Cube(SOLIDS_STANDARD_TESTS<TV>& tests,TV_INT resolution,const RANGE<TV>& range,T stiffness,T poissons_ratio)
{
    GRID<TV> cube_grid(resolution,range);
    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT& obj=tests.Create_Mattress(cube_grid);
    tests.solid_body_collection.Add_Force(Create_Finite_Volume(obj,new COROTATED_FIXED<T,TV::m>(stiffness,poissons_ratio,0)));
    return obj;
}

template<class TV,class T> IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>*
Add_Ground(SOLIDS_STANDARD_TESTS<TV>& tests,T stiffness=(T)1e4,T separation_parameter=(T)1e-4,
    T length_scale=1,T y=0,TV normal=TV::Axis_Vector(1))
{
    RIGID_BODY<TV>& rigid_body=tests.Add_Analytic_Plane(normal,normal*y,10,5);
    IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>* penalty_force=
        new IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>(tests.solid_body_collection.deformable_body_collection.particles,
            rigid_body.implicit_object,stiffness,separation_parameter,length_scale);
     tests.solid_body_collection.Add_Force(penalty_force);
    return penalty_force;
}

template<class T>
void Init_Test(SIMULATION<VECTOR<T,2> >& simulation,STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
}

template<class T>
void Init_Test(SIMULATION<VECTOR<T,3> >& simulation,STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,TV::m> TV_INT;

    int test_number=0,resolution=6,seed=-1;
    T stiffness=(T)1e6,poissons_ratio=(T).3;
    TV nonuniform_scale=TV((T).5,(T).9,(T)1.2);
    parse_args.Extra_Optional(&test_number,"example number","example number to run");
    parse_args.Add("-seed",&seed,"fixed seed","set random seed");
    parse_args.Add("-resolution",&resolution,"res","resolution");
    parse_args.Add("-stiff",&stiffness,"stiffness","constitutive model stiffness");
    parse_args.Add("-pr",&poissons_ratio,"ratio","constitutive model poissons ratio");
    parse_args.Add("-nonuniform_scale",&nonuniform_scale,"scale","scale factor");
    parse_args.Parse();

    SOLIDS_STANDARD_TESTS<TV> tests(stream_type,getenv("PHYSBAM_DATA_DIRECTORY"),simulation.solid_body_collection);
    RANDOM_NUMBERS<T> random;
    if(seed!=-1) random.Set_Seed(seed);
    DEFORMABLE_PARTICLES<TV>& particles=simulation.solid_body_collection.deformable_body_collection.particles;

    switch(test_number){
        case 0:
            Insert_Cube(tests,TV_INT()+resolution,RANGE<TV>::Centered_Box(),stiffness,poissons_ratio);
            break;

        case 1:
            random.Fill_Uniform(Insert_Cube(tests,TV_INT()+resolution,RANGE<TV>::Centered_Box(),stiffness,poissons_ratio).particles.X,-1,1);
            break;

        case 2:
            Insert_Cube(tests,TV_INT()+resolution,RANGE<TV>::Centered_Box(),stiffness,poissons_ratio).particles.X.Fill(TV());
            break;

        case 3:{
            ARRAY_VIEW<TV> X=Insert_Cube(tests,TV_INT()+resolution,RANGE<TV>::Centered_Box(),stiffness,poissons_ratio).particles.X;
            for(int i=0;i<X.m;i++) X(i)*=nonuniform_scale;}
            break;

        case 4:{
            TETRAHEDRALIZED_VOLUME<T>& tv=Insert_Cube(tests,TV_INT()+resolution,RANGE<TV>::Centered_Box(),stiffness,poissons_ratio);
            simulation.solid_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(particles,tv.mesh));
            Add_Ground(tests,(T)1e4,(T)2e-2,(T).1,(T)-2);}
            break;

        case 5:{
            TETRAHEDRALIZED_VOLUME<T>& tv=tests.Create_Tetrahedralized_Volume(tests.data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,1000);
            simulation.solid_body_collection.Add_Force(Create_Finite_Volume(tv,new COROTATED_FIXED<T,TV::m>(stiffness,poissons_ratio,0)));
            simulation.solid_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(particles,tv.mesh));
            Add_Ground(tests,(T)1e4,(T)2e-2,(T).1,(T)0);}
            break;

        default:
            LOG::cout<<"test number not implemented: "<<test_number<<std::endl;
            PHYSBAM_FATAL_ERROR();}
}

template<class TV>
void Integration_Test(int argc,char* argv[],PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef float RW;
    typedef VECTOR<int,TV::m> TV_INT;
    STREAM_TYPE stream_type((RW()));
    std::string output_directory="output";
    bool enforce_definiteness=false;

    int frames=120,steps_per_frame=10;
    T frame_dt=1./24,dt=frame_dt/steps_per_frame;
    SIMULATION<TV> simulation;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection=simulation.solid_body_collection;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    parse_args.Add_Not("-mr",&simulation.nm.use_cg,"use minres instead of cg");
    parse_args.Add("-o",&output_directory,"dir","output directory");
    parse_args.Add("-kry_it",&simulation.nm.max_krylov_iterations,"iter","maximum iterations for Krylov solver");
    parse_args.Add("-kry_tol",&simulation.nm.krylov_tolerance,"tol","tolerance for Krylov solver");
    parse_args.Add("-newton_it",&simulation.nm.max_iterations,"iter","maximum iterations for Newton");
    parse_args.Add("-newton_tol",&simulation.nm.tolerance,"tol","tolerance for Newton");
    parse_args.Add("-newton_cd_tol",&simulation.nm.countdown_tolerance,"tol","tolerance for Newton");
    parse_args.Add("-kry_fail",&simulation.nm.fail_on_krylov_not_converged,"terminate if Krylov solver fails to converge");
    parse_args.Add("-angle_tol",&simulation.nm.angle_tolerance,"tol","gradient descent tolerance");
    parse_args.Add("-dt",&dt,"step","time step size");
    parse_args.Add("-frames",&frames,"frames","number of frames of simulation");
    parse_args.Add("-spf",&steps_per_frame,"steps","number of time steps per frame");
    parse_args.Add("-gss",&simulation.nm.use_golden_section_search,"use golden section search instead of wolfe conditions line search");
    parse_args.Add("-enf_def",&enforce_definiteness,"enforce definiteness in system");
    parse_args.Parse(true);

    frame_dt=dt*steps_per_frame;

    if(simulation.nm.use_golden_section_search) simulation.nm.use_wolfe_search=false;

    LOG::cout<<std::setprecision(16);

    GRID<TV> grid(TV_INT(),RANGE<TV>::Unit_Box());
    VIEWER_WRAPPER<TV> viewer_wrapper(stream_type,grid,output_directory,simulation);
    NM_Flush_State=[&](const char* str){viewer_wrapper.Flush_State(str);};

    Init_Test(simulation,stream_type,parse_args);
    for(int i=0;i<deformable_body_collection.deformables_forces.m;i++)
        deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;
    if(enforce_definiteness) solid_body_collection.Enforce_Definiteness(true);

    solid_body_collection.Update_Simulated_Particles();
    NM_Flush_State("frame %d");

    for(int frame=1;frame<=frames;frame++){
        for(int s=0;s<steps_per_frame;s++)
            simulation.Advance_One_Time_Step_Position(dt);
        NM_Flush_State("frame %d");}

    LOG::Finish_Logging();
}

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool use_2d=false;
    PARSE_ARGS parse_args(argc,argv);
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;
    parse_args.Add("-2d",&use_2d,"Use 2D");
    parse_args.Parse(true);

    if(use_2d)
        Integration_Test<VECTOR<double,2> >(argc,argv,parse_args);
    else
        Integration_Test<VECTOR<double,3> >(argc,argv,parse_args);

    return 0;
}
