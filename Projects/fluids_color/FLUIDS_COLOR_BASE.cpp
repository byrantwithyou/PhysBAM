//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_BOX.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ELLIPSOID.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ROTATE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SPHERE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_TRANSLATE.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include "ANALYTIC_POLYMER_STRESS.h"
#include "ANALYTIC_VELOCITY.h"
#include "FLUIDS_COLOR_BASE.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUIDS_COLOR_BASE<TV>::
FLUIDS_COLOR_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :PLS_FC_EXAMPLE<TV>(stream_type_input),test_number(0),resolution(32),stored_last_frame(0),user_last_frame(false),
    mu0(1),mu1(2),rho0(1),rho1(2),inv_Wi0(1.2),inv_Wi1(1.4),beta0(1.3),beta1(1.7),unit_mu(0),unit_rho(0),unit_st(0),
    unit_p(0),weiss(1),weiss_inv(1),m(1),s(1),kg(1),bc_n(false),bc_d(false),bc_s(false),test_analytic_diff(false),
    refine(1),surface_tension(0),override_rho0(false),override_rho1(false),override_mu0(false),override_mu1(false),
    override_inv_Wi0(false),override_inv_Wi1(false),override_beta0(false),override_beta1(false),
    override_surface_tension(false),use_pls_over_levelset(false),use_levelset_over_pls(false),use_test_output(false),
    analytic_initial_only(false),number_of_threads(1),override_output_directory(false),use_u0(false),
    use_u1(false),use_p0(false),use_p1(false),use_S0(false),use_S1(false)
{
    last_frame=16;
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Add("-restart",&restart,"frame","restart frame");
    parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
    parse_args.Add("-substeps",&write_substeps_level,"level","output-substep level");
    parse_args.Add("-substeps_delay",&substeps_delay_frame,"frame","delay substeps until after this frame");
    parse_args.Add("-dt",&dt,"dt","time step size to use for simulation");
    parse_args.Add("-steps",&time_steps_per_frame,"steps","number of time steps per frame");
    parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
    parse_args.Add("-dump_matrix",&dump_matrix,"dump out system and rhs");
    parse_args.Add("-sparse_dump_matrix",&sparse_dump_matrix,"dump out system (efficiently) and rhs");
    parse_args.Add("-use_pls",&use_pls_over_levelset,"use particle level set");
    parse_args.Add("-use_ls",&use_levelset_over_pls,"use level set method");
    parse_args.Add("-mu0",&mu0,&override_mu0,"viscosity","viscosity for first fluid region");
    parse_args.Add("-mu1",&mu1,&override_mu1,"viscosity","viscosity for second fluid region");
    parse_args.Add("-rho0",&rho0,&override_rho0,"density","density for first fluid region");
    parse_args.Add("-rho1",&rho1,&override_rho1,"density","density for second fluid region");
    parse_args.Add("-inv_Wi0",&inv_Wi0,&override_inv_Wi0,"inv_Wi","1/Wi for first fluid region");
    parse_args.Add("-inv_Wi1",&inv_Wi1,&override_inv_Wi1,"inv_Wi","1/Wi for second fluid region");
    parse_args.Add("-beta0",&beta0,&override_beta0,"beta","stress coefficient for first fluid region");
    parse_args.Add("-beta1",&beta1,&override_beta1,"beta","stress coefficient for second fluid region");
    parse_args.Add("-surface_tension",&surface_tension,&override_surface_tension,"density","density for second fluid region");
    parse_args.Add("-test_output_prefix",&test_output_prefix,&use_test_output,"","prefix to use for test output");
    parse_args.Add("-m",&m,"scale","meter scale");
    parse_args.Add("-s",&s,"scale","second scale");
    parse_args.Add("-kg",&kg,"scale","kilogram scale");
    parse_args.Add("-bc_n",&bc_n,"use Neumann boundary conditions");
    parse_args.Add("-bc_d",&bc_d,"use Dirichlet boundary conditions");
    parse_args.Add("-bc_s",&bc_s,"use slip boundary conditions");
    parse_args.Add("-test_diff",&test_analytic_diff,"test analytic derivatives");
    parse_args.Add_Not("-no_advect",&use_advection,"Disable advection");
    parse_args.Add("-no_solve",&omit_solve,"Disable visocity and pressure solve");
    parse_args.Add_Not("-no_reduced_advect",&use_reduced_advection,"Peform reduced advection");
    parse_args.Add("-refine",&refine,"num","Refine space/time by this factor");
    parse_args.Add("-null_p",&use_p_null_mode,"Assume pressure null mode and project it out");
    parse_args.Add("-threads",&number_of_threads,"threads","Number of threads");
    parse_args.Add("-o",&output_directory,&override_output_directory,"dir","Output directory");
    parse_args.Add("-dump_eigen",&dump_largest_eigenvector,"Dump largest few eigenvectors");
    parse_args.Add("-use_mg",&use_multigrid,"Use multigrid preconditioning");
    parse_args.Add("-mg_levels",&num_multigrid_levels,"levels","Number of multigrid levels");
    parse_args.Add("-max_iter",&max_iter,"number","Number of Krylov iterations");
    parse_args.Add("-solver_tolerance",&solver_tolerance,"tol","Krylov tolerance");
    parse_args.Add("-test_system",&test_system,"Run basic Krylov system tests");
    parse_args.Add_Not("-no_preconditioner",&use_preconditioner,"disable Jacobi preconditioner");
    parse_args.Add("-set_u0",&str_u0,&use_u0,"prog","Override velocity");
    parse_args.Add("-set_u1",&str_u1,&use_u1,"prog","Override velocity");
    parse_args.Add("-set_p0",&str_p0,&use_p0,"prog","Override pressure");
    parse_args.Add("-set_p1",&str_p1,&use_p1,"prog","Override pressure");
    parse_args.Add("-set_S0",&str_S0,&use_S0,"prog","Override polymer stress");
    parse_args.Add("-set_S1",&str_S1,&use_S1,"prog","Override polymer stress");
    parse_args.Parse(true);

#ifdef USE_OPENMP
    omp_set_num_threads(number_of_threads);
#pragma omp parallel
#pragma omp single
    {
        if(omp_get_num_threads()!=number_of_threads) PHYSBAM_FATAL_ERROR();
        LOG::cout<<"Running on "<<number_of_threads<<" threads"<<std::endl;
    }
#endif

    resolution*=refine;
    dt/=refine;
    time_steps_per_frame*=refine;
    stored_last_frame=last_frame;
    if(weiss>1e9)weiss_inv=0;else weiss_inv = 1.0/weiss;
    unit_mu=kg*pow<2-TV::m>(m)/s;
    unit_rho=kg/pow<TV::m>(m);
    unit_st=kg*pow<3-TV::m>(m)/(s*s);
    unit_p=kg*pow<2-TV::m>(m)/(s*s);
    mu0*=unit_mu;
    mu1*=unit_mu;
    rho0*=unit_rho;
    rho1*=unit_rho;
    inv_Wi0/=s;
    inv_Wi1/=s;
    beta0*=unit_p;
    beta1*=unit_p;
    surface_tension*=unit_st;
    dt*=s;
    PHYSBAM_ASSERT(bc_n+bc_d+bc_s<2);
    bc_type=bc_n?NEUMANN:(bc_s?SLIP:DIRICHLET);

    analytic_levelset=0;
    save_pressure=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUIDS_COLOR_BASE<TV>::
~FLUIDS_COLOR_BASE()
{
    delete analytic_levelset;
    analytic_velocity.Delete_Pointers_And_Clean_Memory();
    analytic_pressure.Delete_Pointers_And_Clean_Memory();
    analytic_polymer_stress.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function After_Initialize_Example
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
After_Initialize_Example()
{
    if(!override_output_directory) output_directory=LOG::sprintf("Test_%d",test_number);

    if(use_pls_over_levelset){if(use_level_set_method) use_pls=true;use_level_set_method=false;}
    if(use_levelset_over_pls){if(use_pls) use_level_set_method=true;use_pls=false;}

    if(analytic_velocity.m) number_of_colors=analytic_velocity.m;

    if(use_u0 && analytic_velocity.m>=1){
        delete analytic_velocity(0);
        analytic_velocity(0)=new ANALYTIC_VELOCITY_PROGRAM<TV>(str_u0);}

    if(use_u1 && analytic_velocity.m>=2){
        delete analytic_velocity(1);
        analytic_velocity(1)=new ANALYTIC_VELOCITY_PROGRAM<TV>(str_u1);}

    if(use_p0 && analytic_pressure.m>=1){
        delete analytic_pressure(0);
        analytic_pressure(0)=new ANALYTIC_PRESSURE_PROGRAM<TV>(str_p0);}

    if(use_p1 && analytic_pressure.m>=2){
        delete analytic_pressure(1);
        analytic_pressure(1)=new ANALYTIC_PRESSURE_PROGRAM<TV>(str_p1);}

    if(use_S0 && analytic_polymer_stress.m>=1){
        delete analytic_polymer_stress(0);
        analytic_polymer_stress(0)=new ANALYTIC_POLYMER_STRESS_PROGRAM<TV>(str_S0);}

    if(use_S1 && analytic_polymer_stress.m>=2){
        delete analytic_polymer_stress(1);
        analytic_polymer_stress(1)=new ANALYTIC_POLYMER_STRESS_PROGRAM<TV>(str_S1);}
}
//#####################################################################
// Function Initialize_Common_Example
//#####################################################################
template<class TV> bool FLUIDS_COLOR_BASE<TV>::
Initialize_Common_Example()
{
    typename TV::SPIN spin_count;
    TV vector_count;
    for(int i=0;i<TV::SPIN::m;i++) spin_count(i)=i;
    for(int i=0;i<TV::m;i++) vector_count(i)=i;
    MATRIX<T,TV::m> spin_mat=MATRIX<T,TV::m>::Cross_Product_Matrix(spin_count+1);
    IDENTITY_MATRIX<T,TV::m> id;

    //Tests 0-31: Tests created prior to selection of examples for 2013 JCP paper
    //Tests 101-111: Examples 1-11 in 2013 JCP paper, may coincide with Tests 0-31
    //Tests 201-249: Examples for upcoming polymer stress JCP paper
    //Tests 250-299: Tests using polymer stress, prior to example selection
    switch(test_number){
        case 0:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            Add_Velocity([](auto X,auto t){return TV()+1;});
            Add_Pressure([](auto X,auto t){return 0;});
            use_p_null_mode=true;
            break;
        case 3:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
            auto rot=[=](auto X,auto t){return spin_mat*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(rho0/unit_rho)*rot(X,t).Magnitude_Squared();});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            break;}
        case 4:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
            Add_Velocity([](auto X,auto t){return TV()+1;});
            Add_Pressure([](auto X,auto t){return 0;});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            break;
        case 7:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            Add_Velocity([](auto X,auto t){return X.Dot(TV::Axis_Vector(0))/(t+1)*TV::Axis_Vector(0);});
            Add_Pressure([](auto X,auto t){return 0;});
            omit_solve=true;
            break;
        case 9:
        case 106:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            {
                T x0=(T).2,x1=(T).5,x2=(T).8,v0=1,v2=-1,u_mu0=mu0/unit_mu,u_mu1=mu1/unit_mu;
                T v1=(v0*u_mu0/(x1-x0)+v2*u_mu1/(x2-x1))/(u_mu0/(x1-x0)+u_mu1/(x2-x1));
                MATRIX<T,TV::m> du0,du1;
                du0(1,0)=(v1-v0)/(x1-x0);
                du1(1,0)=(v2-v1)/(x2-x1);
                ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x0,TV::Axis_Vector(0),DIRICHLET,0);
                ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x2,TV::Axis_Vector(0),1,DIRICHLET);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
                Add_Velocity([=](auto X,auto t){return du0*(X-TV::Axis_Vector(0)*x1)+TV::Axis_Vector(1)*v1;});
                Add_Velocity([=](auto X,auto t){return du1*(X-TV::Axis_Vector(0)*x1)+TV::Axis_Vector(1)*v1;});
                Add_Pressure([](auto X,auto t){return 0;});
                Add_Pressure([](auto X,auto t){return 0;});
                use_p_null_mode=true;
//                    if(test_number==106) use_level_set_method=true;
            }
            break;
        case 11:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
            auto rot=[=](auto X,auto t){return spin_mat*(X-(vector_count*(T).2+(T).6));};
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(rho0/unit_rho)*rot(X,t).Magnitude_Squared();});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            break;}
        case 12:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            {
                T x0=(T).2,x1=(T).8,g=9.8,a=rho0/unit_rho*g/(2*mu0/unit_mu);
                ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x0,TV::Axis_Vector(0),DIRICHLET,0);
                ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),DIRICHLET,DIRICHLET);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);

                Add_Velocity([=](auto X,auto t){auto x=X.Dot(TV::Axis_Vector(0));return TV::Axis_Vector(1)*a*((x-(x0+x1))*x+x0*x1);});
                Add_Pressure([](auto X,auto t){return 0;});
                use_p_null_mode=true;
            }
            break;
        case 13:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            {
                T x0=(T).2,x1=(T).5,x2=(T).8,v0=1,v1=1,v2=-1;
                MATRIX<T,TV::m> du0,du1;
                du0(1,0)=(v1-v0)/(x1-x0);
                du1(1,0)=(v2-v1)/(x2-x1);
                ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x0,TV::Axis_Vector(0),DIRICHLET,0);
                ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x2,TV::Axis_Vector(0),1,DIRICHLET);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
                Add_Velocity([=](auto X,auto t){return du0*(X-TV::Axis_Vector(0)*x1)+TV::Axis_Vector(1)*v1;});
                Add_Velocity([=](auto X,auto t){return du1*(X-TV::Axis_Vector(0)*x1)+TV::Axis_Vector(1)*v1;});
                Add_Pressure([](auto X,auto t){return 0;});
                Add_Pressure([](auto X,auto t){return 0;});
                use_p_null_mode=true;
            }
            break;
        case 14:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            {
                T x0=(T).2,x1=(T).5,x2=(T).8,v0=1,v2=-1;
                ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x0,TV::Axis_Vector(0),DIRICHLET,0);
                ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x2,TV::Axis_Vector(0),1,DIRICHLET);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
                Add_Velocity([=](auto X,auto t){return TV::Axis_Vector(1)*v0;});
                Add_Velocity([=](auto X,auto t){return TV::Axis_Vector(1)*v2;});
                Add_Pressure([](auto X,auto t){return 0;});
                Add_Pressure([](auto X,auto t){return 0;});
                use_discontinuous_velocity=true;
                use_p_null_mode=true;
            }
            break;
        case 15:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            {
                T x0=(T).2,x1=(T).5,x2=(T).8;
                ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),x0,DIRICHLET,0);
                ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),x2,1,DIRICHLET);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),x1,0,1)))->Add(ab)->Add(cd);
                Add_Velocity([](auto X,auto t){auto x=X.Dot(TV::Axis_Vector(0));return TV::Axis_Vector(1)*(((T).5*x+(T).2)*x+(T).3);});
                auto rot=[=](auto X,auto t){return spin_mat*(X-(TV()+(T).1));};
                Add_Velocity(rot);
                Add_Pressure([](auto X,auto t){return 0;});
                Add_Pressure([=](auto X,auto t){return (T).5*(rho1/unit_rho)*rot(X,t).Magnitude_Squared();});
                use_discontinuous_velocity=true;
                use_p_null_mode=true;
            }
            break;
        case 17:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            {
                TV vel=vector_count*(T).3+(T).2;
                analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4),vel);
                Add_Velocity([=](auto X,auto t){return vel+1;});
                Add_Pressure([](auto X,auto t){return 0;});
                if(bc_type!=NEUMANN) use_p_null_mode=true;
            }
            break;
        case 18:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            {
                analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4),vector_count*(T).1+(T).1);
                TV trans_vel=vector_count*(T)-.9+(T).5;

                auto rot=[=](auto X,auto t){return spin_mat*(X-trans_vel*t-(TV()+1));};
                Add_Velocity([=](auto X,auto t){return rot(X,t)+trans_vel;});
                Add_Pressure([rho0=rho0,unit_rho=unit_rho,rot=rot](auto X,auto t){return (T).5*(rho0/unit_rho)*rot(X,t).Magnitude_Squared();});
                if(bc_type!=NEUMANN) use_p_null_mode=true;
            }
            break;
        case 20:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            {
                TV vel=vector_count*(T)-.7+(T).2;
                analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,1,0),vel);
                Add_Velocity([=](auto X,auto t){return vel;});
                Add_Velocity([=](auto X,auto t){return vel;});
                Add_Pressure([](auto X,auto t){return 0;});
                Add_Pressure([](auto X,auto t){return 0;});
                use_p_null_mode=true;
                use_level_set_method=true;
            }
            break;
        case 21:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            {
                TV vel=TV::Axis_Vector(0);
                T x0=(T).2,x1=(T).5,x2=(T).8;
                ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x0,TV::Axis_Vector(0),DIRICHLET,0);
                ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x2,TV::Axis_Vector(0),1,DIRICHLET);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(0),0,1),vel)))->Add(ab)->Add(cd);
                Add_Velocity([=](auto X,auto t){return vel;});
                Add_Velocity([=](auto X,auto t){return vel;});
                Add_Pressure([](auto X,auto t){return 0;});
                Add_Pressure([](auto X,auto t){return 0;});
                use_p_null_mode=true;
                use_level_set_method=true;
            }
            break;
        case 22:
        case 107:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m/100,true);
            {
                T radius=(T).003,curvature=(TV::m-1)/radius;
                if(!override_surface_tension) surface_tension=(T)0.07197*unit_st;
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).005,radius,0,1);
                Add_Velocity([](auto X,auto t){return TV();});
                Add_Pressure([=](auto X,auto t){return (surface_tension*curvature)/unit_p;});
                Add_Velocity([](auto X,auto t){return TV();});
                Add_Pressure([](auto X,auto t){return 0;});
                use_p_null_mode=true;
                use_level_set_method=true;
            }
            break;
        case 24:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
            Add_Velocity([](auto X,auto t){return TV();});
            Add_Pressure([](auto X,auto t){return 1;});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            break;}
        case 25:{//A circle being rotated around inside a larger circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);

            ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(vector_count*(T)-.2+(T).5,.1,1,0),spin_count+1,TV()+(T).5);
            ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),-4,-4);
            analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).4,0,1)))->Add(ab)->Add(cd);
            auto rot=[=](auto X,auto t){return spin_mat*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(rho0/unit_rho)*rot(X,t).Magnitude_Squared();});
            Add_Pressure([rho1=rho1,unit_rho=unit_rho,rot=rot](auto X,auto t){return (T).5*(rho1/unit_rho)*rot(X,t).Magnitude_Squared();});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            //use_level_set_method=true;

            break;}
        case 26:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            TV vel=vector_count*(T).3+(T).2;
            analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_BOX<TV>(TV()+(T).2,TV()+(T).45,1,0),spin_count*0,TV()+(T).5),vel);
            Add_Velocity([=](auto X,auto t){return vel;});
            Add_Velocity([=](auto X,auto t){return vel;});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Pressure([](auto X,auto t){return 0;});
            use_p_null_mode=true;
            use_level_set_method=true;
            break;}
        case 28:{//Test 25 but using the level set method
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(vector_count*(T)-.0+(T).5,.1,1,0),spin_count+1,TV()+(T).5);
            ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),-4,-4);
            analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).4,0,1)))->Add(ab)->Add(cd);
            auto rot=[=](auto X,auto t){return spin_mat*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(rho0/unit_rho)*rot(X,t).Magnitude_Squared();});
            Add_Pressure([rho1=rho1,unit_rho=unit_rho,rot=rot](auto X,auto t){return (T).5*(rho1/unit_rho)*rot(X,t).Magnitude_Squared();});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_level_set_method=true;
                
            break;}
        case 31:{//Test to illustrate issues with extrapolation
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            TV vel=TV::Axis_Vector(0);
            T x0=(T).2,x1=(T).8;
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x0,TV::Axis_Vector(0),DIRICHLET,0);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),DIRICHLET,DIRICHLET);
            analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
            Add_Velocity([=](auto X,auto t){return vel;});
            Add_Pressure([](auto X,auto t){return 0;});
            use_p_null_mode=true;
            break;
        }
        case 104:
        case 105:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*(T)pi*m,true);
            ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),.8*(T)pi,1,SLIP);
            analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV::Axis_Vector(0)*.2*pi,(T).2*(T)pi,0,1)))->Add(ab)->Add(cd);
            MATRIX<T,TV::m> du0;for(int i=0;i<TV::m;i++)du0(i,i)=-1;du0(1,1)+=TV::m;
            Add_Velocity([=](auto X,auto t){return du0*(X-TV::Axis_Vector(0)*.2*pi);});
            auto rot=[=](auto X,auto t){return spin_mat*X;};
            Add_Velocity(rot);
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Pressure([=](auto X,auto t){return (T).5*(rho1/unit_rho)*rot(X,t).Magnitude_Squared();});
            use_discontinuous_velocity=true;
            use_p_null_mode=true;
            break;
        }
        case 108:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            TV r=TV()+(T).4;
            r(0)=(T).7;
            if(!override_mu0) mu0=3*unit_mu;
            if(!override_mu1) mu1=1*unit_mu;
            if(!override_rho0) rho0=0.001*unit_rho;
            if(!override_rho1) rho1=0.002*unit_rho;
            if(!override_surface_tension) surface_tension=(T)10*unit_st;
            analytic_levelset=new ANALYTIC_LEVELSET_ELLIPSOID<TV>(TV(),r,0,1);
            Add_Velocity([](auto X,auto t){return TV();});
            Add_Velocity([](auto X,auto t){return TV();});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Pressure([](auto X,auto t){return 0;});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_pls=true;
            analytic_initial_only=true;
            break;
        }
        case 250:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            Add_Velocity([](auto X,auto t){return TV()+1;});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Polymer_Stress([](auto X,auto t){return SM()+1;});
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }
        case 251:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
            Add_Velocity([](auto X,auto t){return TV()+1;});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Polymer_Stress([](auto X,auto t){return SM()+1;});
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }
        case 252:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);    
            T x0=(T).2,x1=(T).5,x2=(T).8;TV a;a(0)=rho0/unit_rho;
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x0,TV::Axis_Vector(0),-4,0);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x2,TV::Axis_Vector(0),0,-4);
            analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
            Add_Velocity([](auto X,auto t){return TV()+1;});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Polymer_Stress([=](auto X,auto t){return id*(X.Dot(a)+(T)1);});
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }
        case 255:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
            auto rot=[=](auto X,auto t){return spin_mat*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(rho0/unit_rho)*rot(X,t).Magnitude_Squared();});
            Add_Polymer_Stress([=](auto X,auto t){return id*(exp(-weiss_inv*t)+(T)1);});
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }
        case 258:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            Add_Velocity([](auto X,auto t){return TV()+1;});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Polymer_Stress(
                [=](auto X,auto t)
                {
                    auto Z=(T)(2*pi)*X;
                    auto W=(cos(Z)+sin(Z))*((T)1+exp(-t));
                    return Outer_Product(W)+id;
                });
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }
        case 259:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            Add_Velocity([](auto X,auto t){return sin((T)(2*pi)*X*3);});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Polymer_Stress(
                [=](auto X,auto t)
                {
                    auto Z=(T)(2*pi)*X;
                    auto W=(cos(Z)+sin(Z))*((T)1+exp(-t));
                    return Outer_Product(W)+id;
                });
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }

        default: return false;}
    return true;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Write_Output_Files(const int frame)
{
    BASE::Write_Output_Files(frame);

    if(use_test_output){
        std::string file=LOG::sprintf("%s/%s-%03d.txt",output_directory.c_str(),test_output_prefix.c_str(),frame);
        OCTAVE_OUTPUT<T> oo(file.c_str());
        ARRAY<T> u;
        for(int c=0;c<face_velocities.m;c++){
            ARRAY_VIEW<T> a(face_velocities(c).buffer_size,face_velocities(c).base_pointer);
            u.Append_Elements(a);}
        oo.Write("u",u);}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Initialize()
{
    if(analytic_levelset && analytic_velocity.m) Analytic_Test();
    else
        switch(test_number){
            default: PHYSBAM_FATAL_ERROR("Missing test number");}

    PHYSBAM_ASSERT(rho.m==number_of_colors && mu.m==number_of_colors);
    if(user_last_frame) last_frame=stored_last_frame;
}
//#####################################################################
// Function Set_Level_Set
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Set_Level_Set(T time)
{
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        int c=-4;
        T p=analytic_levelset->phi(it.Location()/m,time/s,c)*m;
        levelset_color.phi(it.index)=abs(p);
        levelset_color.color(it.index)=c==-4?bc_type:c;}
    Fill_Levelsets_From_Levelset_Color();
}
//#####################################################################
// Function Level_Set_Error
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Level_Set_Error(T time)
{
    if(!analytic_levelset || analytic_initial_only) return;
    ARRAY<VECTOR<T,3> > colors;
    colors.Append(VECTOR<T,3>((T).25,(T).25,(T).25));
    colors.Append(VECTOR<T,3>((T).5,(T).5,(T).5));
    colors.Append(VECTOR<T,3>(1,1,1));
    colors.Append(VECTOR<T,3>(1,0,0));
    colors.Append(VECTOR<T,3>(0,1,0));
    colors.Append(VECTOR<T,3>(0,0,1));
    colors.Append(VECTOR<T,3>(1,1,0));
    colors.Append(VECTOR<T,3>(0,1,1));
    colors.Append(VECTOR<T,3>(1,0,1));

    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        int c=levelset_color.color(it.index);
        T p=levelset_color.phi(it.index);
        Add_Debug_Particle(it.Location(),colors(c+3));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("level set",0,1);
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        int c=-4;
        T p=analytic_levelset->phi(it.Location()/m,time/s,c)*m;
        if(levelset_color.color(it.index)!=(c==-4?bc_type:c)){
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(levelset_color.phi(it.index))+abs(p));}
        else{
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,0));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(levelset_color.phi(it.index)-p));}}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("level set error",0,1);
}
//#####################################################################
// Function Velocity_Error
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Velocity_Error(T time)
{
    if(!analytic_velocity.m || analytic_initial_only) return;
    T u_inf=0,u_2=0,p_inf=0,p_2=0,p_ave=0,a=0,b=0,pa=0,pb=0,l_inf=0,l_2=0;
    int num_u=0,num_p=0,num_l=0;
        
    ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > > saved_face_velocities(face_velocities);
    face_velocities*=(T)0;
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        int c=levelset_color.Color(it.Location());
        if(c<0) continue;
        T A=saved_face_velocities(c)(it.Full_Index()),B=analytic_velocity(c)->u(it.Location()/m,time/s)(it.Axis())*m/s;
        a=max(a,abs(A));
        b=max(b,abs(B));
        u_inf=max(u_inf,abs(A-B));
        u_2+=sqr(A-B);
        num_u++;
        face_velocities(c)(it.Full_Index())=A-B;}
    if(num_u) u_2/=num_u;
    u_2=sqrt(u_2);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("velocity error",0,1);
    face_velocities=saved_face_velocities;
    if(use_p_null_mode){
        int cnt=0;
        for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
//                if(levelset_color.Phi(it.Location())<grid.dX.Max()/2) continue;
            int c=levelset_color.Color(it.Location());
            if(c<0) continue;
            T A=pressure(it.index),B=analytic_pressure(c)->p(it.Location()/m,time/s)*unit_p,D=A-B;
            p_ave+=D;
            cnt++;}
        if(cnt) p_ave/=cnt;}
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
//            if(levelset_color.Phi(it.Location())<grid.dX.Max()/2) continue;
        int c=levelset_color.Color(it.Location());
        if(c<0) continue;
        T A=pressure(it.index),B=analytic_pressure(c)->p(it.Location()/m,time/s)*unit_p,D=A-B-p_ave;
        pa=max(pa,abs(A));
        pb=max(pb,abs(B));
        p_inf=max(p_inf,abs(D));
        p_2+=sqr(D);
        num_p++;
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(D));}
    if(num_p) p_2/=num_p;
    p_2=sqrt(p_2);

    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        int c=-4;
        T p=analytic_levelset->phi(it.Location()/m,time/s,c)*m;
        if(abs(p)>grid.dX.Max()*3) continue;
        T e=levelset_color.color(it.index)!=(c==-4?bc_type:c)?abs(levelset_color.phi(it.index))+abs(p):abs(levelset_color.phi(it.index)-p);
        l_inf=max(l_inf,e);
        l_2+=sqr(e);
        num_l++;}
    if(num_l) l_2/=num_l;
    l_2=sqrt(l_2);

    LOG::printf("max_error %-22.16g %-22.16g %-22.16g %-22.16g  p %-22.16g %-22.16g %-22.16g %-22.16g  l %-22.16g %-22.16g", u_inf, u_2, a, b, p_inf, p_2, pa, pb, l_inf, l_2);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("pressure error",0,1);

    if(use_polymer_stress){
        T max_S=0;
        for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
            int c=levelset_color.Color(it.Location());
            if(c<0) continue;
            SYMMETRIC_MATRIX<T,TV::m> A=polymer_stress(c)(it.index),B=analytic_polymer_stress(c)->S(it.Location()/m,time/s),D=A-B;
            T norm=D.Frobenius_Norm();
            max_S=max(max_S,norm);
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,norm);}
        LOG::printf("S error: %-22.16g",max_S);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("polymer stress error",0,1);}
    LOG::printf("\n");
}
//#####################################################################
// Function Dump_Analytic_Levelset
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Dump_Analytic_Levelset(T time)
{
    ARRAY<VECTOR<T,3> > colors;
    colors.Append(VECTOR<T,3>((T).25,(T).25,(T).25));
    colors.Append(VECTOR<T,3>((T).5,(T).5,(T).5));
    colors.Append(VECTOR<T,3>(1,1,1));
    colors.Append(VECTOR<T,3>(1,0,0));
    colors.Append(VECTOR<T,3>(0,1,0));
    colors.Append(VECTOR<T,3>(0,0,1));
    colors.Append(VECTOR<T,3>(1,1,0));
    colors.Append(VECTOR<T,3>(0,1,1));
    colors.Append(VECTOR<T,3>(1,0,1));
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        int c=-4;
        T p=analytic_levelset->phi(it.Location()/m,time/s,c)*m;
        if(c==-4) c=bc_type;
        Add_Debug_Particle(it.Location(),colors(c+3));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("analytic level set (phi)",0,1);
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        int c=-4;
        T p=analytic_levelset->phi(it.Location()/m,time/s,c)*m;
        (void)p;
        if(c==-4) c=bc_type;
        Add_Debug_Particle(it.Location(),colors(c+3));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,analytic_levelset->N(it.Location()/m,time/s,c));}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("analytic level set (N)",0,1);
}
//#####################################################################
// Function Get_Initial_Velocities
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Get_Initial_Velocities(T time)
{
    if(analytic_levelset && analytic_velocity.m)
        for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
            int c=face_color(it.Full_Index());
            if(c<0) continue;
            face_velocities(c)(it.Full_Index())=analytic_velocity(c)->u(it.Location()/m,time/s)(it.Axis())*m/s;}
}
//#####################################################################
// Function Get_Initial_Polymer_Stresses
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Get_Initial_Polymer_Stresses(T time)
{
    //Right now we are filling in for all colors. We may not want to do that in the future.
    if(analytic_levelset && analytic_velocity.m)
        for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
            for(int c=0;c<analytic_velocity.m;c++)
                polymer_stress(c)(it.index)=analytic_polymer_stress(c)->S(it.Location()/m,time/s);}
}
//#####################################################################
// Function Analytic_Test
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Analytic_Test()
{
    mu.Append(mu0);
    rho.Append(rho0);
    inv_Wi.Append(inv_Wi0);
    polymer_stress_coefficient.Append(beta0);
    for(int i=1;i<number_of_colors;i++){
        mu.Append(mu1);
        rho.Append(rho1);
        inv_Wi.Append(inv_Wi1);
        polymer_stress_coefficient.Append(beta1);}

    Set_Level_Set(0);
    Dump_Analytic_Levelset(0);

    if(test_analytic_diff){
        analytic_levelset->Test(grid.domain);
        RANDOM_NUMBERS<T> rand;
        TV X;
        int c=-4;
        do{
            X=rand.Get_Uniform_Vector(grid.domain);
            analytic_levelset->phi(X/m,0,c);}
        while(c<0);
        analytic_velocity(c)->Test(X);
        analytic_polymer_stress(c)->Test(X);}
}
//#####################################################################
// Function Begin_Time_Step
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
Begin_Time_Step(const T time)
{
    if(analytic_velocity.m && analytic_levelset && !use_level_set_method && !use_pls && !analytic_initial_only)
        Set_Level_Set(time+dt);
}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class TV> void FLUIDS_COLOR_BASE<TV>::
End_Time_Step(const T time)
{
    Level_Set_Error(time);
    Velocity_Error(time);

    T max_u=0;
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        int c=levelset_color.Color(it.Location());
        if(c<0) continue;
        max_u=std::max(max_u,abs(face_velocities(c)(it.Full_Index())));}
    LOG::cout<<"max u "<<max_u<<std::endl;
}
//#####################################################################
// Function Stress
//#####################################################################
template<class TV> auto FLUIDS_COLOR_BASE<TV>::
Stress(const TV& X,int color,T time) -> MATRIX<T,TV::m>
{
    T p=analytic_pressure(color)->p(X/m,time/s)*unit_p;
    MATRIX<T,TV::m> du=analytic_velocity(color)->du(X/m,time/s)/s;
    return (du+du.Transposed())*mu(color)-p;
}
//#####################################################################
// Function Polymer_Stress
//#####################################################################
template<class TV> auto FLUIDS_COLOR_BASE<TV>::
Polymer_Stress(const TV& X,int color,T time) -> SYMMETRIC_MATRIX<T,TV::m>
{
    PHYSBAM_ASSERT(use_polymer_stress);
    return analytic_polymer_stress(color)->S(X/m,time/s);
}
//#####################################################################
// Function Polymer_Stress_Forcing_Term
//#####################################################################
template<class TV> auto FLUIDS_COLOR_BASE<TV>::
Polymer_Stress_Forcing_Term(const TV& X,int color,T time) -> SYMMETRIC_MATRIX<T,TV::m>
{
    PHYSBAM_ASSERT(use_polymer_stress);
    TV Z=X/m;
    time/=s;
    ANALYTIC_POLYMER_STRESS<TV>& as=*analytic_polymer_stress(color);
    ANALYTIC_VELOCITY<TV>& au=*analytic_velocity(color);
    TV u=au.u(Z,time);
    MATRIX<T,TV::m> du=au.du(Z,time);
    SYMMETRIC_MATRIX<T,TV::m> S=as.S(Z,time);
    SYMMETRIC_MATRIX<T,TV::m> f=as.dSdt(Z,time);
    if(use_advection) f+=Contract<2>(as.dS(Z,time),u);
    f-=(du*S).Twice_Symmetric_Part()+inv_Wi(color)*s*((T)1-S);
    return f/s;
}
//#####################################################################
// Function Jump_Interface_Condition
//#####################################################################
template<class TV> TV FLUIDS_COLOR_BASE<TV>::
Jump_Interface_Condition(const TV& X,int color0,int color1,T time)
{
    Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));

    if(surface_tension){
        LEVELSET<TV>& ls=*particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.levelsets(color1);
        static T last_time=-1;
        if(last_time!=time){
            last_time=time;
            ls.Compute_Curvature(time);
            ls.Compute_Normals(time);}

        T k=CUBIC_MN_INTERPOLATION_UNIFORM<TV,T>().Clamped_To_Array(grid,*ls.curvature,X);
        TV n=CUBIC_MN_INTERPOLATION_UNIFORM<TV,TV>().Clamped_To_Array(grid,*ls.normals,X).Normalized();
        Add_Debug_Particle(X,VECTOR<T,3>(0,0,1));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,k*surface_tension*n);
        return k*surface_tension*n;}

    if(analytic_velocity.m && analytic_levelset && !analytic_initial_only){
        MATRIX<T,TV::m> jump_stress=Stress(X,color1,time); if(use_polymer_stress) jump_stress+=polymer_stress_coefficient(color1)*Polymer_Stress(X,color1,time);
        if(color0>=0){ jump_stress-=Stress(X,color0,time); if(use_polymer_stress) jump_stress-=polymer_stress_coefficient(color0)*Polymer_Stress(X,color0,time);}
        TV n=analytic_levelset->N(X/m,time/s,color1);
        return jump_stress*n;}

    return TV();
}
//#####################################################################
// Function Volume_Force
//#####################################################################
template<class TV> TV FLUIDS_COLOR_BASE<TV>::
Volume_Force(const TV& X,int color,T time)
{
    if(analytic_velocity.m && analytic_levelset && !analytic_initial_only){
        ANALYTIC_VELOCITY<TV>* av=analytic_velocity(color);
        ANALYTIC_PRESSURE<TV>* ap=analytic_pressure(color);
        T rh=rho(color)/unit_rho;
        T mh=mu(color)/unit_mu;
        TV f=rh*av->dudt(X/m,time/s)+ap->dp(X/m,time/s)-mh*av->Lu(X/m,time/s);
        if(use_advection) f+=rh*av->du(X/m,time/s)*av->u(X/m,time/s);
        if(use_polymer_stress){
            T beta=polymer_stress_coefficient(color)/unit_p;
            f-=beta*Contract<1,2>(analytic_polymer_stress(color)->dS(X/m,time/s));}
        return f*unit_p/m;}
    return gravity;
}
//#####################################################################
// Function Velocity_Jump
//#####################################################################
template<class TV> TV FLUIDS_COLOR_BASE<TV>::
Velocity_Jump(const TV& X,int color0,int color1,T time)
{
    Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
    if(analytic_velocity.m && analytic_levelset && !analytic_initial_only){
        TV jump_u=analytic_velocity(color1)->u(X/m,time/s);
        if(color0>=0) jump_u-=analytic_velocity(color0)->u(X/m,time/s);
        return jump_u*m/s;}
    return TV();
}
//#####################################################################
}
