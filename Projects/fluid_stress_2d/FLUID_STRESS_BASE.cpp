//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_BOX.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ELLIPSOID.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ROTATE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SCALE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SPHERE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_TRANSLATE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_VORTEX.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include "ANALYTIC_POLYMER_STRESS.h"
#include "ANALYTIC_VELOCITY.h"
#include "FLUID_STRESS_BASE.h"
#include "STRESS_EXAMPLE.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_STRESS_BASE<TV>::
FLUID_STRESS_BASE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
    :STRESS_EXAMPLE<TV>(stream_type),test_number(0),resolution(32),stored_last_frame(0),user_last_frame(false),
    unit_mu(0),unit_rho(0),unit_st(0),unit_p(0),weiss(1),weiss_inv(1),m(1),s(1),kg(1),
    bc_n(false),bc_d(false),test_analytic_diff(false),refine(1),analytic_initial_only(false),
    number_of_threads(1),override_output_directory(false)
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
    parse_args.Add("-m",&m,"scale","meter scale");
    parse_args.Add("-s",&s,"scale","second scale");
    parse_args.Add("-kg",&kg,"scale","kilogram scale");
    parse_args.Add("-bc_n",&bc_n,"use Neumann boundary conditions");
    parse_args.Add("-bc_d",&bc_d,"use Dirichlet boundary conditions");
    parse_args.Add("-bc_s",&bc_s,"use slip boundary conditions");
    parse_args.Add("-test_diff",&test_analytic_diff,"test analytic derivatives");
    parse_args.Add("-refine",&refine,"num","Refine space/time by this factor");
    parse_args.Add("-threads",&number_of_threads,"threads","Number of threads");
    parse_args.Add("-o",&output_directory,&override_output_directory,"dir","Output directory");
    parse_args.Parse(true);

    resolution*=refine;
    dt/=refine;
    time_steps_per_frame*=refine;
    stored_last_frame=last_frame;
    if(weiss>1e9)weiss_inv=0;else weiss_inv = 1.0/weiss;
    unit_mu=kg*pow<2-TV::m>(m)/s;
    unit_rho=kg/pow<TV::m>(m);
    unit_st=kg*pow<3-TV::m>(m)/(s*s);
    unit_p=kg*pow<2-TV::m>(m)/(s*s);
    dt*=s;
    PHYSBAM_ASSERT(bc_n+bc_d+bc_s<2);
    bc_type=bc_n?NEUMANN:(bc_s?SLIP:DIRICHLET);
    
    analytic_levelset=0;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_STRESS_BASE<TV>::
~FLUID_STRESS_BASE()
{
    delete analytic_levelset;
    delete analytic_velocity;
    delete analytic_polymer_stress;
}
//#####################################################################
// Function After_Initialize_Example
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
After_Initialize_Example()
{
    if(!override_output_directory) output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number);
}
//#####################################################################
// Function Initialize_Common_Example
//#####################################################################
template<class TV> bool FLUID_STRESS_BASE<TV>::
Initialize_Common_Example()
{
    typename TV::SPIN spin_count;
    typedef SYMMETRIC_MATRIX<T,TV::m> T_MAT;
    TV vector_count;
    for(int i=0;i<TV::SPIN::m;i++) spin_count(i)=i;
    for(int i=0;i<TV::m;i++) vector_count(i)=i;

    switch(test_number){
        case 0:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            analytic_velocity=new ANALYTIC_VELOCITY_CONST<TV>(TV()+1);
            analytic_polymer_stress=new ANALYTIC_POLYMER_STRESS_CONST<TV>();
            break;
        case 1:
            struct ANALYTIC_POLYMER_1:public ANALYTIC_POLYMER_STRESS<TV>
            {
                ANALYTIC_POLYMER_1(){}
                ~ANALYTIC_POLYMER_1() {}
                virtual T_MAT S(const TV& X,T t) const {return T_MAT::Outer_Product(sin(X));}
                virtual T_MAT dSdX(const TV& X,T t,int a) const
                {return T_MAT::Symmetric_Outer_Product(sin(X),TV::Axis_Vector(a)*cos(X(a)));}
                virtual T_MAT dSdt(const TV& X,T t) const {return T_MAT();}
            };
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*pi*m),true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            analytic_velocity=new ANALYTIC_VELOCITY_CONST<TV>(TV(.3,.3));
            analytic_polymer_stress=new ANALYTIC_POLYMER_STRESS_TRANSLATE<TV>(new ANALYTIC_POLYMER_1,TV(.3,.3));
            break;
        case 2:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),.8,1,0);
            analytic_velocity=new ANALYTIC_VELOCITY_ROTATION<TV>(TV(),typename TV::SPIN()+1,1);
            analytic_polymer_stress=new ANALYTIC_POLYMER_STRESS_ROTATION<TV>(new ANALYTIC_POLYMER_1,TV(),typename TV::SPIN()+1);
            break;
        case 3:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),.8,1,0);
            analytic_velocity=new ANALYTIC_VELOCITY_ROTATION<TV>(TV(),typename TV::SPIN()+1,1);
            analytic_polymer_stress=new ANALYTIC_POLYMER_1;
            break;
        case 4:
            struct ANALYTIC_POLYMER_4:public ANALYTIC_POLYMER_STRESS<TV>
            {
                ANALYTIC_POLYMER_4(){}
                ~ANALYTIC_POLYMER_4() {}
                virtual T_MAT S(const TV& X,T t) const {return T_MAT::Outer_Product(sin(X))*(2+sin(t));}
                virtual T_MAT dSdX(const TV& X,T t,int a) const
                {return T_MAT::Symmetric_Outer_Product(sin(X),TV::Axis_Vector(a)*cos(X(a)))*(2+sin(t));}
                virtual T_MAT dSdt(const TV& X,T t) const {return T_MAT::Outer_Product(sin(X))*cos(t);}
            };
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),.8,1,0);
            analytic_velocity=new ANALYTIC_VELOCITY_ROTATION<TV>(TV(),typename TV::SPIN()+1,1);
            analytic_polymer_stress=new ANALYTIC_POLYMER_4;
            break;
        case 5:
            struct ANALYTIC_VELOCITY_5:public ANALYTIC_VELOCITY<TV>
            {
                typedef typename TV::SCALAR T;
                TV c;
                typename TV::SPIN w;
                ANALYTIC_VELOCITY_5(TV cc,typename TV::SPIN ww,T rho): c(cc),w(ww){}
                virtual TV u(const TV& X,T t) const {return w.Cross(X-c)*(1+sin(t));}
                virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return MATRIX<T,TV::m>::Cross_Product_Matrix(w)*(1+sin(t));}
            };
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),.8,1,0);
            analytic_velocity=new ANALYTIC_VELOCITY_5(TV(),typename TV::SPIN()+1,1);
            analytic_polymer_stress=new ANALYTIC_POLYMER_4;
            break;
        default: return false;}
    return true;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Write_Output_Files(const int frame)
{
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Initialize()
{
    PHYSBAM_ASSERT(analytic_levelset && analytic_velocity && analytic_polymer_stress);
    Analytic_Test();

    if(user_last_frame) last_frame=stored_last_frame;
}
//#####################################################################
// Function Set_Level_Set
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Set_Level_Set(T time)
{
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
        levelset.phi(it.index)=analytic_levelset->phi2(it.Location()/m,time/s)*m;
}
//#####################################################################
// Function Velocity_Error
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Stress_Error(T time)
{
    PHYSBAM_ASSERT(analytic_polymer_stress);
    if(analytic_initial_only) return;
    T S_inf=0,S_2=0,a=0,b=0;
    int num_S=0;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        T p=analytic_levelset->phi2(it.Location()/m,time/s)*m;
        if(p>0) continue;
        SYMMETRIC_MATRIX<T,TV::m> B=analytic_polymer_stress->S(it.Location()/m,time/s)*unit_p;
        SYMMETRIC_MATRIX<T,TV::m> A=polymer_stress(it.index),D=A-B;
        num_S++;
        S_2+=D.Frobenius_Norm_Squared();
        T S_max=D.Max_Abs();
        S_inf=max(S_max,S_inf);
        a=max(A.Max_Abs(),a);
        b=max(B.Max_Abs(),b);
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,S_max);}
    if(num_S) S_2/=num_S*TV::m*TV::m;
    S_2=sqrt(S_2);
    LOG::printf("max_error %-22.16g %-22.16g %-22.16g %-22.16g\n", S_inf, S_2, a, b);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("stress error",0,1);
}
//#####################################################################
// Function Dump_Analytic_Levelset
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Dump_Analytic_Levelset(T time)
{
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        T p=analytic_levelset->phi2(it.Location()/m,time/s)*m;
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p>0,p<=0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("analytic level set (phi)",0,1);
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        T p=analytic_levelset->phi2(it.Location()/m,time/s)*m;
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p>0,p<=0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,analytic_levelset->N2(it.Location()/m,time/s));}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("analytic level set (N)",0,1);
}
//#####################################################################
// Function Get_Velocities
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Get_Velocities(T time)
{
    if(analytic_levelset && analytic_velocity)
        for(FACE_ITERATOR<TV> it(grid,number_of_ghost_cells);it.Valid();it.Next())
            face_velocities(it.Full_Index())=analytic_velocity->u(it.Location()/m,time)(it.Axis())*m/s;
}
//#####################################################################
// Function Get_Initial_Polymer_Stresses
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Get_Initial_Polymer_Stresses()
{
    if(analytic_levelset && analytic_polymer_stress)
        for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
            polymer_stress(it.index)=analytic_polymer_stress->S(it.Location()/m,0)*unit_p;}
}
//#####################################################################
// Function Analytic_Test
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Analytic_Test()
{
    Set_Level_Set(0);
    Dump_Analytic_Levelset(0);

    if(test_analytic_diff){
        analytic_levelset->Test(grid.domain);
        RANDOM_NUMBERS<T> rand;
        TV X;
        do X=rand.Get_Uniform_Vector(grid.domain);
        while(analytic_levelset->phi2(X/m,0)>0);
        analytic_velocity->Test(X);
        analytic_polymer_stress->Test(X);}
}
//#####################################################################
// Function Begin_Time_Step
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
Begin_Time_Step(const T time)
{
}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class TV> void FLUID_STRESS_BASE<TV>::
End_Time_Step(const T time)
{
    Stress_Error(time);
}
//#####################################################################
// Function Polymer_Stress
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> FLUID_STRESS_BASE<TV>::
Polymer_Stress(const TV& X,T time)
{
    return analytic_polymer_stress->S(X/m,time/s)*unit_p;
}
//#####################################################################
// Function Polymer_Stress_Forcing_Term
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> FLUID_STRESS_BASE<TV>::
Polymer_Stress_Forcing_Term(const TV& X,T time)
{
    SYMMETRIC_MATRIX<T,TV::m> m;
    TV u=analytic_velocity->u(X,time);
    for(int i=0;i<TV::m;i++) m+=u(i)*analytic_polymer_stress->dSdX(X,time,i);
    return m+analytic_polymer_stress->dSdt(X,time);
}
//#####################################################################
// Function Volume_Force
//#####################################################################
template<class TV> TV FLUID_STRESS_BASE<TV>::
Volume_Force(const TV& X,T time)
{
    return TV();
}
//#####################################################################
template class FLUID_STRESS_BASE<VECTOR<float,2> >;
template class FLUID_STRESS_BASE<VECTOR<double,2> >;
}
