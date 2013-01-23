//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR__
#define __FLUIDS_COLOR__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ROTATE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SCALE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SPHERE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_TRANSLATE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_VORTEX.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace PhysBAM{

template<class TV>
class FLUIDS_COLOR:public PLS_FC_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef PLS_FC_EXAMPLE<TV> BASE;

public:
    using BASE::grid;using BASE::output_directory;using BASE::domain_boundary;using BASE::face_velocities;
    using BASE::write_substeps_level;using BASE::restart;using BASE::last_frame;using BASE::use_level_set_method;using BASE::use_pls;
    using BASE::dt;using BASE::levelset_color;using BASE::mu;using BASE::rho;using BASE::dump_matrix;using BASE::number_of_colors;
    using BASE::use_advection;using BASE::use_reduced_advection;using BASE::omit_solve;using BASE::use_discontinuous_velocity;
    using BASE::time_steps_per_frame;using BASE::use_p_null_mode;using BASE::Fill_Levelsets_From_Levelset_Color;

    enum WORKAROUND{SLIP=-3,DIRICHLET=-2,NEUMANN=-1}; // From CELL_DOMAIN_INTERFACE_COLOR

    int test_number;
    int resolution;
    int stored_last_frame;
    bool user_last_frame;
    T mu0,mu1;
    T rho0,rho1;
    T m,s,kg;
    int bc_type;
    bool bc_n,bc_d,bc_s;
    bool test_analytic_diff;
    bool no_advection;
    int refine;
    static T Large_Phi() {return 1000;}

    struct ANALYTIC_VELOCITY
    {
        virtual ~ANALYTIC_VELOCITY(){}
        virtual TV u(const TV& X,T t) const=0;
        virtual MATRIX<T,2> du(const TV& X,T t) const=0;
        virtual T p(const TV& X,T t) const=0;
        virtual TV F(const TV& X,T t) const=0;
    };
    ARRAY<ANALYTIC_VELOCITY*> analytic_velocity;
    ANALYTIC_LEVELSET<TV>* analytic_levelset;

    FLUIDS_COLOR(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
        :PLS_FC_EXAMPLE<TV>(stream_type),test_number(0),resolution(32),stored_last_frame(0),user_last_frame(false),mu0(1),mu1(2),
        rho0(1),rho1(2),m(1),s(1),kg(1),bc_n(false),bc_d(false),bc_s(false),test_analytic_diff(false),no_advection(false),refine(1)
    {
        last_frame=16;
        int number_of_threads=1;
        parse_args.Extra(&test_number,"example number","example number to run");
        parse_args.Add("-restart",&restart,"frame","restart frame");
        parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
        parse_args.Add("-substep",&write_substeps_level,"level","output-substep level");
        parse_args.Add("-dt",&dt,"dt","time step size to use for simulation");
        parse_args.Add("-steps",&time_steps_per_frame,"steps","number of time steps per frame");
        parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
        parse_args.Add("-dump_matrix",&dump_matrix,"dump out system and rhs");
        parse_args.Add("-mu0",&mu0,"viscosity","viscosity for first fluid region");
        parse_args.Add("-mu1",&mu1,"viscosity","viscosity for second fluid region");
        parse_args.Add("-rho0",&rho0,"density","density for first fluid region");
        parse_args.Add("-rho1",&rho1,"density","density for second fluid region");
        parse_args.Add("-m",&m,"scale","meter scale");
        parse_args.Add("-s",&s,"scale","second scale");
        parse_args.Add("-kg",&kg,"scale","kilogram scale");
        parse_args.Add("-bc_n",&bc_n,"use Neumann boundary conditions");
        parse_args.Add("-bc_d",&bc_d,"use Dirichlet boundary conditions");
        parse_args.Add("-bc_s",&bc_s,"use slip boundary conditions");
        parse_args.Add("-test_diff",&test_analytic_diff,"test analytic derivatives");
        parse_args.Add("-no_advect",&no_advection,"Disable advection");
        parse_args.Add("-no_solve",&omit_solve,"Disable visocity and pressure solve");
        parse_args.Add("-reduced_advect",&use_reduced_advection,"Peform reduced advection");
        parse_args.Add("-refine",&refine,"num","Refine space/time by this factor");
        parse_args.Add("-null_p",&use_p_null_mode,"Assume pressure null mode and project it out");
        parse_args.Add("-threads",&number_of_threads,"threads","Number of threads");
        parse_args.Parse();

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
        mu0*=kg/s;
        mu1*=kg/s;
        rho0*=kg/sqr(m);
        rho1*=kg/sqr(m);
        dt*=s;
        PHYSBAM_ASSERT(bc_n+bc_d+bc_s<2);
        bc_type=bc_n?NEUMANN:(bc_s?SLIP:DIRICHLET);
        use_advection=!no_advection;

        analytic_levelset=0;

        output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number);
        switch(test_number){
            case 0:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>;
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV()+1));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 1:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>;
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 2:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 3:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION(TV()+(T).5,rho0));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 4:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV()+1));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 5:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 6:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION(TV()+(T).5,rho0));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 7:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>;
                analytic_velocity.Append(new ANALYTIC_VELOCITY_RAREFACTION);
                omit_solve=true;
                break;
            case 8:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>;
                analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg),TV((T).2,(T).5)));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 9:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    T x0=(T).2,x1=(T).5,x2=(T).8,v0=1,v2=-1;
                    T v1=(v0*mu0/(x1-x0)+v2*mu1/(x2-x1))/(mu0/(x1-x0)+mu1/(x2-x1));
                    MATRIX<T,2> du0(0,(v1-v0)/(x1-x0),0,0),du1(0,(v2-v1)/(x2-x1),0,0);
                    ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV(x0,0),TV(1,0),DIRICHLET,0);
                    ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV(x2,0),TV(1,0),1,DIRICHLET);
                    analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV(x1,0),TV(1,0),0,1)))->Add(ab)->Add(cd);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE(TV(x1,0),TV(0,v1),du0,rho0*sqr(m)/kg));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE(TV(x1,0),TV(0,v1),du1,rho1*sqr(m)/kg));
                    use_p_null_mode=true;
                }
                break;
            case 10:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV()+1));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 11:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION(TV((T).6,(T).8),rho0*sqr(m)/kg));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 12:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    T x0=(T).2,x1=(T).8,g=9.8,a=rho0*sqr(m)*g/(2*mu0*s);
                    ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV(x0,0),TV(1,0),DIRICHLET,0);
                    ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),DIRICHLET);
                    analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV(x1,0),TV(1,0),0,1)))->Add(ab)->Add(cd);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_QUADRATIC_X(a,-a*(x0+x1),a*x0*x1,mu0*s/kg));
                    use_p_null_mode=true;
                }
                break;
            case 13:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    T x0=(T).2,x1=(T).5,x2=(T).8,v0=1,v1=1,v2=-1;
                    MATRIX<T,2> du0(0,(v1-v0)/(x1-x0),0,0),du1(0,(v2-v1)/(x2-x1),0,0);
                    ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV(x0,0),TV(1,0),DIRICHLET,0);
                    ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV(x2,0),TV(1,0),1,DIRICHLET);
                    analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV(x1,0),TV(1,0),0,1)))->Add(ab)->Add(cd);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE(TV(x1,0),TV(0,v1),du0,rho0*sqr(m)/kg));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE(TV(x1,0),TV(0,v1),du1,rho1*sqr(m)/kg));
                    use_p_null_mode=true;
                }
                break;
            case 14:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    T x0=(T).2,x1=(T).5,x2=(T).8,v0=1,v2=-1;
                    ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV(x0,0),TV(1,0),DIRICHLET,0);
                    ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV(x2,0),TV(1,0),1,DIRICHLET);
                    analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV(x1,0),TV(1,0),0,1)))->Add(ab)->Add(cd);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV(0,v0)));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV(0,v2)));
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
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_QUADRATIC_X((T).5,(T).2,(T).3,mu0*s/kg));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION(TV((T).1,(T).1),rho1*sqr(m)/kg));
                    use_discontinuous_velocity=true;
                    use_p_null_mode=true;
                }
                break;
            case 16:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                {
                    TV vel((T).2,(T).5);
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_VORTEX<TV>((T).2),vel);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg),vel));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                }
                break;
            case 17:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    TV vel((T).2,(T).5);
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3),vel);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE(new ANALYTIC_VELOCITY_CONST(TV()+1),vel));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                }
                break;
            case 18:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3),TV(.1,.2));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE(new ANALYTIC_VELOCITY_ROTATION(TV()+1,rho0*sqr(m)/kg),TV(0.5,-0.4)));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                }
                break;
            case 19:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                {
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_SCALE<TV>(new ANALYTIC_LEVELSET_VORTEX<TV>((T).2),-2),4),TV(9,.4));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg),TV(.5,-.2)));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                }
                break;
            case 20:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    TV vel((T).2,(T).5);
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,1,0),vel);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(vel));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(vel));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                    use_level_set_method=true;
                }
                break;
            case 21:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    TV vel(1,0);
                    T x0=(T).2,x1=(T).5,x2=(T).8;
                    ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV(x0,0),TV(1,0),DIRICHLET,0);
                    ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV(x2,0),TV(1,0),1,DIRICHLET);
                    analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV(x1,0),TV(1,0),0,1),vel)))->Add(ab)->Add(cd);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(vel));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(vel));
                    use_p_null_mode=true;
                    use_level_set_method=true;
                }
                break;
            default: PHYSBAM_FATAL_ERROR("Missing test number");}

        if(analytic_velocity.m) number_of_colors=analytic_velocity.m;
    }

    ~FLUIDS_COLOR()
    {
        delete analytic_levelset;
        analytic_velocity.Delete_Pointers_And_Clean_Memory();
    }

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);}

    void Initialize()
    {
        if(analytic_levelset && analytic_velocity.m) Analytic_Test();
        else
            switch(test_number){
                default: PHYSBAM_FATAL_ERROR("Missing test number");}

        PHYSBAM_ASSERT(rho.m==number_of_colors && mu.m==number_of_colors);
        if(user_last_frame) last_frame=stored_last_frame;
    }

    void Set_Level_Set(T time)
    {
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=-4;
            T p=analytic_levelset->phi(it.Location()/m,time,c)*m;
            levelset_color.phi(it.index)=abs(p);
            levelset_color.color(it.index)=c==-4?bc_type:c;}
        Fill_Levelsets_From_Levelset_Color();
    }

    void Level_Set_Error(T time)
    {
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=levelset_color.color(it.index);
            T p=levelset_color.phi(it.index);
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(c<0,c>=0,p<0));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("level set",0,1);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=-4;
            T p=analytic_levelset->phi(it.Location()/m,time,c)*m;
            if(levelset_color.color(it.index)!=(c==-4?bc_type:c)){
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(levelset_color.phi(it.index))+abs(p));}
            else{
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(levelset_color.phi(it.index)-p));}}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("level set error",0,1);
    }

    void Dump_Analytic_Levelset(T time)
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
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=-4;
            T p=analytic_levelset->phi(it.Location()/m,time,c)*m;
            if(c==-4) c=bc_type;
            Add_Debug_Particle(it.Location(),colors(c+3));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("analytic level set (phi)",0,1);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=-4;
            T p=analytic_levelset->phi(it.Location()/m,time,c)*m;
            (void)p;
            if(c==-4) c=bc_type;
            Add_Debug_Particle(it.Location(),colors(c+3));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,analytic_levelset->N(it.Location()/m,time,c));}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("analytic level set (N)",0,1);
    }

    void Analytic_Test()
    {
        mu.Append(mu0);
        rho.Append(rho0);
        for(int i=1;i<number_of_colors;i++){
            mu.Append(mu1);
            rho.Append(rho1);}
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
            int c=-4;
            analytic_levelset->phi(it.Location()/m,0,c);
            if(c<0) continue;
            face_velocities(c)(it.Full_Index())=analytic_velocity(c)->u(it.Location()/m,0)(it.Axis())*m/s;}

        Set_Level_Set(0);
        Dump_Analytic_Levelset(0);

        if(test_analytic_diff){
            RANDOM_NUMBERS<T> rand;
            TV X,dX;
            T e=1e-6,t=rand.Get_Uniform_Number(0,1);
            rand.Fill_Uniform(dX,-e,e);
            int c=-4;
            do{
                X=rand.Get_Uniform_Vector(grid.domain);
                analytic_levelset->phi(X/m,0,c);}
            while(c<0);
            TV u0=analytic_velocity(c)->u(X/m,t)*m/s,u1=analytic_velocity(c)->u((X+dX)/m,t)*m/s;
            MATRIX<T,2> du0=analytic_velocity(c)->du(X/m,t)/s,du1=analytic_velocity(c)->du((X+dX)/m,t)/s;
            T erru=((du0+du1)*dX/2-(u1-u0)).Magnitude()/e;
            int c0,c1;
            T l0=analytic_levelset->phi(X/m,t,c0)*m,l1=analytic_levelset->phi((X+dX)/m,t,c1)*m;
            if(c0>=0) l0=-l0;
            if(c1>=0) l1=-l1;
            TV dl0=analytic_levelset->N(X/m,t,c),dl1=analytic_levelset->N((X+dX)/m,t,c);
            T errl=abs((dl0+dl1).Dot(dX)/2-(l1-l0))/e;
            LOG::cout<<"analytic diff test "<<erru<<"  "<<errl<<std::endl;}
    }

    struct ANALYTIC_VELOCITY_CONST:public ANALYTIC_VELOCITY
    {
        TV au;
        ANALYTIC_VELOCITY_CONST(TV v): au(v){}
        virtual TV u(const TV& X,T t) const {return au;}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>();}
        virtual T p(const TV& X,T t) const {return 0;}
        virtual TV F(const TV& X,T t) const {return TV();}
    };

    struct ANALYTIC_VELOCITY_ROTATION:public ANALYTIC_VELOCITY
    {
        TV c;
        T rho;
        ANALYTIC_VELOCITY_ROTATION(TV cc,T rho): c(cc),rho(rho){}
        virtual TV u(const TV& X,T t) const {return (X-c).Orthogonal_Vector();}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>(0,1,-1,0);}
        virtual T p(const TV& X,T t) const {return (T).5*rho*(X-c).Magnitude_Squared();}
        virtual TV F(const TV& X,T t) const {return TV();}
    };

    struct ANALYTIC_VELOCITY_VORTEX:public ANALYTIC_VELOCITY
    {
        T nu,rho;
        ANALYTIC_VELOCITY_VORTEX(T mu,T rho): nu(mu/rho),rho(rho){}
        virtual TV u(const TV& X,T t) const
        {return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y))*exp(-2*nu*t);}
        virtual MATRIX<T,2> du(const TV& X,T t) const
        {T c=cos(X.x)*cos(X.y),s=sin(X.x)*sin(X.y);return MATRIX<T,2>(c,s,-s,-c)*exp(-2*nu*t);}
        virtual T p(const TV& X,T t) const {return (T).25*rho*(cos(2*X.x)+cos(2*X.y))*exp(-4*nu*t);}
        virtual TV F(const TV& X,T t) const {return TV();}
    };

    struct ANALYTIC_VELOCITY_RAREFACTION:public ANALYTIC_VELOCITY
    {
        ANALYTIC_VELOCITY_RAREFACTION(){}
        virtual TV u(const TV& X,T t) const {return TV(X.x,0)/(t+1);}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>(1,0,0,0)/(t+1);}
        virtual T p(const TV& X,T t) const {return 0;}
        virtual TV F(const TV& X,T t) const {return TV();}
    };

    struct ANALYTIC_VELOCITY_ROTATION_DECAY:public ANALYTIC_VELOCITY
    {
        ANALYTIC_VELOCITY_ROTATION_DECAY(){}
        virtual TV u(const TV& X,T t) const {return X.Orthogonal_Vector()*2*exp(-X.Magnitude_Squared());}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>(-2*X.x*X.y,1-2*sqr(X.y),-1+2*sqr(X.x),2*X.x*X.y)*2*exp(-X.Magnitude_Squared());}
        virtual T p(const TV& X,T t) const {return 0;}
        virtual TV F(const TV& X,T t) const {return TV();}
    };

    struct ANALYTIC_VELOCITY_AFFINE:public ANALYTIC_VELOCITY
    {
        TV v0;
        MATRIX<T,2> du0;
        T rho;
        ANALYTIC_VELOCITY_AFFINE(const TV& x0,const TV& v0,const MATRIX<T,2>& du0,T rho): v0(v0-du0*x0),du0(du0) {}
        virtual TV u(const TV& X,T t) const {return du0*X+v0;}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return du0;}
        virtual T p(const TV& X,T t) const {return 0;}
        virtual TV F(const TV& X,T t) const {return rho*du0*(du0*X+v0);}
    };

    struct ANALYTIC_VELOCITY_QUADRATIC_X:public ANALYTIC_VELOCITY
    {
        T a,b,c,mu;
        ANALYTIC_VELOCITY_QUADRATIC_X(T a,T b,T c,T mu): a(a),b(b),c(c),mu(mu) {}
        virtual TV u(const TV& X,T t) const {return TV(0,(a*X.x+b)*X.x+c);}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>(0,2*a*X.x+b,0,0);}
        virtual T p(const TV& X,T t) const {return 0;}
        virtual TV F(const TV& X,T t) const {return TV(0,-2*a*mu);}
    };

    struct ANALYTIC_VELOCITY_TRANSLATE:public ANALYTIC_VELOCITY
    {
        ANALYTIC_VELOCITY* av;
        TV vel;

        ANALYTIC_VELOCITY_TRANSLATE(ANALYTIC_VELOCITY* av,const TV& vel): av(av),vel(vel) {}
        ~ANALYTIC_VELOCITY_TRANSLATE() {delete av;}
        virtual TV u(const TV& X,T t) const {return av->u(X-vel*t,t)+vel;}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return av->du(X-vel*t,t);}
        virtual T p(const TV& X,T t) const {return av->p(X-vel*t,t);}
        virtual TV F(const TV& X,T t) const {return av->F(X-vel*t,t);}
    };

    void Begin_Time_Step(const T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset && !use_level_set_method && !use_pls)
            Set_Level_Set(time+dt);
    }

    void End_Time_Step(const T time) PHYSBAM_OVERRIDE
    {
        Level_Set_Error(time);
        if(analytic_velocity.m){
            T max_error=0,a=0,b=0;
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                int c=levelset_color.Color(it.Location());
                if(c<0) continue;
                T A=face_velocities(c)(it.Full_Index()),B=analytic_velocity(c)->u(it.Location()/m,time/s)(it.Axis())*m/s;
                a=max(a,abs(A));
                b=max(b,abs(B));
                max_error=max(max_error,abs(A-B));
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(A-B));}
            LOG::cout<<"max_error "<<max_error<<"  "<<a<<"  "<<b<<std::endl;}
    }

    TV Dirichlet_Boundary_Condition(const TV& X,int bc_color,int fluid_color,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
        if(analytic_velocity.m) return analytic_velocity(fluid_color)->u(X/m,time/s)*m/s;
        return TV();
    }

    TV Neumann_Boundary_Condition(const TV& X,int bc_color,int fluid_color,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));
        if(analytic_velocity.m && analytic_levelset){
            MATRIX<T,2> du=analytic_velocity(fluid_color)->du(X/m,time/s)/s;
            TV n=analytic_levelset->N(X/m,time/s,fluid_color);
            T p=analytic_velocity(fluid_color)->p(X/m,time/s)*kg/(s*s*m);
            return (du+du.Transposed())*n*mu(fluid_color)-p*n;}
        return TV();
    }

    TV Jump_Interface_Condition(const TV& X,int color0,int color1,T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset){
            MATRIX<T,2> du0=analytic_velocity(color0)->du(X/m,time/s)/s,du1=analytic_velocity(color1)->du(X/m,time/s)/s;
            T p0=analytic_velocity(color0)->p(X/m,time/s)*kg/(s*s*m),p1=analytic_velocity(color1)->p(X/m,time/s)*kg/(s*s*m);
            MATRIX<T,2> stress0=(du0+du0.Transposed())*mu(color0)-p0,stress1=(du1+du1.Transposed())*mu(color1)-p1;
            TV n=analytic_levelset->N(X/m,time/s,color1);
            return (stress0-stress1)*n;}
        return TV();
    }

    TV Volume_Force(const TV& X,int color,T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset)
            return analytic_velocity(color)->F(X/m,time/s)*kg/(m*s*s);
        return TV();
    }

    TV Velocity_Jump(const TV& X,int color0,int color1,T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset)
            return (analytic_velocity(color1)->u(X/m,time/s)-analytic_velocity(color0)->u(X/m,time/s))*m/s;
        return TV();
    }

//#####################################################################
};
}

#endif
