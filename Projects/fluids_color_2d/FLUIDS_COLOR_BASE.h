//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR_BASE__
#define __FLUIDS_COLOR_BASE__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_BOX.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ELLIPSOID.h>
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
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace PhysBAM{

template<class TV> class FLUIDS_COLOR;

template<class TV>
struct ANALYTIC_VELOCITY
{
    typedef typename TV::SCALAR T;
    virtual ~ANALYTIC_VELOCITY(){}
    virtual TV u(const TV& X,T t) const=0;
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const=0;
    virtual T p(const TV& X,T t) const=0;
    virtual TV F(const TV& X,T t) const=0;
};

template<class TV>
struct ANALYTIC_VELOCITY_CONST:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV au;
    T const_p;
    ANALYTIC_VELOCITY_CONST(TV v): au(v),const_p(0) {}
    virtual TV u(const TV& X,T t) const {return au;}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return MATRIX<T,TV::m>();}
    virtual T p(const TV& X,T t) const {return const_p;}
    virtual TV F(const TV& X,T t) const {return TV();}
};

template<class TV>
struct ANALYTIC_VELOCITY_ROTATION:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV c;
    typename TV::SPIN w;
    T rho;
    ANALYTIC_VELOCITY_ROTATION(TV cc,typename TV::SPIN ww,T rho): c(cc),w(ww),rho(rho){}
    virtual TV u(const TV& X,T t) const {return w.Cross(X-c);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return MATRIX<T,TV::m>::Cross_Product_Matrix(w);}
    virtual T p(const TV& X,T t) const {return (T).5*rho*u(X,t).Magnitude_Squared();}
    virtual TV F(const TV& X,T t) const {return TV();}
};

template<class TV>
struct ANALYTIC_VELOCITY_RAREFACTION:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_VELOCITY_RAREFACTION(){}
    virtual TV u(const TV& X,T t) const {return X.x/(t+1)*TV::Axis_Vector(0);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {MATRIX<T,TV::m> M;M(0,0)=1/(t+1);return M;}
    virtual T p(const TV& X,T t) const {return 0;}
    virtual TV F(const TV& X,T t) const {return TV();}
};

template<class TV>
struct ANALYTIC_VELOCITY_AFFINE:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    TV v0;
    MATRIX<T,TV::m> du0;
    T rho;
    ANALYTIC_VELOCITY_AFFINE(const TV& x0,const TV& v0,const MATRIX<T,TV::m>& du0,T rho): v0(v0-du0*x0),du0(du0) {}
    virtual TV u(const TV& X,T t) const {return du0*X+v0;}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return du0;}
    virtual T p(const TV& X,T t) const {return 0;}
    virtual TV F(const TV& X,T t) const {return rho*du0*(du0*X+v0);}
};

template<class TV>
struct ANALYTIC_VELOCITY_QUADRATIC_X:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    T a,b,c,mu;
    ANALYTIC_VELOCITY_QUADRATIC_X(T a,T b,T c,T mu): a(a),b(b),c(c),mu(mu) {}
    virtual TV u(const TV& X,T t) const {return TV::Axis_Vector(1)*((a*X.x+b)*X.x+c);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {MATRIX<T,TV::m> M;M(1,0)=2*a*X.x+b;return M;}
    virtual T p(const TV& X,T t) const {return 0;}
    virtual TV F(const TV& X,T t) const {return -2*a*mu*TV::Axis_Vector(1);}
};

template<class TV>
struct ANALYTIC_VELOCITY_TRANSLATE:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_VELOCITY<TV>* av;
    TV vel;
    ANALYTIC_VELOCITY_TRANSLATE(ANALYTIC_VELOCITY<TV>* av,const TV& vel): av(av),vel(vel) {}
    ~ANALYTIC_VELOCITY_TRANSLATE() {delete av;}
    virtual TV u(const TV& X,T t) const {return av->u(X-vel*t,t)+vel;}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const {return av->du(X-vel*t,t);}
    virtual T p(const TV& X,T t) const {return av->p(X-vel*t,t);}
    virtual TV F(const TV& X,T t) const {return av->F(X-vel*t,t);}
};

template<class TV>
class FLUIDS_COLOR_BASE:public PLS_FC_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef PLS_FC_EXAMPLE<TV> BASE;

public:
    using BASE::grid;using BASE::output_directory;using BASE::face_velocities;
    using BASE::write_substeps_level;using BASE::restart;using BASE::last_frame;using BASE::use_level_set_method;using BASE::use_pls;
    using BASE::dt;using BASE::levelset_color;using BASE::mu;using BASE::rho;using BASE::dump_matrix;using BASE::number_of_colors;
    using BASE::use_advection;using BASE::use_reduced_advection;using BASE::omit_solve;using BASE::use_discontinuous_velocity;
    using BASE::time_steps_per_frame;using BASE::use_p_null_mode;using BASE::Fill_Levelsets_From_Levelset_Color;
    using BASE::particle_levelset_evolution_multiple;using BASE::face_color;using BASE::substeps_delay_frame;

    enum WORKAROUND{SLIP=-3,DIRICHLET=-2,NEUMANN=-1}; // From CELL_DOMAIN_INTERFACE_COLOR

    int test_number;
    int resolution;
    int stored_last_frame;
    bool user_last_frame;
    T mu0,mu1;
    T rho0,rho1;
    T unit_mu,unit_rho,unit_st,unit_p;
    T m,s,kg;
    int bc_type;
    bool bc_n,bc_d,bc_s;
    bool test_analytic_diff;
    bool no_advection;
    int refine;
    static T Large_Phi() {return 1000;}
    T surface_tension;
    bool override_rho0;
    bool override_rho1;
    bool override_mu0;
    bool override_mu1;
    bool use_pls_over_levelset;

    ARRAY<ANALYTIC_VELOCITY<TV>*> analytic_velocity,initial_analytic_velocity;
    ANALYTIC_LEVELSET<TV>* analytic_levelset;
    bool analytic_initial_only;
    int number_of_threads;
    bool override_output_directory;

    FLUIDS_COLOR_BASE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
        :PLS_FC_EXAMPLE<TV>(stream_type),test_number(0),resolution(32),stored_last_frame(0),user_last_frame(false),mu0(1),mu1(2),rho0(1),
        rho1(2),unit_mu(0),unit_rho(0),unit_st(0),unit_p(0),m(1),s(1),kg(1),bc_n(false),bc_d(false),bc_s(false),test_analytic_diff(false),
        no_advection(false),refine(1),surface_tension(0),override_rho0(false),override_rho1(false),override_mu0(false),override_mu1(false),
        use_pls_over_levelset(false),analytic_initial_only(false),number_of_threads(1),override_output_directory(false)
    {
        last_frame=16;
        parse_args.Extra(&test_number,"example number","example number to run");
        parse_args.Add("-restart",&restart,"frame","restart frame");
        parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
        parse_args.Add("-substep",&write_substeps_level,"level","output-substep level");
        parse_args.Add("-substeps_delay",&substeps_delay_frame,"frame","delay substeps until after this frame");
        parse_args.Add("-dt",&dt,"dt","time step size to use for simulation");
        parse_args.Add("-steps",&time_steps_per_frame,"steps","number of time steps per frame");
        parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
        parse_args.Add("-dump_matrix",&dump_matrix,"dump out system and rhs");
        parse_args.Add("-use_pls",&use_pls_over_levelset,"use_particle_level_set");
        parse_args.Add("-mu0",&mu0,&override_mu0,"viscosity","viscosity for first fluid region");
        parse_args.Add("-mu1",&mu1,&override_mu1,"viscosity","viscosity for second fluid region");
        parse_args.Add("-rho0",&rho0,&override_rho0,"density","density for first fluid region");
        parse_args.Add("-rho1",&rho1,&override_rho1,"density","density for second fluid region");
        parse_args.Add("-m",&m,"scale","meter scale");
        parse_args.Add("-s",&s,"scale","second scale");
        parse_args.Add("-kg",&kg,"scale","kilogram scale");
        parse_args.Add("-bc_n",&bc_n,"use Neumann boundary conditions");
        parse_args.Add("-bc_d",&bc_d,"use Dirichlet boundary conditions");
        parse_args.Add("-bc_s",&bc_s,"use slip boundary conditions");
        parse_args.Add("-test_diff",&test_analytic_diff,"test analytic derivatives");
        parse_args.Add("-no_advect",&no_advection,"Disable advection");
        parse_args.Add("-no_solve",&omit_solve,"Disable visocity and pressure solve");
        parse_args.Add_Not("-no_reduced_advect",&use_reduced_advection,"Peform reduced advection");
        parse_args.Add("-refine",&refine,"num","Refine space/time by this factor");
        parse_args.Add("-null_p",&use_p_null_mode,"Assume pressure null mode and project it out");
        parse_args.Add("-threads",&number_of_threads,"threads","Number of threads");
        parse_args.Add("-o",&output_directory,&override_output_directory,"dir","Output directory");
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
        unit_mu=kg*pow<TV::m-2>(m)/s;
        unit_rho=kg/pow<TV::m>(m);
        unit_st=kg*pow<3-TV::m>(m)/(s*s);
        unit_p=kg*pow<2-TV::m>(m)/(s*s);
        mu0*=unit_mu;
        mu1*=unit_mu;
        rho0*=unit_rho;
        rho1*=unit_rho;
        surface_tension*=unit_st;
        dt*=s;
        PHYSBAM_ASSERT(bc_n+bc_d+bc_s<2);
        bc_type=bc_n?NEUMANN:(bc_s?SLIP:DIRICHLET);
        use_advection=!no_advection;

        analytic_levelset=0;
    }

    ~FLUIDS_COLOR_BASE()
    {
        delete analytic_levelset;
        analytic_velocity.Delete_Pointers_And_Clean_Memory();
    }

    void After_Initialize_Example()
    {
        if(!override_output_directory) output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number);

        if(use_pls_over_levelset){use_pls=use_level_set_method;use_level_set_method=false;}

        if(analytic_velocity.m) number_of_colors=analytic_velocity.m;
    }

    bool Initialize_Common_Example()
    {
        typename TV::SPIN spin_count;
        TV vector_count;
        for(int i=0;i<TV::SPIN::m;i++) spin_count(i)=i;
        for(int i=0;i<TV::m;i++) vector_count(i)=i;

        switch(test_number){
            case 0:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV()+1));
                use_p_null_mode=true;
                break;
            case 3:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).5,spin_count+1,rho0/unit_rho));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 4:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV()+1));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 7:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_RAREFACTION<TV>);
                omit_solve=true;
                break;
            case 9:
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
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(1)*v1,du0,rho0/unit_rho));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(1)*v1,du1,rho1/unit_rho));
                    use_p_null_mode=true;
                }
                break;
            case 11:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(vector_count*(T).2+(T).6,spin_count+1,rho0/unit_rho));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 12:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    T x0=(T).2,x1=(T).8,g=9.8,a=rho0/unit_rho*g/(2*mu0/unit_mu);
                    ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x0,TV::Axis_Vector(0),DIRICHLET,0);
                    ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),DIRICHLET,DIRICHLET);
                    analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_QUADRATIC_X<TV>(a,-a*(x0+x1),a*x0*x1,mu0/unit_mu));
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
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(1)*v1,du0,rho0/unit_rho));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE<TV>(TV::Axis_Vector(0)*x1,TV::Axis_Vector(1)*v1,du1,rho1/unit_rho));
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
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV::Axis_Vector(1)*v0));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV::Axis_Vector(1)*v2));
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
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_QUADRATIC_X<TV>((T).5,(T).2,(T).3,mu0/unit_mu));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).1,spin_count+1,rho1/unit_rho));
                    use_discontinuous_velocity=true;
                    use_p_null_mode=true;
                }
                break;
            case 17:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    TV vel=vector_count*(T).3+(T).2;
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4),vel);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE<TV>(new ANALYTIC_VELOCITY_CONST<TV>(TV()+1),vel));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                }
                break;
            case 18:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4),vector_count*(T).1+(T).1);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE<TV>(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+1,spin_count+1,rho0/unit_rho),vector_count*(T)-.9+(T).5));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                }
                break;
            case 20:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    TV vel=vector_count*(T)-.7+(T).2;
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,1,0),vel);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(vel));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(vel));
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
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(vel));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(vel));
                    use_p_null_mode=true;
                    use_level_set_method=true;
                }
                break;
            case 22:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,1);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV()));
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV()));
                surface_tension=unit_st;
                use_p_null_mode=true;
                use_level_set_method=true;
                break;
            case 24:{
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
                ANALYTIC_VELOCITY_CONST<TV>* a=new ANALYTIC_VELOCITY_CONST<TV>(TV());
                a->const_p=1;
                analytic_velocity.Append(a);
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;}
            case 25:{//A circle being rotated around inside a larger circle
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);

                ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(vector_count*(T)-.2+(T).5,.1,1,0),spin_count+1,TV()+(T).5);
                ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),-4,-4);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).4,0,1)))->Add(ab)->Add(cd);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).5,spin_count+1,rho0/unit_rho));
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).5,spin_count+1,rho0/unit_rho));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                //use_level_set_method=true;

                break;}
            case 26:{
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                TV vel=vector_count*(T).3+(T).2;
                //analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_BOX<TV>(TV()+(T).2,TV()+(T).45,1,0),vel);
                analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_BOX<TV>(TV()+(T).2,TV()+(T).45,1,0),spin_count*0,TV()+(T).5),vel);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(vel));
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(vel));
                use_p_null_mode=true;
                use_level_set_method=true;
                break;}
            case 28:{//Test 25 but using the level set method
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                
                ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(vector_count*(T)-.1+(T).5,.1,1,0),spin_count+1,TV()+(T).5);
                ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),-4,-4);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).4,0,1)))->Add(ab)->Add(cd);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).5,spin_count+1,rho0/unit_rho));
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).5,spin_count+1,rho0/unit_rho));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                use_level_set_method=true;
                
                break;}
            default: return false;}
        return true;
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
            T p=analytic_levelset->phi(it.Location()/m,time/s,c)*m;
            levelset_color.phi(it.index)=abs(p);
            levelset_color.color(it.index)=c==-4?bc_type:c;}
        Fill_Levelsets_From_Levelset_Color();
    }

    void Level_Set_Error(T time)
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

        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=levelset_color.color(it.index);
            T p=levelset_color.phi(it.index);
            Add_Debug_Particle(it.Location(),colors(c+3));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("level set",0,1);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
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

    void Velocity_Error(T time)
    {
        if(!analytic_velocity.m || analytic_initial_only) return;
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
        LOG::cout<<"max_error "<<max_error<<"  "<<a<<"  "<<b<<std::endl;
        PHYSBAM_DEBUG_WRITE_SUBSTEP("velocity error",0,1);
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
            T p=analytic_levelset->phi(it.Location()/m,time/s,c)*m;
            if(c==-4) c=bc_type;
            Add_Debug_Particle(it.Location(),colors(c+3));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("analytic level set (phi)",0,1);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=-4;
            T p=analytic_levelset->phi(it.Location()/m,time/s,c)*m;
            (void)p;
            if(c==-4) c=bc_type;
            Add_Debug_Particle(it.Location(),colors(c+3));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,analytic_levelset->N(it.Location()/m,time/s,c));}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("analytic level set (N)",0,1);
    }

    void Get_Initial_Velocities()
    {
        if(analytic_levelset && analytic_velocity.m)
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                int c=face_color(it.Full_Index());
                if(c<0) continue;
                face_velocities(c)(it.Full_Index())=analytic_velocity(c)->u(it.Location()/m,0)(it.Axis())*m/s;}
    }

    void Analytic_Test()
    {
        mu.Append(mu0);
        rho.Append(rho0);
        for(int i=1;i<number_of_colors;i++){
            mu.Append(mu1);
            rho.Append(rho1);}

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
            MATRIX<T,TV::m> du0=analytic_velocity(c)->du(X/m,t)/s,du1=analytic_velocity(c)->du((X+dX)/m,t)/s;
            T erru=((du0+du1)*dX/2-(u1-u0)).Magnitude()/e;
            int c0,c1;
            T l0=analytic_levelset->phi(X/m,t/s,c0)*m,l1=analytic_levelset->phi((X+dX)/m,t/s,c1)*m;
            if(c0>=0) l0=-l0;
            if(c1>=0) l1=-l1;
            TV dl0=analytic_levelset->N(X/m,t/s,c),dl1=analytic_levelset->N((X+dX)/m,t/s,c);
            T errl=abs((dl0+dl1).Dot(dX)/2-(l1-l0))/e;
            LOG::cout<<"analytic diff test "<<erru<<"  "<<errl<<std::endl;}
    }

    void Begin_Time_Step(const T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset && !use_level_set_method && !use_pls && !analytic_initial_only)
            Set_Level_Set(time+dt);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
            for (int c=0;c<face_velocities.m;c++){
            face_velocities(c)(it.Full_Index())=analytic_velocity(c)->u(it.Location()/m,time/s)(it.Axis())*m/s;
              }}

    }

    void End_Time_Step(const T time) PHYSBAM_OVERRIDE
    {
        Level_Set_Error(time);
        Velocity_Error(time);
    }

    MATRIX<T,TV::m> Stress(const TV& X,int color,T time)
    {
        T p=analytic_velocity(color)->p(X/m,time/s)*unit_p;
        MATRIX<T,TV::m> du=analytic_velocity(color)->du(X/m,time/s)/s;
        return (du+du.Transposed())*mu(color)-p;
    }

    TV Jump_Interface_Condition(const TV& X,int color0,int color1,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));

        if(surface_tension){
            T k=particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.levelsets(color1)->Compute_Curvature(X);
            TV n=particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.levelsets(color1)->Normal(X);
            Add_Debug_Particle(X,VECTOR<T,3>(0,0,1));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,k*surface_tension*n);
            return k*surface_tension*n;}

        if(analytic_velocity.m && analytic_levelset && !analytic_initial_only){
            MATRIX<T,TV::m> jump_stress=Stress(X,color1,time);
            if(color0>=0) jump_stress-=Stress(X,color0,time);
            TV n=analytic_levelset->N(X/m,time/s,color1);
            return jump_stress*n;}

        return TV();
    }

    TV Volume_Force(const TV& X,int color,T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset && !analytic_initial_only)
            return analytic_velocity(color)->F(X/m,time/s)*kg/(m*s*s);
        return TV();
    }

    TV Velocity_Jump(const TV& X,int color0,int color1,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
        if(analytic_velocity.m && analytic_levelset && !analytic_initial_only){
            TV jump_u=analytic_velocity(color1)->u(X/m,time/s);
            if(color0>=0) jump_u-=analytic_velocity(color0)->u(X/m,time/s);
            return jump_u*m/s;}
        return TV();
    }

//#####################################################################
};
}

#endif
