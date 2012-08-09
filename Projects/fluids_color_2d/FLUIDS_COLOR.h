//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR__
#define __FLUIDS_COLOR__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>

namespace PhysBAM{

template<class TV>
class FLUIDS_COLOR:public PLS_FC_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef PLS_FC_EXAMPLE<TV> BASE;

public:
    using BASE::grid;using BASE::output_directory;using BASE::domain_boundary;using BASE::face_velocities;
    using BASE::particle_levelset_evolution;using BASE::write_substeps_level;using BASE::restart;using BASE::last_frame;
    using BASE::dt;using BASE::levelset_color;using BASE::mu;using BASE::rho;using BASE::dump_matrix;
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

    struct ANALYTIC_VELOCITY
    {
        virtual TV u(const TV& X,T t) const=0;
        virtual MATRIX<T,2> du(const TV& X,T t) const=0;
        virtual T p(const TV& X,T t) const=0;
    }* analytic_velocity;

    struct ANALYTIC_LEVELSET
    {
        virtual T phi(const TV& X,T t,int& c) const=0;
        virtual TV N(const TV& X,T t) const=0;
    }* analytic_levelset;

    FLUIDS_COLOR(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
        :PLS_FC_EXAMPLE<TV>(stream_type),test_number(0),resolution(32),stored_last_frame(0),user_last_frame(false),mu0(1),mu1(2),
        rho0(1),rho1(2),m(1),s(1),kg(1),bc_n(false),bc_d(false),bc_s(false)
    {
        last_frame=16;
        parse_args.Set_Extra_Arguments(1,"<test-number>");
        parse_args.Add("-restart",&restart,"frame","restart frame");
        parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
        parse_args.Add("-substep",&write_substeps_level,"level","output-substep level");
        parse_args.Add("-dt",&dt,"dt","time step size to use for simulation");
        parse_args.Add("-steps",&this->time_steps_per_frame,"steps","number of time steps per frame");
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
        parse_args.Parse();
        if(!STRING_UTILITIES::String_To_Value(parse_args.Extra_Arg(0),test_number)) throw VALUE_ERROR("The argument is not an integer.");

        stored_last_frame=last_frame;
        mu0*=kg/s;
        mu1*=kg/s;
        rho0*=kg/sqr(m);
        rho1*=kg/sqr(m);
        dt*=s;
        PHYSBAM_ASSERT(bc_n+bc_d+bc_s<2);
        bc_type=bc_n?NEUMANN:(bc_s?SLIP:DIRICHLET);

        analytic_velocity=0;
        analytic_levelset=0;

        output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number);
        switch(test_number){
            case 0: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);break;
            case 1: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);break;
            case 2: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);break;
            case 3: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);break;
            case 4: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);break;
            case 5: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);break;
            default: PHYSBAM_FATAL_ERROR("Missing test number");}
    }

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);}

    void Initialize()
    {
        switch(test_number){
            case 0: Analytic_Test(new ANALYTIC_LEVELSET_PERIODIC,new ANALYTIC_VELOCITY_CONST(TV()+1),1);break;
            case 1: Analytic_Test(new ANALYTIC_LEVELSET_PERIODIC,new ANALYTIC_VELOCITY_VORTEX(mu0,rho0),1);break;
            case 2: Analytic_Test(new ANALYTIC_LEVELSET_VORTEX((T).2),new ANALYTIC_VELOCITY_VORTEX(mu0,rho0),1);break;
            case 3: Analytic_Test(new ANALYTIC_LEVELSET_CIRCLE(TV()+(T).5,(T).3),new ANALYTIC_VELOCITY_ROTATION(TV()+(T).5),1);break;
            case 4: Analytic_Test(new ANALYTIC_LEVELSET_CIRCLE(TV()+(T).5,(T).3),new ANALYTIC_VELOCITY_CONST(TV()+1),1);break;
            case 5: Analytic_Test(new ANALYTIC_LEVELSET_CIRCLE(TV()+(T).5,(T).3),new ANALYTIC_VELOCITY_VORTEX(mu0,rho0),1);break;
            default: PHYSBAM_FATAL_ERROR("Missing test number");}

        if(user_last_frame) last_frame=stored_last_frame;
    }

    void Analytic_Test(ANALYTIC_LEVELSET* ls,ANALYTIC_VELOCITY* v,int colors)
    {
        analytic_levelset=ls;
        analytic_velocity=v;
        mu.Append(mu0);
        rho.Append(rho0);
        if(colors>1){
            mu.Append(mu1);
            rho.Append(rho1);}
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next())
            face_velocities(it.Full_Index())=analytic_velocity->u(it.Location()/m,0)(it.Axis())*m/s;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=0;
            T p=analytic_levelset->phi(it.Location()/m,0,c)*m;
            levelset_color.phi(it.index)=abs(p);
            levelset_color.color(it.index)=c==-4?bc_type:c;}

        if(test_analytic_diff){
            RANDOM_NUMBERS<T> rand;
            TV X=rand.Get_Uniform_Vector(grid.domain),dX;
            T e=1e-6,t=rand.Get_Uniform_Number(0,1);
            rand.Fill_Uniform(dX,-e,e);
            TV u0=analytic_velocity->u(X/m,t)*m/s,u1=analytic_velocity->u((X+dX)/m,t)*m/s;
            MATRIX<T,2> du0=analytic_velocity->du(X/m,t)/s,du1=analytic_velocity->du((X+dX)/m,t)/s;
            T erru=((du0+du1)*dX/2-(u1-u0)).Magnitude()/e;
            int c0,c1;
            T l0=analytic_levelset->phi(X/m,t,c0)*m,l1=analytic_levelset->phi((X+dX)/m,t,c1)*m;
            if(c0>=0) l0=-l0;
            if(c1>=0) l1=-l1;
            TV dl0=analytic_levelset->N(X/m,t),dl1=analytic_levelset->N((X+dX)/m,t);
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
    };

    struct ANALYTIC_VELOCITY_ROTATION:public ANALYTIC_VELOCITY
    {
        TV c;
        ANALYTIC_VELOCITY_ROTATION(TV cc): c(cc){}
        virtual TV u(const TV& X,T t) const {return (X-c).Orthogonal_Vector();}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>(0,1,-1,0);}
        virtual T p(const TV& X,T t) const {return 0;}
    };

    struct ANALYTIC_VELOCITY_VORTEX:public ANALYTIC_VELOCITY
    {
        T nu;
        ANALYTIC_VELOCITY_VORTEX(T mu,T pho): nu(mu/pho){}
        virtual TV u(const TV& X,T t) const
        {return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y))*exp(-2*nu*t);}
        virtual MATRIX<T,2> du(const TV& X,T t) const
        {T c=cos(X.x)*cos(X.y),s=sin(X.x)*sin(X.y);return MATRIX<T,2>(c,s,-s,-c)*exp(-2*nu*t);}
        virtual T p(const TV& X,T t) const {return 0;}
    };

    struct ANALYTIC_LEVELSET_PERIODIC:public ANALYTIC_LEVELSET
    {
        virtual T phi(const TV& X,T t,int& c) const {c=0;return 1;}
        virtual TV N(const TV& X,T t) const {return TV(1,0);}
    };

    struct ANALYTIC_LEVELSET_CIRCLE:public ANALYTIC_LEVELSET
    {
        TV cen;
        T r;
        ANALYTIC_LEVELSET_CIRCLE(TV cc,T rr): cen(cc),r(rr){}
        virtual T phi(const TV& X,T t,int& c) const {T p=(X-cen).Magnitude()-r;c=p>=0?-4:0;return abs(p);}
        virtual TV N(const TV& X,T t) const {return (X-cen).Normalized();}
    };

    struct ANALYTIC_LEVELSET_VORTEX:public ANALYTIC_LEVELSET
    {
        T k;
        ANALYTIC_LEVELSET_VORTEX(T kk): k(kk){}
        virtual T phi(const TV& X,T t,int& c) const
        {
            TV w=(X-pi/2).Normalized()+pi/2;
            VECTOR<T,3> z(w.x,w.y,0);
            for(int i=1;i<100;i++){
                T cx=cos(z.x),cy=cos(z.y),sx=sin(z.x),sy=sin(z.y);
                VECTOR<T,3> G(-2*X.x+2*z.x+z.z*cx*sy,-2*X.y+2*z.y+z.z*sx*cy,-k+sx*sy);
                MATRIX<T,3> H(2-z.z*sx*sy,z.z*cx*cy,cx*sy,z.z*cx*cy,2-z.z*sx*sy,sx*cy,cx*sy,sx*cy,0);
                z-=H.Solve_Linear_System(G);
                if(G.Magnitude_Squared()<1e-25) break;}
            c=((T).2-sin(X.x)*sin(X.y))>0?-4:0;
            return (z.Remove_Index(2)-X).Magnitude();
        }
        virtual TV N(const TV& X,T t) const
        {return TV(cos(X.x)*sin(X.y),sin(X.x)*cos(X.y)).Normalized();}
    };

    void Begin_Time_Step(const T time) PHYSBAM_OVERRIDE {}

    void End_Time_Step(const T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity){
            T max_error=0,a=0,b=0;
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                if(levelset_color.Color(it.Location())<0) continue;
                T A=face_velocities(it.Full_Index()),B=analytic_velocity->u(it.Location()/m,time)(it.Axis())*m/s;
                a=max(a,abs(A));
                b=max(b,abs(B));
                max_error=max(max_error,abs(A-B));
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(A-B));}
            LOG::cout<<"max_error "<<max_error<<"  "<<a<<"  "<<b<<std::endl;}
    }

    TV Dirichlet_Boundary_Condition(const TV& X,int bc_color,int fluid_color,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X/m,VECTOR<T,3>(1,0,0));
        if(analytic_velocity) return analytic_velocity->u(X/m,time)*m/s;
        return TV();
    }

    TV Neumann_Boundary_Condition(const TV& X,int bc_color,int fluid_color,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X/m,VECTOR<T,3>(0,1,0));
        if(analytic_velocity && analytic_levelset){
            MATRIX<T,2> du=analytic_velocity->du(X/m,time)/s;
            TV n=analytic_levelset->N(X/m,time);
            T p=analytic_velocity->p(X/m,time)*kg/(s*s*m);
            return (du+du.Transposed())*n*mu(fluid_color)+p*n;}
        return TV();
    }

//#####################################################################
};
}

#endif
