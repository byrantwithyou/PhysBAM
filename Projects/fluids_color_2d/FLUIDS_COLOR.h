//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR__
#define __FLUIDS_COLOR__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
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

        output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number);
        switch(test_number){
            case 0: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);break;
            case 1: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);break;
            case 2: grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);break;
            default: PHYSBAM_FATAL_ERROR("Missing test number");}
    }

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);}

    void Initialize()
    {
        switch(test_number){
            case 0: Uniform_Translation_Periodic();break;
            case 1: Taylor_Green_Vortex_Periodic();break;
            case 2: Taylor_Green_Vortex_BC((T).2);break;
            default: PHYSBAM_FATAL_ERROR("Missing test number");}

        if(user_last_frame) last_frame=stored_last_frame;
    }

    void Uniform_Translation_Periodic()
    {
        mu.Append(mu0);
        rho.Append(rho0);
        face_velocities.Fill(1);
        levelset_color.phi.Fill(-1);
        levelset_color.color.Fill(0);
    }

    TV Taylor_Green_Vortex_Velocity(const TV& X,T t) const
    {
        return TV(sin(X.x/m)*cos(X.y/m),-cos(X.x/m)*sin(X.y/m))*exp(-2*mu0/rho0*t/(m*m))*m/s;
    }

    T Taylor_Green_Vortex_Isocontour(const TV& X) const
    {
        return sin(X.x/m)*sin(X.y/m);
    }

    TV Taylor_Green_Vortex_Normal(const TV& X) const
    {
        return Taylor_Green_Vortex_Velocity(X,0).Orthogonal_Vector();
    }

    T Taylor_Green_Vortex_Levelset(const TV& X,T k) const
    {
        TV w=(X-pi/2).Normalized()+pi/2;
        VECTOR<T,3> z(w.x,w.y,0);
        for(int i=1;i<100;i++){
            T cx=cos(z.x),cy=cos(z.y),sx=sin(z.x),sy=sin(z.y);
            VECTOR<T,3> G(-2*X.x+2*z.x+z.z*cx*sy,-2*X.y+2*z.y+z.z*sx*cy,-k+sx*sy);
            MATRIX<T,3> H(2-z.z*sx*sy,z.z*cx*cy,cx*sy,z.z*cx*cy,2-z.z*sx*sy,sx*cy,cx*sy,sx*cy,0);
            z-=H.Solve_Linear_System(G);
            if(G.Magnitude_Squared()<1e-25) break;}
        T sign=((T).2-sin(X.x)*sin(X.y))>0?1:-1;
        return (z.Remove_Index(2)-X).Magnitude()*sign;
    }

    void Taylor_Green_Vortex_Periodic()
    {
        mu.Append(mu0);
        rho.Append(rho0);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next())
            face_velocities(it.Full_Index())=Taylor_Green_Vortex_Velocity(it.Location(),0)(it.Axis());
        levelset_color.phi.Fill(-1);
        levelset_color.color.Fill(0);
    }

    void Taylor_Green_Vortex_BC(T contour)
    {
        mu.Append(mu0);
        rho.Append(rho0);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next())
            face_velocities(it.Full_Index())=Taylor_Green_Vortex_Velocity(it.Location(),0)(it.Axis());
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            T p=Taylor_Green_Vortex_Levelset(it.Location(),contour);
            levelset_color.phi(it.index)=abs(p);
            levelset_color.color(it.index)=p>0?bc_type:0;}
    }

    void Begin_Time_Step(const T time) PHYSBAM_OVERRIDE {}

    void End_Time_Step(const T time) PHYSBAM_OVERRIDE
    {
        if(test_number==1 || test_number==2){
            T max_error=0,a=0,b=0;
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                if(levelset_color.Color(it.Location())<0) continue;
                T A=face_velocities(it.Full_Index()),B=Taylor_Green_Vortex_Velocity(it.Location(),time)(it.Axis());
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
        if(test_number==2) return Taylor_Green_Vortex_Velocity(X,time);
        return TV();
    }

    TV Neumann_Boundary_Condition(const TV& X,int bc_color,int fluid_color,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));
        if(test_number==2){
            TV N=Taylor_Green_Vortex_Normal(X);
            N.y*=-1;
            return 2*sin(X.x/m)*sin(X.y/m)*N;}
        return TV();
    }

//#####################################################################
};
}

#endif
