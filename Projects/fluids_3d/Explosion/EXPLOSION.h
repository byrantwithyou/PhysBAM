//#####################################################################
// Copyright 2002-2004, Ron Fedkiw, Frank Losasso, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXPLOSION
//##################################################################### 
#ifndef __EXPLOSION__
#define __EXPLOSION__

#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template <class T,class RW=T>
class EXPLOSION:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID<TV>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::last_frame;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::abort_when_dt_below;

    T source_vorticity_magnitude;
    T particle_vorticity_minimum_density;
    RANDOM_NR3 random;
    VORTEX_PARTICLE_EVOLUTION_3D<T,RW> vortex_particle_evolution;

    T rho,rho_bottom,rho_top,buoyancy_constant,thermal_buoyancy_constant;T explosion_divergence;T explosion_end_time;
    BOX_3D<T> source_domain;
    bool use_vortex_particles;
    T particle_radius;

    EXPLOSION()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,fluids_parameters.SMOKE)
    {
        //fluids_parameters.grid.Initialize(51,76,51,0,1,0,(T)1.5,0,1);
        PARAMETER_LIST parameters;parameters.Read("sim.params");
//        VECTOR<int,3> nodes=parameters.Get_Parameter_And_Set_If_Undefined("nodes",VECTOR<int,3>(51,76,51));
//        VECTOR<int,3> nodes=parameters.Get_Parameter_And_Set_If_Undefined("nodes",VECTOR<int,3>(51,76,51));
        VECTOR<int,3> nodes=parameters.Get_Parameter_And_Set_If_Undefined("nodes",VECTOR<int,3>(101,151,101));
        fluids_parameters.use_maccormack_semi_lagrangian_advection=parameters.Get_Parameter_And_Set_If_Undefined("use_maccormack_semi_lagrangian_advection",false);
        fluids_parameters.grid.Initialize(nodes.x,nodes.y,nodes.z,0,1,0,1.5,0,1);
        first_frame=0;
        last_frame=parameters.Get_Parameter_And_Set_If_Undefined("last_frame",(int)500);
        restart_frame=parameters.Get_Parameter_And_Set_If_Undefined("restart_frame",(int)0);
        restart=parameters.Get_Parameter_And_Set_If_Undefined("restart",(bool)false);
        frame_rate=60;
        fluids_parameters.cfl=.9;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.use_vorticity_confinement=parameters.Get_Parameter_And_Set_If_Undefined("use_vorticity_confinement",false);
        fluids_parameters.confinement_parameter=parameters.Get_Parameter_And_Set_If_Undefined("confinement_parameter",(T).25);
        fluids_parameters.kolmogorov=(T)0;
        rho=5;rho_bottom=1;rho_top=(T).65;buoyancy_constant=0;fluids_parameters.gravity=0;
        fluids_parameters.temperature_products=3000;
        fluids_parameters.write_velocity=true;fluids_parameters.use_body_force=true;
        fluids_parameters.write_debug_data=true;//solids_parameters.write=false;solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=false;
        source_domain=BOX_3D<T>((T).45,(T).55,(T)0,(T).1,(T).45,(T).55);
        fluids_parameters.Initialize_Domain_Boundary_Conditions();
        output_directory=parameters.Get_Parameter_And_Set_If_Undefined("output_directory",(std::string)"Explosion/output");

        explosion_divergence=100;
        explosion_end_time=(T).5;
        fluids_parameters.use_non_zero_divergence=true;
        thermal_buoyancy_constant=.001;

        // setup vortex particles
        fluids_parameters.use_body_force=true;
        use_vortex_particles=parameters.Get_Parameter_And_Set_If_Undefined("use_vortex_particles",false);
        vortex_particle_evolution.particle_confinement_parameter=(T)1.00;
        vortex_particle_evolution.renormalize_vorticity_after_stretching_tilting=true;
        particle_vorticity_minimum_density=(T).9;
        source_vorticity_magnitude=(T)200;
        particle_radius=(T).02;
        vortex_particle_evolution.Initialize(fluids_parameters.grid);
        vortex_particle_evolution.force_scaling=parameters.Get_Parameter_And_Set_If_Undefined("vortex_force_scale",(T)1e-4);

        parameters.Write("ran.param");
    }

    ~EXPLOSION()
    {}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluids_parameters.incompressible->projection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(T)0; 
    // tests.Initial_Velocity(iterator.Location())[iterator.Axis()];
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    if(time<=explosion_end_time) for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++) if(source_domain.Lazy_Inside(grid.X(i,j,ij))){
            fluids_parameters.density_container.density(i,j,ij)=rho;fluids_parameters.temperature_container.temperature(i,j,ij)=fluids_parameters.temperature_products;}
    // keep density >= 0 and T >=0
    for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++){
        fluids_parameters.density_container.density(i,j,ij)=max((T)0,fluids_parameters.density_container.density(i,j,ij));
        fluids_parameters.temperature_container.temperature(i,j,ij)=max((T)fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature(i,j,ij));}
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Get_Divergence(ARRAY<T,VECTOR<int,3> >& divergence,const T dt,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    T expansion=2*explosion_divergence*sin(2*time/explosion_end_time)/exp(2*time/explosion_end_time);
    for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++) if(source_domain.Lazy_Inside(grid.X(i,j,ij))) divergence(i,j,ij)=expansion;
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
void Get_Body_Force(T_FACE_ARRAYS_SCALAR& force,const T dt,const T time) PHYSBAM_OVERRIDE
{ 
    LOG::Push_Scope("Vortex","Vortex Particle Solve");
    GRID<TV>& grid=fluids_parameters.grid;
    T_FACE_ARRAYS_SCALAR::copy(T(),force);
    // Buoyancy
    DENSITY_CONTAINER<GRID<TV> >& density=fluids_parameters.density_container;
    TEMPERATURE_CONTAINER<GRID<TV> >& temperature=fluids_parameters.temperature_container;

    ARRAY<T,VECTOR<int,3> > density_ghost(grid,3,false),temperature_ghost(grid,3,false);
    density.boundary->Fill_Ghost_Cells_Cell(grid,density.density,density_ghost,time);
    temperature.boundary->Fill_Ghost_Cells_Cell(grid,temperature.temperature,temperature_ghost,time);

    GRID<TV> v_grid=grid.Get_Y_Face_Grid();
    for(int i=1;i<=v_grid.m;i++) for(int j=1;j<=v_grid.n;j++) for(int ij=1;ij<=v_grid.mn;ij++){ // y-direction forces only
        T rho_atm=rho_bottom+(rho_top-rho_bottom)*(grid.y(j)-grid.ymin)/(grid.ymax-grid.ymin);
        if(density_ghost(i,j,ij)>.05){
            T density_difference=(T).5*(density_ghost(i,j-1,ij)+density_ghost(i,j,ij))-rho_atm;
            T temperature_difference=(T).5*(temperature_ghost(i,j-1,ij)+temperature_ghost(i,j,ij))-fluids_parameters.temperature_container.ambient_temperature;
            if(density_difference>0||temperature_difference>0)force.v(i,j,ij)=thermal_buoyancy_constant*temperature_difference-buoyancy_constant*density_difference;}}

    if(use_vortex_particles){
        // Vortex particles
        T_FACE_ARRAYS_SCALAR face_velocities_ghost(grid,3,false);
        fluids_parameters.incompressible->boundary->Fill_Ghost_Cells_Face(grid,fluids_parameters.incompressible->projection.face_velocities,face_velocities_ghost,time);
        // Compute force
        vortex_particle_evolution.Compute_Body_Force(face_velocities_ghost,force,dt,time);
        // Add particles
        if(time <= explosion_end_time){
            int add_count=0;VORTICITY_PARTICLES<T,VECTOR<T,3> >& vorticity_particles=vortex_particle_evolution.vorticity_particles;
            VECTOR<T,3> cell_upper=(T).5*VECTOR<T,3>(grid.dx,grid.dy,grid.dz),cell_lower=-cell_upper;
            for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++)
                if(source_domain.Lazy_Inside(grid.X(i,j,ij)) && time>(T)1/24 && random.Get_Uniform_Number((T)0,(T)1)<(T).005){
                    add_count++;
                    int particle_id=vorticity_particles.array_collection->Add_Element(); 
                    vorticity_particles.radius(particle_id)=particle_radius;
                    vorticity_particles.X(particle_id)=grid.X(i,j,ij)+random.Get_Uniform_Vector(cell_lower,cell_upper);
                    vorticity_particles.vorticity(particle_id)=(T)source_vorticity_magnitude*VECTOR<T,3>::Cross_Product(VECTOR<T,3>(0,1,0),(vorticity_particles.X(particle_id)-source_domain.Center()).Normalized()).Normalized();}
            printf("Added %d particles\n",add_count);}
        // Advance Particles
        vortex_particle_evolution.Euler_Step(face_velocities_ghost,dt,time);}
    LOG::Pop_Scope();
}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Read_Output_Files_Fluids(frame);
    if(use_vortex_particles) vortex_particle_evolution.Read_Output_Files(output_directory,frame);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Write_Output_Files(frame);
    if(use_vortex_particles) vortex_particle_evolution.Write_Output_Files(output_directory,frame);
}
//#####################################################################
};    
}
#endif


