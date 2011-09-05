#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOAT
//#####################################################################
#ifndef __BOAT__
#define __BOAT__

#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Interpolation/BSPLINE.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_ROTATION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_RLE_3D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
namespace PhysBAM{

template<class T_input>
class BOAT:public SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_3D<T_input> >
{
    typedef T_input T;
    typedef RLE_GRID_3D<T> T_GRID;
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID> BASE;
    using BASE::fluids_parameters;using BASE::output_directory;using BASE::write_frame_title;using BASE::particle_levelset;using BASE::incompressible;
    using BASE::data_directory;using BASE::Adjust_Phi_With_Sources;using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::solids_parameters;
    using BASE::refine_all_water;using BASE::frame_rate;using BASE::first_frame;using BASE::last_frame;using BASE::initial_time;using BASE::use_incompressible_cfl;
    using BASE::restart;using BASE::restart_frame;using BASE::Time_At_Frame;using BASE::ground_j;using BASE::stream_type;

    struct FORCE_FROM_AGE:public NONLINEAR_FUNCTION<T(T)>
    {
        T initial_force,final_force,start_age,end_age;

        T operator()(const T age) const PHYSBAM_OVERRIDE
        {if(age<=start_age) return initial_force;
        if(age>=end_age) return final_force;
        T t=(age-start_age)/(end_age-start_age);
        return (1-t)*initial_force+t*final_force;}
    };

    T_GRID& grid;
    VORTEX_PARTICLE_EVOLUTION_RLE_3D<T> vortex_particle_evolution;
    T vorticity_magnitude;
    TV vorticity_direction;
    RANDOM_NUMBERS<T> random;
    T seeding_rate,unseeded_particles;
    FORCE_FROM_AGE force_from_age;

    T initial_height;
    ARRAY<BOX<TV> > sources;
    ARRAY<MATRIX<T,4> > world_to_source;
    ARRAY<TV> source_velocity;

    struct BOAT_PATH
    {
        T speed,amplitude,period,frequency;
        T period_distance;
        INTERPOLATION_CURVE<T,TV> position,velocity;

        BOAT_PATH(const T speed_input,const T amplitude_input,const T period_input)
            :speed(speed_input),amplitude(amplitude_input),period(period_input)
        {
            frequency=2*(T)pi/period;
            int n=1000;
            T start_time=0,end_time=0;
            for(int i=-5;i<=n+5;i++){
                T s=period*i/n;
                TV X(speed*s,0,amplitude*sin(frequency*s));
                T t=0;
                if(position.control_points.m){
                    T last_t=position.control_points.Last().t;
                    TV last_X=position.control_points.Last().value;
                    t=last_t+(last_X-X).Magnitude()/speed;}
                if(!i) start_time=t;else if(i==n) end_time=t;
                LOG::cout<<"i "<<i<<", X = "<<X<<std::endl;
                position.Add_Control_Point(t,X);}
            for(int i=1;i<=position.control_points.m;i++) position.control_points(i).t-=start_time;
            period_distance=speed*period;
            period=end_time-start_time;

            for(int i=1;i<=position.control_points.m;i++){
                T t=position.control_points(i).t;
                velocity.Add_Control_Point(t,position.Derivative(fmod(t+period,period)));}

            for(int i=0;i<=n;i++){
                T t=period*i/n;
                LOG::cout<<"i "<<i<<", t "<<t<<std::endl;
                LOG::cout<<"  t = "<<Position(t)<<std::endl;
                LOG::cout<<"  v = "<<Velocity(t)<<", s "<<Velocity(t).Magnitude()<<std::endl;
                LOG::cout<<"  r = "<<Orientation(t)<<std::endl;
                LOG::cout<<"  w = "<<Angular_Velocity(t)<<std::endl;}
        }

        TV Position(const T time) const
        {T t=fmod(time,period);
        return position.Value(t)+TV(period_distance*(time-t)/period,0,0);}

        TV Velocity(const T time) const
        {return velocity.Value(fmod(time,period));}

        ROTATION<TV> Orientation(const T time) const
        {T z=Position(time).z/amplitude;
        T z_dot=Velocity(time).z/Velocity(0).z;
        //TV direction=Velocity(time).Normalized();
        //ROTATION<TV> turn=ROTATION<TV>::Rotation_Quaternion(TV(1,0,0),direction);
        TV max_rotation=ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),Velocity(0).Normalized()).Rotation_Vector();
        ROTATION<TV> turn=ROTATION<TV>::From_Rotation_Vector(z_dot*max_rotation);
        ROTATION<TV> roll=ROTATION<TV>::From_Rotation_Vector(TV(-(T)pi/180*20*z,0,0));
        return turn*roll;}

        TV Angular_Velocity(const T time) const
        {T dt=(T)1e-3,one_over_dt=1/dt;
        return one_over_dt*(Orientation(time+dt)*Orientation(time).Inverse()).Rotation_Vector();}
    };

    FRAME<TV> boat_start;
    TV boat_velocity;
    RIGID_BODY<TV>* boat_rigid_body;
    BOAT_PATH boat_path;

    BOAT(const STREAM_TYPE stream_type,const int resolution,const T optical_depth)
        :SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_3D<T> >(stream_type),grid(*fluids_parameters.grid),vortex_particle_evolution(grid),initial_height((T).5),boat_rigid_body(0),boat_path(2,(T).5,2)
    {
        LOG::cout<<"running boat at resolution "<<resolution<<", optical depth "<<optical_depth<<std::endl;

        // frames
        frame_rate=96;
        first_frame=0;
        last_frame=2000;

        // set up the standard fluid environment
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_debug_data=true;
        fluids_parameters.write_ghost_values=true;fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.incompressible_iterations=40;
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.second_order_cut_cell_method=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_vorticity_confinement_fuel=false;

        // override standard fluid environment
        fluids_parameters.incompressible_iterations=200;

        // setup grid
        grid.Set_Uniform_Grid(GRID<TV>(150*resolution+1,25*resolution+1,30*resolution+1,BOX<TV>(0,15,0,(T)2.5,0,3))); // hopefully final domain
        //grid.Set_Uniform_Grid(GRID<TV>(75*resolution+1,25*resolution+1,30*resolution+1,BOX<TV>(0,(T)7.5,0,(T)2.5,0,3))); // full domain
        //grid.Set_Uniform_Grid(GRID<TV>(35*resolution+1,25*resolution+1,30*resolution+1,BOX<TV>(0,(T)3.5,0,(T)2.5,0,3))); // half domain
        //grid.Set_Uniform_Grid(GRID<TV>(15*resolution+1,25*resolution+1,10*resolution+1,BOX<TV>(0,(T)1.5,0,(T)2.5,1,2))); // tiny test domain
        //grid.Set_Uniform_Grid(GRID<TV>(15*resolution+1,25*resolution+1,30*resolution+1,BOX<TV>(0,(T)1.5,0,(T)2.5,0,3))); // steering test domain
        LOG::cout<<"uniform_grid = "<<grid.uniform_grid<<std::endl;

        // set bandwidths
        grid.Set_Positive_Bandwidth_In_Cells(3);
        if(optical_depth<5*grid.Minimum_Edge_Length()){
            LOG::cout<<"Warning: optical depth too small, switching to bandwidth of 5 cells."<<std::endl;
            grid.Set_Negative_Bandwidth_In_Cells(5);}
        else grid.Set_Negative_Bandwidth_From_Optical_Depth(optical_depth);
        grid.Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile();

        // set output directory
        std::string base_name="boat";
        output_directory=STRING_UTILITIES::string_sprintf("Boat/%s_resolution_%d_%d_%d_depth_%g",base_name.c_str(),(grid.uniform_grid.counts.x-1),(grid.uniform_grid.counts.y-1),(grid.uniform_grid.counts.z-1),optical_depth);
        if(refine_all_water) output_directory+="_uniform";

        RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

        // setup boat
        int boat=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/Thin_Shells/boat_hires",(T).5);
        boat_start.t=TV(-(T).5,(T).5,(T)1.5);
        boat_start.r=ROTATION<TV>((T).2,TV(0,0,1))*ROTATION<TV>(-(T)pi/2,TV(0,1,0));
        boat_velocity=TV(2,0,0);
        boat_rigid_body=&rigid_body_collection.Rigid_Body(boat);
        boat_rigid_body->Set_Frame(boat_start);
        boat_rigid_body->Update_Bounding_Box();
        boat_start.t.x+=(T).075-boat_rigid_body->Axis_Aligned_Bounding_Box().min_corner.x;
        Update_Boat(0);

        // add ground
        int ground=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/ground",(T).08);
        rigid_body_collection.rigid_body_particle.X(ground)=grid.uniform_grid.Domain().Center();
        rigid_body_collection.rigid_body_particle.X(ground).y=0;
        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Get_Collision_Geometry(ground)->refine_nearby_fluid=false;

        // nonstandard settings
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;
        fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;
        //particle_levelset.levelset.clamp_phi_with_collision_bodies=false;
        particle_levelset.use_removed_negative_particles_in_long_cells=true;
        fluids_parameters.collidable_phi_replacement_value=grid.positive_bandwidth;
        use_incompressible_cfl=false;

        // set coupling parameters
        fluids_parameters.Use_Fluid_Coupling_Defaults();

        // vortex particle settings
        fluids_parameters.use_body_force=true;
        vorticity_magnitude=400; // base
        //vorticity_magnitude=300;
        vorticity_direction=TV(1,0,0);
        vortex_particle_evolution.Set_Radius((T).05);
        vortex_particle_evolution.force_scaling=(T)5e-5; // base
        //vortex_particle_evolution.force_scaling=(T)3.75e-5;
        vortex_particle_evolution.renormalize_vorticity_after_stretching_tilting=true;
        seeding_rate=200;
        unseeded_particles=0;

        boat_rigid_body->Update_Bounding_Box();
        T boat_length_time=boat_rigid_body->Axis_Aligned_Bounding_Box().Edge_Lengths().x/boat_velocity.x;

        // vortex particle force scaling ramp
/*
        force_from_age.initial_force=(T)5e-5;
        force_from_age.final_force=(T)2.5e-5;
        force_from_age.start_age=(T).5*boat_length_time;force_from_age.end_age=(T)1.5*boat_length_time;
*/
        force_from_age.initial_force=(T)5e-5;
        force_from_age.final_force=(T)1.25e-5;
        force_from_age.start_age=(T).5*boat_length_time;force_from_age.end_age=(T)2*boat_length_time;
        LOG::cout<<"force scaling ramp"<<std::endl;
        LOG::cout<<"  force = "<<force_from_age.initial_force<<" to "<<force_from_age.final_force<<std::endl;
        LOG::cout<<"  age= "<<force_from_age.start_age<<" to "<<force_from_age.end_age<<std::endl;
        vortex_particle_evolution.Set_Force_Scaling_From_Age(force_from_age);

        fluids_parameters.solve_neumann_regions=false;
    }

    ~BOAT()
    {}

    // Unused callbacks
    void Initialize_Velocities() PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Modify_Grid_After_Rebuild(T_GRID& new_grid,const T time) PHYSBAM_OVERRIDE {}
    void Transfer_Extra_State(const T_GRID& new_grid) PHYSBAM_OVERRIDE {}

    T Initial_Ground(const TV_HORIZONTAL& X) const PHYSBAM_OVERRIDE
    {return 0;}

    T Initial_Phi(const TV& X) const PHYSBAM_OVERRIDE
    {T phi=X.y-initial_height;
    phi=max(phi,-boat_rigid_body->Implicit_Geometry_Extended_Value(X));
    return phi;}

    T Initial_Phi_Object(const TV& X) const PHYSBAM_OVERRIDE
    {return 1;}

    void Initialize_Grids_Extra() PHYSBAM_OVERRIDE
    {vortex_particle_evolution.Initialize_Grid();}

    void Update_Boat(const T time)
    {bool simple_path=false;
    if(simple_path){
        boat_rigid_body->X()=boat_start.t+time*boat_velocity;
        boat_rigid_body->Twist().linear=boat_velocity;}
    else{
        boat_rigid_body->X()=boat_path.Position(time)+boat_start.t;
        boat_rigid_body->Rotation()=boat_path.Orientation(time)*boat_start.r;
        boat_rigid_body->Twist().linear=boat_path.Velocity(time);
        boat_rigid_body->Twist().angular=boat_path.Angular_Velocity(time);}}

    void Construct_Levelsets_For_Objects(const T time)
    {Update_Boat(time);
    LOG::cout<<"boat position = "<<boat_rigid_body->X()<<", velocity = "<<boat_rigid_body->velocity<<std::endl;
    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Construct_Levelsets_For_Objects(time);}

    void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Update_Fluid_Parameters(dt,time);/*tests->Update_Sources(time);*/}

    bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
    {if(time==initial_time) return false; // skip initial call
    for(int s=1;s<=sources.m;s++)Adjust_Phi_With_Source(sources(s),world_to_source(s),true);return false;}

    void Get_Source_Reseed_Mask(ARRAY<bool>*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
    {for(int s=1;s<=sources.m;s++)Get_Source_Reseed_Mask(sources(s),world_to_source(s),cell_centered_mask,s==1);}

    void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
    {for(int s=1;s<=sources.m;s++)Get_Source_Velocities(sources(s),world_to_source(s),source_velocity(s));}

    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
    {T rigid_dt_denominator=TV::Dot_Product(abs(boat_velocity),grid.uniform_grid.One_Over_DX());
    dt=min(dt,1/rigid_dt_denominator);}

    void Get_Cell_Should_Be_Long(ARRAY<bool>& cell_should_be_long,const T time) const PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Get_Cell_Should_Be_Long(cell_should_be_long,time);
    for(int s=1;s<=sources.m;s++)Get_Source_Cell_Should_Be_Long(sources(s),world_to_source(s),cell_should_be_long,true);}

    struct RIGID_BODY_IMPLICIT_SURFACE:public IMPLICIT_OBJECT<TV>
    {
        const RIGID_BODY<TV>& rigid_body;
        RIGID_BODY_IMPLICIT_SURFACE(const RIGID_BODY<TV>& rigid_body_input)
            :rigid_body(rigid_body_input)
        {}

        T Extended_Phi(const TV& X) const
        {return rigid_body.Implicit_Geometry_Extended_Value(X);}
    };

    void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Read_Output_Files_Fluids(frame);
    vortex_particle_evolution.Read_Output_Files(stream_type,output_directory,frame);}

    void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Write_Output_Files(frame);
    vortex_particle_evolution.Write_Output_Files(stream_type,output_directory,frame);}

//#####################################################################
    void Get_Body_Force(ARRAY<T>& force,const T dt,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Get_Body_Force
//#####################################################################
template<class T> void BOAT<T>::
Get_Body_Force(ARRAY<T>& force,const T dt,const T time)
{
    LOG::SCOPE scope("VORTEX","evolving vortex particles");
    VORTICITY_PARTICLES<TV>& vorticity_particles=vortex_particle_evolution.vorticity_particles;
    ARRAY_VIEW<T>& vorticity_particle_age=*vorticity_particles.array_collection->template Get_Array<T >(ATTRIBUTE_ID_AGE);

    ARRAY<T> V_ghost(grid.number_of_faces,false);
    incompressible.boundary->Fill_Ghost_Cells_Face(grid,incompressible.V,V_ghost,time);

    // compute force and advance particles
    force.Fill(0);
    vortex_particle_evolution.Compute_Body_Force(V_ghost,force,dt,time);
    vortex_particle_evolution.Euler_Step(V_ghost,dt,time);
    vorticity_particle_age+=dt;

    // delete particles that are outside the water or too close to the true domain boundary
    BOX<TV> safe_domain=grid.mpi_grid?grid.mpi_grid->global_uniform_grid.Domain():grid.domain;
    safe_domain.Change_Size(-vortex_particle_evolution.radius);
    for(int p=vorticity_particles.array_collection->Size();p>=1;p--)if(safe_domain.Lazy_Outside(vorticity_particles.X(p)) || particle_levelset.levelset.Phi(vorticity_particles.X(p))>0)
        vorticity_particles.array_collection->Delete_Element(p);

/*
    static bool first_time=true;
    TV center=safe_domain.Center();center.x+=grid.uniform_grid.dx/2;center.z+=grid.uniform_grid.dx/2;center.y=initial_height-(T).05;
    if(first_time){
        int p=vorticity_particles.array_collection->Add_Element();
        vorticity_particles.vorticity(p)=vorticity_magnitude*TV(0,1,0);}
    if(vorticity_particles.array_collection->Size()) vorticity_particles.X(1)=center;
    if(vorticity_particles.array_collection->Size()) vorticity_particles.vorticity(1)=vorticity_magnitude*TV(0,1,0);
    first_time=false;
*/

    // figure out where to seed
    BOX<TV> motor(-(T).4,-(T).3,-(T).2,-(T).05,-(T).2,(T).2);

    // seed particles with rejection sampling...to get a location
    unseeded_particles+=seeding_rate*dt;    
    int iterations=1000;
    while(unseeded_particles>=1){
        if(!iterations--){LOG::cout<<"ran out of seeding iterations"<<std::endl;break;}
        TV X=boat_rigid_body->X()+(boat_rigid_body->Rotation()*boat_start.r.Inverse()).Rotate(random.Get_Uniform_Vector(motor));
        if(particle_levelset.levelset.Phi(X)<0 && boat_rigid_body->Implicit_Geometry_Extended_Value(X)>0){
            int p=vorticity_particles.array_collection->Add_Element();
            vorticity_particles.X(p)=X;
            int sign=random.Get_Uniform_Integer(0,1);if(sign==0)sign=-1;
            //vorticity_particles.vorticity(p)=vorticity_magnitude*vorticity_direction;
            vorticity_particles.vorticity(p)=vorticity_magnitude*random.template Get_Direction<TV>();
            vorticity_particles.vorticity(p)=vorticity_magnitude*random.template Get_Direction<TV>();
            vorticity_particle_age(p)=0;
            unseeded_particles--;
            LOG::cout<<"added vortex particle, X = "<<X<<std::endl;}}
    LOG::cout<<"total vortex particles = "<<vortex_particle_evolution.vorticity_particles.array_collection->Size()<<std::endl;
}
//#####################################################################
}
#endif
#endif
