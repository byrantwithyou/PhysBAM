//#####################################################################
// Copyright 2006-2007, Nipun Kwatra, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPRAY_AND_FOAM
//#####################################################################
#ifndef __SPRAY_AND_FOAM__
#define __SPRAY_AND_FOAM__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_PARTICLES.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <cstdio>
namespace PhysBAM{

template<class T_input,class T_GRID=GRID<VECTOR<T_input,3> > >
class SPRAY_AND_FOAM:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>
{
    typedef T_input T;
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;using BASE::output_directory;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::stream_type;using BASE::Get_Object_Velocities;using BASE::Adjust_Particle_For_Objects; // silence -Woverloaded-virtual
    using BASE::solid_body_collection;using BASE::test_number;using BASE::resolution;

    bool initial_frame;
    std::string input_water_directory;
    T interpolated_water_time;
    T_ARRAYS_SCALAR interpolated_water_phi;
    T_FACE_ARRAYS_SCALAR interpolated_water_V;

    int current_low_frame;
    GRID<TV> low_frame_grid;
    T_ARRAYS_SCALAR low_frame_phi;
    T_FACE_ARRAYS_SCALAR low_frame_V;

    int current_high_frame;
    GRID<TV> high_frame_grid;
    T_ARRAYS_SCALAR high_frame_phi;
    T_FACE_ARRAYS_SCALAR high_frame_V;

    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> last_frame_removed_particles;
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> new_removed_particles;
    ARRAY<T> new_removed_particles_spawn_time;
    T time_particles_were_spawned_last;

    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> spray_particles;
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> foam_particles;

    bool use_log_removed_particle_times;
    int particle_number_amplification;
    T life_foam;

    T_ARRAYS_SCALAR object_phi;
    T_FACE_ARRAYS_SCALAR object_V;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    int lighthouse,beach,cove;

    RANGE<VECTOR<T,3> > kill_particle_box;
    RANDOM_NUMBERS<T> random;

    SPRAY_AND_FOAM(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>(stream_type,0,fluids_parameters.SMOKE),rigid_body_collection(solid_body_collection.rigid_body_collection)
    {
    }

    virtual ~SPRAY_AND_FOAM()
    {}


//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    // CALL THIS BEFORE THE CLASS IS CONSTRUCTED
    //PARTICLE_LEVELSET_REMOVED_PARTICLES<T,TV>::Attribute_Map().Store(PARTICLE_LEVELSET_REMOVED_PARTICLES<T,TV>::Attribute_Map().id);
//        fluids_parameters.grid->Initialize(97,33,65,0,60,0,20,0,40);
//        fluids_parameters.grid->Initialize(151,51,101,0,60,0,20,0,40);
//        fluids_parameters.grid->Initialize(240,40,160,-60,60,0,20,0,80);
    int cells=1*resolution;
    fluids_parameters.grid->Initialize(14*cells+1,3*cells+1,8*cells+1,-80,60,0,30,0,80,true);
//        fluids_parameters.grid->Initialize(385,129,129,0,60,0,20,0,20);
    first_frame=0;
    last_frame=3000;
    frame_rate=36;
    fluids_parameters.cfl=3.5;
    fluids_parameters.use_density=true;
    fluids_parameters.use_temperature=false;
    fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[2][1]=false;fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T)10;
    fluids_parameters.gravity=(T)0;
    write_output_files=true;fluids_parameters.write_debug_data=true;
    fluids_parameters.fluid_boundary=new BOUNDARY_OPEN_WATER<T_GRID>((T).5,!fluids_parameters.domain_walls[0][0],!fluids_parameters.domain_walls[0][1],!fluids_parameters.domain_walls[1][0],
        !fluids_parameters.domain_walls[1][1],!fluids_parameters.domain_walls[2][0],!fluids_parameters.domain_walls[2][1]);
    fluids_parameters.Initialize_Domain_Boundary_Conditions();

    fluids_parameters.kolmogorov=(T).1;

    initial_frame=true;
    output_directory=STRING_UTILITIES::string_sprintf("Spray_And_Foam/Secondary_Resolution_%d_%d_%d",(fluids_parameters.grid->counts.x-1),(fluids_parameters.grid->counts.y-1),(fluids_parameters.grid->counts.z-1));

    //input_water_directory="/disk2/spray_backup/lowres_lighthouse_sph/Test_3_Lighthouse_Resolution_240_40_160/";
    //input_water_directory="/solver/vol3/spray/merges/Test_2_Lighthouse_Resolution_560_120_320/";
    input_water_directory="/n/redgreen/disk2/spray/Lighthouse_Hires_TWC_UD_50FIXEDPIC_EO_25BP_MERGED";
    interpolated_water_phi.Resize(fluids_parameters.grid->Domain_Indices(3));
    interpolated_water_V.Resize(*fluids_parameters.grid,3);

    current_low_frame=INT_MAX;current_high_frame=INT_MAX;

    use_log_removed_particle_times=true;
    time_particles_were_spawned_last=0;
    particle_number_amplification=10; // ~100,000 particles per amplification
    life_foam=5;
    kill_particle_box=RANGE<VECTOR<T,3> >(58,61,0,30,0,80);
    random.Set_Seed(1);

    /*
      lighthouse=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/lighthouse",(T)2,true,true,false);
      rigid_body_collection.Rigid_Body(lighthouse).frame.t=TV((T)25,(T)-2,(T)30);rigid_body_collection.Rigid_Body(lighthouse).is_kinematic=true;
      beach=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/icecube",(T)100,true,true,false);
      rigid_body_collection.Rigid_Body(beach).frame.t=TV((T)-25,(T)-41,(T)35);rigid_body_collection.Rigid_Body(beach).is_kinematic=true;
      rigid_body_collection.Rigid_Body(beach).frame.r=QUATERNION<T>::From_Euler_Angles(0,0,-.1);
    */
    lighthouse=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/lighthouse",(T)2,true,true,false);
    rigid_body_collection.rigid_body_particle.frame(lighthouse).t=TV((T)25,(T)-2,(T)30);rigid_body_collection.Rigid_Body(lighthouse).Is_Kinematic()=true;
    cove=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/cove",(T)1,true,true,false);
    rigid_body_collection.rigid_body_particle.frame(cove).t=TV((T)-80,(T)11,(T)40);rigid_body_collection.Rigid_Body(cove).Is_Kinematic()=true;
    /*lighthouse=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/lighthouse",(T)2,true,true,false);
      cove=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/cove",(T)1,true,true,false);
      rigid_body_collection.rigid_body_particle(cove)->Frame().t=TV((T)-80,(T)11,(T)40);rigid_body_collection.rigid_body_particle(cove)->is_kinematic=true;
    */
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);
    object_phi.Resize(fluids_parameters.grid->Domain_Indices(3));object_V.Resize(*fluids_parameters.grid,3);
    for(int i=object_phi.domain.min_corner.x;i<=object_phi.domain.max_corner.x;i++)for(int j=object_phi.domain.min_corner.y;j<=object_phi.domain.max_corner.y;j++)for(int ij=object_phi.domain.min_corner.z;ij<=object_phi.domain.max_corner.z;ij++){
        TV X=fluids_parameters.grid->X(i,j,ij);
        object_phi(i,j,ij)=min(rigid_body_collection.Rigid_Body(lighthouse).Implicit_Geometry_Extended_Value(X),rigid_body_collection.Rigid_Body(cove).Implicit_Geometry_Extended_Value(X));
        //object_phi(i,j,ij)=min(rigid_body_collection.rigid_body_particle(lighthouse)->Implicit_Geometry_Extended_Value(X),rigid_body_collection.rigid_body_particle(beach)->Implicit_Geometry_Extended_Value(X));
    }
    LEVELSET_3D<GRID<TV> > object_levelset(*fluids_parameters.grid,object_phi);
    object_levelset.Fast_Marching_Method();

    Load_Simulated_Water(first_frame);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Limit_Dt
//#####################################################################
/*void Limit_Dt(T& dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::Limit_Dt(dt,time);
    dt=min(dt,(T)1/(5*frame_rate));
}*/
//#####################################################################
// Function Get_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::Set_Dirichlet_Boundary_Conditions(time);
    T_FACE_ARRAYS_BOOL& psi_N=fluids_parameters.incompressible->projection.elliptic_solver->psi_N;
    T_FACE_ARRAYS_SCALAR& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
        TV_INT cell1=iterator.First_Cell_Index();TV_INT cell2=iterator.Second_Cell_Index();
        if(interpolated_water_phi(cell1)<0||interpolated_water_phi(cell2)<0){
            psi_N.Component(axis)(face)=true;face_velocities(axis,face)=interpolated_water_V(axis,face);}}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::Preprocess_Frame(frame);

    Load_Simulated_Water(frame);

    Find_Removed_Particles_That_Will_Appear_This_Frame(-1,frame);
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    // keep density >= 0 and T >=0
    for(int i=0;i<fluids_parameters.grid->counts.x;i++)for(int j=0;j<fluids_parameters.grid->counts.y;j++)for(int ij=0;ij<fluids_parameters.grid->counts.z;ij++){
        if(interpolated_water_phi(i,j,ij)<=0)
            fluids_parameters.density_container.density(i,j,ij)=0;
        else
            fluids_parameters.density_container.density(i,j,ij)=max((T)0,fluids_parameters.density_container.density(i,j,ij));}
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time) PHYSBAM_OVERRIDE
{
    Compute_Interpolated_Water(time);
    PHYSBAM_FATAL_ERROR();
#if 0    
    fluids_parameters.Extrapolate_Velocity_Into_Object(interpolated_water_phi,interpolated_water_V,3,true,time,false);
#endif
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
        if(axis==1) fluid_collection.incompressible_fluid_collection.face_velocities(axis,face)=10;
        else fluid_collection.incompressible_fluid_collection.face_velocities(axis,face)=0;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    RANGE<VECTOR<T,3> > source_domain=RANGE<VECTOR<T,3> >(-61,-57,-1,21,-1,81);
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        TV_INT face=iterator.Face_Index();int axis=iterator.Axis();TV X=fluids_parameters.grid->Face(axis,face);
        if(source_domain.Lazy_Inside(X)){
            //fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(face)=true;
            if(axis==1) face_velocities(axis,face)=10;
            else face_velocities(axis,face)=0;}}
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(initial_frame){initial_frame=false;return;}
    static int unique_id=1;

    //LOG::Push_Scope("PART UPDATE","Spray/Foam Particle Update");

    Compute_Interpolated_Water(time);
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::Update_Fluid_Parameters(dt,time);

/*
    RANGE<VECTOR<T,3> > source_domain(15,17,15,17,19,21);
    for(int i=0;i<fluids_parameters.grid->m;i++)for(int j=0;j<fluids_parameters.grid->n;j++)for(int ij=0;ij<fluids_parameters.grid->mn;ij++)if(source_domain.Lazy_Inside(fluids_parameters.grid->X(i,j,ij))){
        fluids_parameters.density_container.density(i,j,ij)=1;}
*/

    T density_conversion_speed=(T).29;
    LOG::Time(STRING_UTILITIES::string_sprintf("adding density from %d spray particles",spray_particles.array_collection->Size()));
    for(int k=0;k<spray_particles.array_collection->Size();k++){
        T b=dt*density_conversion_speed;
        T_GRID &grid=*fluids_parameters.grid;
        T one_over_total_particle_cell_weight=0;
        TV radius_vector=(spray_particles.radius(k)+grid.Minimum_Edge_Length())*TV::All_Ones_Vector();
        T one_over_radius_squared=1/sqr(spray_particles.radius(k));
        RANGE<TV_INT> particle_cells(grid.Clamp_To_Cell(spray_particles.X(k)-radius_vector),grid.Clamp_To_Cell(spray_particles.X(k)+radius_vector));

        for(CELL_ITERATOR iterator(grid,particle_cells);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();
            TV X_minus_Xp=grid.Center(cell)-spray_particles.X(k);T distance_squared=X_minus_Xp.Magnitude_Squared();
            T weight=1-distance_squared*one_over_radius_squared;
            if(weight>0) one_over_total_particle_cell_weight+=weight;}
        if(one_over_total_particle_cell_weight) one_over_total_particle_cell_weight=1/one_over_total_particle_cell_weight;
        for(CELL_ITERATOR iterator(grid,particle_cells);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();
            TV X_minus_Xp=grid.Center(cell)-spray_particles.X(k);T distance_squared=X_minus_Xp.Magnitude_Squared();
            T weight=1-distance_squared*one_over_radius_squared;
            if(weight>0) fluids_parameters.density_container.density(cell)+=b*weight*one_over_total_particle_cell_weight;}
        if(!one_over_total_particle_cell_weight){
            TV_INT i=fluids_parameters.grid->Closest_Node(spray_particles.X(k));
            fluids_parameters.density_container.density(i)+=b;}

        spray_particles.radius(k)-=b;
        if(spray_particles.radius(k)<=0)
            spray_particles.array_collection->Add_To_Deletion_List(k);}
    LOG::Time("deleting spray particles on deletion list");
    spray_particles.array_collection->Delete_Elements_On_Deletion_List();
    LOG::cout<<"kept "<<spray_particles.array_collection->Size()<<" spray particles"<<std::endl;

    // spawn new particles
    LOG::Time("spawning new particles");
    int number_of_particles_to_add=0;
    T particle_velocity_low=(T).1;
    T particle_velocity_high=(T)9.5;
    int particle_number_amplification_low=0;
    int particle_number_amplification_high=7;
    ARRAY<int> new_removed_particles_amplification_factor(new_removed_particles.array_collection->Size());
    for(int i=0;i<new_removed_particles.array_collection->Size();i++)
        if(new_removed_particles_spawn_time(i)<=time+dt&&new_removed_particles_spawn_time(i)>time && !kill_particle_box.Lazy_Inside(new_removed_particles.X(i))){
            T particle_velcocity_magnitude=new_removed_particles.V(i).Magnitude();
            particle_velcocity_magnitude=clamp(particle_velcocity_magnitude,particle_velocity_low,particle_velocity_high);
            T particle_velocity_scale=(particle_velcocity_magnitude-particle_velocity_low)/(particle_velocity_high-particle_velocity_low);
            int amplification_factor=particle_number_amplification_low+(int)(particle_velocity_scale*(particle_number_amplification_high-particle_number_amplification_low));
            new_removed_particles_amplification_factor(i)=amplification_factor;
            number_of_particles_to_add+=amplification_factor;}

    int particle_to_add=spray_particles.array_collection->Size()+1;
    LOG::cout<<"spawning "<<number_of_particles_to_add<<" new particles"<<"of "<<new_removed_particles.array_collection->Size()*particle_number_amplification<<std::endl;
    spray_particles.array_collection->Add_Elements(number_of_particles_to_add);
    ARRAY_VIEW<int>& id=*spray_particles.array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID);
    for(int i=0;i<new_removed_particles.array_collection->Size();i++)if(new_removed_particles_spawn_time(i)<=time+dt&&new_removed_particles_spawn_time(i)>time  && !kill_particle_box.Lazy_Inside(new_removed_particles.X(i))){
        for(int k=0;k<new_removed_particles_amplification_factor(i);k++){
            int index=particle_to_add++;
            spray_particles.X(index)=new_removed_particles.X(i)+random.Get_Uniform_Vector(TV(-(T).1,-(T).1,-(T).1)*fluids_parameters.grid->dX.Max(),TV((T).1,(T).1,(T).1)*fluids_parameters.grid->dX.Max());
            spray_particles.V(index)=TV(new_removed_particles.V(i).x*random.Get_Uniform_Number((T).8,(T)1.2),new_removed_particles.V(i).y*random.Get_Uniform_Number((T).8,(T)1.2),new_removed_particles.V(i).z*random.Get_Uniform_Number((T).8,(T)1.2));
            spray_particles.radius(index)=random.Get_Uniform_Number((T).05,(T).9);
            spray_particles.age(index)=0;
            id(index)=unique_id++;}}
    LOG::cout<<"added "<<number_of_particles_to_add<<" out of "<<new_removed_particles.array_collection->Size()<<" new spray particles this frame"<<std::endl;

    LEVELSET_3D<GRID<TV> > interpolated_water_levelset(*fluids_parameters.grid,interpolated_water_phi);
    ARRAY<int> removed_particles_amplification_factor(last_frame_removed_particles.array_collection->Size());
    LOG::Time("spawning new particles from all negative removed particles");
    number_of_particles_to_add=0;
    /*
    //spawn dependent on phi
    T particle_phi_low=0;
    T particle_phi_high=2;
    int particle_number_spawn_low=0;
    int particle_number_spawn_high=1;
    for(int i=0;i<last_frame_removed_particles.array_collection->Size();i++)
        if(!kill_particle_box.Lazy_Inside(last_frame_removed_particles.X(i))){
            T particle_phi=interpolated_water_levelset.Phi(last_frame_removed_particles.X(i));
            particle_phi=clamp(particle_phi,particle_phi_low,particle_phi_high);
            T particle_number_spawn_scale=(particle_phi-particle_phi_low)/(particle_phi_high-particle_phi_low);
            int amplification_factor=particle_number_spawn_low+(int)(particle_number_spawn_scale*(particle_number_spawn_high-particle_number_spawn_low));
            amplification_factor*=random.Get_Uniform_Integer(0,amplification_factor*2);
            removed_particles_amplification_factor(i)=amplification_factor;
            number_of_particles_to_add+=amplification_factor;}
    */
    //spawn dependent on velocities
    T particle_number_spawn_low=0;
    T particle_number_spawn_high=(T).5;
    for(int i=0;i<last_frame_removed_particles.array_collection->Size();i++)
        if(!kill_particle_box.Lazy_Inside(last_frame_removed_particles.X(i))){
            T particle_velcocity_magnitude=last_frame_removed_particles.V(i).Magnitude();
            particle_velcocity_magnitude=clamp(particle_velcocity_magnitude,particle_velocity_low,particle_velocity_high);
            T particle_velocity_scale=(particle_velcocity_magnitude-particle_velocity_low)/(particle_velocity_high-particle_velocity_low);
            T amplification_factor=particle_number_spawn_low+particle_velocity_scale*(particle_number_spawn_high-particle_number_spawn_low);
            amplification_factor*=(dt*frame_rate);
            //amplification_factor*=random.Get_Uniform_Integer(0,amplification_factor*2);
            //removed_particles_amplification_factor(i)=amplification_factor;
            //number_of_particles_to_add+=amplification_factor;
            if(random.Get_Uniform_Number(0.0,1.0)<amplification_factor){
                removed_particles_amplification_factor(i)=1;
                number_of_particles_to_add++;}}

    particle_to_add=spray_particles.array_collection->Size()+1;
    spray_particles.array_collection->Add_Elements(number_of_particles_to_add);
    for(int i=0;i<last_frame_removed_particles.array_collection->Size();i++)if(!kill_particle_box.Lazy_Inside(last_frame_removed_particles.X(i))){
        for(int k=0;k<removed_particles_amplification_factor(i);k++){
            int index=particle_to_add++;
            spray_particles.X(index)=last_frame_removed_particles.X(i)+random.Get_Uniform_Vector(TV(-(T).1,-(T).1,-(T).1)*fluids_parameters.grid->dX.Max(),TV((T).1,(T).1,(T).1)*fluids_parameters.grid->dX.Max());
            spray_particles.V(index)=TV(last_frame_removed_particles.V(i).x*random.Get_Uniform_Number((T).8,(T)1.2),last_frame_removed_particles.V(i).y*random.Get_Uniform_Number((T).8,(T)1.2),last_frame_removed_particles.V(i).z*random.Get_Uniform_Number((T).8,(T)1.2));
            spray_particles.radius(index)=random.Get_Uniform_Number((T).05,(T).9);
            spray_particles.age(index)=0;
            id(index)=unique_id++;}}
    LOG::cout<<"added "<<number_of_particles_to_add<<" new spray particles"<<std::endl;

    LOG::Time("converting spray particles to foam");
    ARRAY<bool> convert(spray_particles.array_collection->Size());int number_to_convert=0;
    for(int i=spray_particles.array_collection->Size();i>=1;i--){
        if(interpolated_water_levelset.Phi(spray_particles.X(i))<-.2*fluids_parameters.grid->dX.Max()){
            if(random.Get_Uniform_Integer(1,8)==1){//to sample 1/8th times.
                convert(i)=true;number_to_convert++;}}}
    LOG::cout<<"marked "<<number_to_convert<<" particles for conversion to foam"<<std::endl;
    foam_particles.array_collection->Preallocate(foam_particles.array_collection->Size()+number_to_convert);
    for(int i=spray_particles.array_collection->Size();i>=1;i--)
    if(convert(i)) foam_particles.array_collection->Take(*spray_particles.array_collection,i);

    for(int i=spray_particles.array_collection->Size();i>=1;i--)
        if(interpolated_water_levelset.Phi(spray_particles.X(i))<-.2*fluids_parameters.grid->dX.Max())
            spray_particles.array_collection->Add_To_Deletion_List(i);
    spray_particles.array_collection->Delete_Elements_On_Deletion_List();

    LOG::Time(STRING_UTILITIES::string_sprintf("simulating %d spray particles",spray_particles.array_collection->Size()));
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
    LEVELSET_3D<GRID<TV> > object_levelset(*fluids_parameters.grid,object_phi);
    for(int i=0;i<spray_particles.array_collection->Size();i++){
        T b=(T)1-(T)std::pow((T)max(spray_particles.radius(i),(T)0),dt);
        spray_particles.V(i)=(b)*interpolation.Clamped_To_Array_Face(*fluids_parameters.grid,fluid_collection.incompressible_fluid_collection.face_velocities,spray_particles.X(i))+(1-b)*spray_particles.V(i);
        spray_particles.V(i).y-=dt*(T)9.8;
        T collision_distance=(T).5*fluids_parameters.grid->Minimum_Edge_Length();
        Adjust_Particle_For_Objects(spray_particles.X(i),spray_particles.V(i),collision_distance,PARTICLE_LEVELSET_REMOVED_NEGATIVE,dt,time);
        spray_particles.X(i)+=spray_particles.V(i)*dt;
        spray_particles.age(i)+=dt;}
    LOG::Time(STRING_UTILITIES::string_sprintf("simulating %d foam particles",foam_particles.array_collection->Size()));
    for(int i=0;i<foam_particles.array_collection->Size();i++){
        T b=(T)1-(T)std::pow((T)max(foam_particles.radius(i),(T)0),dt);// TODO: because radius doesn't change, should we change this?
        foam_particles.V(i)=(b)*interpolation.Clamped_To_Array_Face(*fluids_parameters.grid,interpolated_water_V,foam_particles.X(i))+(1-b)*foam_particles.V(i);
        b=(T).25;
        foam_particles.V(i)=(b)*interpolation.Clamped_To_Array_Face(*fluids_parameters.grid,fluid_collection.incompressible_fluid_collection.face_velocities,foam_particles.X(i))+(1-b)*foam_particles.V(i);
        T collision_distance=(T).5*fluids_parameters.grid->Minimum_Edge_Length();
        Adjust_Particle_For_Objects(foam_particles.X(i),foam_particles.V(i),collision_distance,PARTICLE_LEVELSET_REMOVED_NEGATIVE,dt,time);
        foam_particles.X(i)+=foam_particles.V(i)*dt;
        T phi=interpolated_water_levelset.Phi(foam_particles.X(i));
        int j=1;
        while(fabs(phi)>(T)1e-2){
            j++;
            TV normal=interpolated_water_levelset.Normal(foam_particles.X(i));
            foam_particles.X(i)-=normal*phi;
            phi=interpolated_water_levelset.Phi(foam_particles.X(i));
            if(j>5)break;}
        foam_particles.radius(i)-=dt*(T).33; // TODO: find a good rule for killing foam particles
        foam_particles.age(i)+=dt;
        if(foam_particles.radius(i)<0||foam_particles.age(i)>life_foam)
            foam_particles.array_collection->Add_To_Deletion_List(i);}
    LOG::Time("deleting foam particles on deletion list");
    foam_particles.array_collection->Delete_Elements_On_Deletion_List();
    LOG::cout<<"kept "<<foam_particles.array_collection->Size()<<" foam particles"<<std::endl;

    // remove bad particles
    LOG::Time("removing particles outside domain");

    LOG::Time("removing a small amount of density");
    for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();fluids_parameters.density_container.density(cell)=max((T)0,fluids_parameters.density_container.density(cell)-dt);}

    //LOG::Pop_Scope();
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(TV& X,TV& V,const T collision_distance,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) const // TODO: PHYSBAM_OVERRIDE Problem
{
    TV X_old=X,X_new=X_old+dt*V;
    int aggegate_id;T phi_old,phi;TV N,V_object;COLLISION_GEOMETRY_ID body_id;
    if(fluids_parameters.collision_bodies_affecting_fluid->Get_Body_Penetration(X_old,X_new,collision_distance,dt,body_id,aggegate_id,phi_old,phi,N,V_object)){
        X_new+=(max((T)0,-phi+collision_distance))*N;
        T relative_normal_velocity=TV::Dot_Product(V-V_object,N);if(relative_normal_velocity<0) V-=relative_normal_velocity*N;
        X=X_new-dt*V;
        if(fluids_parameters.collision_bodies_affecting_fluid->Any_Crossover(X_old,X_new,dt)) return false;}
    return true;
}
//#####################################################################
// Function Load_Simulated_Water
//#####################################################################
void Load_Simulated_Water(int frame)
{
    LOG::cout<<"Loading simulated water at frame: "<<frame<<std::endl;

    if(current_high_frame==frame){// the frame is already in the high frame... simply copy it
        LOG::cout<<"Copying water from high frame into low frame"<<std::endl;
        T_ARRAYS_SCALAR::Exchange_Arrays(low_frame_phi,high_frame_phi);
        T_FACE_ARRAYS_SCALAR::Exchange_Arrays(low_frame_V,high_frame_V);
        exchange(current_low_frame,current_high_frame);}
    else{
        LOG::cout<<"Loading low frame: "<<frame<<std::endl;
        LEVELSET_3D<GRID<TV> > low_frame_levelset(low_frame_grid,low_frame_phi);current_low_frame=frame;
        FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/levelset.%d",input_water_directory.c_str(),frame),low_frame_levelset);
        FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/mac_velocities.%d",input_water_directory.c_str(),frame),low_frame_V);}

    LOG::cout<<"Loading high frame: "<<frame+1<<std::endl;
    LEVELSET_3D<GRID<TV> > high_frame_levelset(high_frame_grid,high_frame_phi);current_high_frame=frame+1;
    FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/levelset.%d",input_water_directory.c_str(),(frame+1)),high_frame_levelset);
    FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/mac_velocities.%d",input_water_directory.c_str(),(frame+1)),high_frame_V);
}
//#####################################################################
// Function Compute_Interpolated_Water
//#####################################################################
void Compute_Interpolated_Water(T time)
{
    if(fabs(interpolated_water_time-time)<(T)1e-5)return;

    LEVELSET_3D<GRID<TV> > interpolated_water_levelset(*fluids_parameters.grid,interpolated_water_phi);
    LEVELSET_3D<GRID<TV> > object_levelset(*fluids_parameters.grid,object_phi);
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
    // convert time to frame
    T fractional_frame=time*frame_rate;
    // check if it's on a frame boundary
    if(fabs(floor(fractional_frame+(T).5)-fractional_frame)<(T)1e-5){
        LOG::cout<<"Interpolating single frame"<<std::endl;
        for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();TV X=fluids_parameters.grid->Center(cell);
            TV clamped_location=low_frame_grid.Clamp(X);
            T water_phi=interpolation.Clamped_To_Array_Cell(low_frame_grid,low_frame_phi,clamped_location);
            T object_phi=object_levelset.Extended_Phi(X);
            if(water_phi<object_phi) interpolated_water_phi(cell)=water_phi;
            else interpolated_water_phi(cell)=object_phi;};
        for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            TV_INT face=iterator.Face_Index();int axis=iterator.Axis();TV X=fluids_parameters.grid->Face(axis,face);
            TV clamped_location=low_frame_grid.Clamp(X);
            T water_phi=interpolation.Clamped_To_Array_Cell(low_frame_grid,low_frame_phi,clamped_location);
            T object_phi=object_levelset.Extended_Phi(X);
            if(water_phi<object_phi) interpolated_water_V(axis,face)=interpolation.Clamped_To_Array_Face_Component(axis,low_frame_grid,low_frame_V,clamped_location);
            else interpolated_water_V(axis,face)=0;}}
    else{
        int low_frame=(int)floor(fractional_frame);
        int high_frame=(int)ceil(fractional_frame);
        T low_weight=(high_frame-fractional_frame);
        T high_weight=(fractional_frame-low_frame);
        assert(fabs(1-(low_weight+high_weight))<(T)1e-5);

        LOG::cout<<"Interpolating multiple frame "<<fractional_frame<<" weights are:" <<low_weight<<" "<<high_weight<<std::endl;
        for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();TV X=fluids_parameters.grid->Center(cell);
            TV clamped_location=low_frame_grid.Clamp(X);
            T water_phi=low_weight*interpolation.Clamped_To_Array_Cell(low_frame_grid,low_frame_phi,clamped_location)+
                high_weight*interpolation.Clamped_To_Array_Cell(high_frame_grid,high_frame_phi,clamped_location);
            T object_phi=object_levelset.Extended_Phi(X);
            if(water_phi<object_phi) interpolated_water_phi(cell)=water_phi;
            else interpolated_water_phi(cell)=object_phi;}
        for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            TV_INT face=iterator.Face_Index();int axis=iterator.Axis();TV X=fluids_parameters.grid->Face(axis,face);
            TV clamped_location=low_frame_grid.Clamp(X);
            T water_phi=low_weight*interpolation.Clamped_To_Array_Cell(low_frame_grid,low_frame_phi,clamped_location)+
                high_weight*interpolation.Clamped_To_Array_Cell(high_frame_grid,high_frame_phi,clamped_location);
            T object_phi=object_levelset.Extended_Phi(X);
            if(water_phi<object_phi) interpolated_water_V(axis,face)=low_weight*interpolation.Clamped_To_Array_Face_Component(axis,low_frame_grid,low_frame_V,clamped_location)+
                high_weight*interpolation.Clamped_To_Array_Face_Component(axis,high_frame_grid,high_frame_V,clamped_location);
            else interpolated_water_V(axis,face)=0;}}
    interpolated_water_levelset.Fast_Marching_Method(time,3);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Find_Removed_Particles_That_Will_Appear_This_Frame(int particle_type,const int frame)
{
    new_removed_particles.array_collection->Clean_Memory();
    new_removed_particles_spawn_time.Clean_Memory();
    T frame_plus_one_time=(frame+1)/frame_rate;
    if(!use_log_removed_particle_times){
        // use the last frames particles to find out which particles are 'new'
        last_frame_removed_particles.array_collection->Delete_All_Elements();
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> > particles_array;
        if(particle_type==1) FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/removed_positive_particles.%d",input_water_directory.c_str(),frame),particles_array);
        else FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/removed_negative_particles.%d",input_water_directory.c_str(),frame),particles_array);
        for(int i=particles_array.domain.min_corner.x;i<=particles_array.domain.max_corner.x;i++)for(int j=particles_array.domain.min_corner.y;j<=particles_array.domain.max_corner.y;j++)for(int ij=particles_array.domain.min_corner.z;ij<=particles_array.domain.max_corner.z;ij++){
            if(particles_array(i,j,ij)){
                last_frame_removed_particles.array_collection->Append(*particles_array(i,j,ij)->array_collection);
                delete particles_array(i,j,ij);}}

        HASHTABLE<int,int> hashtable(last_frame_removed_particles.array_collection->Size());
        ARRAY_VIEW<int>& id=*last_frame_removed_particles.array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID);
        for(int i=0;i<last_frame_removed_particles.array_collection->Size();i++){
            hashtable.Insert(id(i),i);}

        // read in this frames particles
        if(particle_type==1) FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/removed_positive_particles.%d",input_water_directory.c_str(),(frame+1)),particles_array);
        else FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/removed_negative_particles.%d",input_water_directory.c_str(),(frame+1)),particles_array);
        for(int i=particles_array.domain.min_corner.x;i<=particles_array.domain.max_corner.x;i++)for(int j=particles_array.domain.min_corner.y;j<=particles_array.domain.max_corner.y;j++)for(int ij=particles_array.domain.min_corner.z;ij<=particles_array.domain.max_corner.z;ij++){
            if(particles_array(i,j,ij)){
                for(int k=0;k<particles_array(i,j,ij)->array_collection->Size();k++){
                    int id=(*particles_array(i,j,ij)->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID))(k);
                    if(!hashtable.Contains(id))
                      new_removed_particles.array_collection->Append(*particles_array(i,j,ij)->array_collection,k);}
                      delete particles_array(i,j,ij);}}
        LOG::cout<<"----There are "<<new_removed_particles.array_collection->Size()<<" new escaped particles in frame "<<frame<<std::endl;

        // take the particle at frame n+1, rewind using levelset at that time to the levelset position
        new_removed_particles_spawn_time.Resize(new_removed_particles.array_collection->Size());
        LEVELSET_3D<GRID<TV> > high_frame_levelset(high_frame_grid,high_frame_phi);
        LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;

        for(int i=0;i<new_removed_particles.array_collection->Size();i++){
            TV X=new_removed_particles.X(i);
            T phi=high_frame_levelset.Phi(X);
            TV particle_V=new_removed_particles.V(i);
            TV levelset_V=interpolation.Clamped_To_Array_Face(high_frame_grid,high_frame_V,X);
            TV relative_V=particle_V-levelset_V; // to compensate for the levelset moving
            T relative_V_magnitude=relative_V.Magnitude();
            if(relative_V_magnitude==0){
                new_removed_particles_spawn_time(i)=frame_plus_one_time;
                continue;}
            T time_from_levelset=clamp(phi/relative_V_magnitude,(T)0,(T)1/frame_rate);
            new_removed_particles_spawn_time(i)=frame_plus_one_time-time_from_levelset;
            new_removed_particles.X(i)=X-time_from_levelset*particle_V;}}
    else{
        ARRAY<PAIR<int,T> > removed_particle_times;
        FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/removed_particle_times.%d",input_water_directory.c_str(),frame),removed_particle_times);

        HASHTABLE<int,int> hashtable(removed_particle_times.m);
        for(int i=0;i<removed_particle_times.m;i++)
            hashtable.Insert(removed_particle_times(i).x,i);

        // read in this frames particles
        last_frame_removed_particles.array_collection->Delete_All_Elements();
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> > particles_array;
        if(particle_type==1) FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/removed_positive_particles.%d",input_water_directory.c_str(),(frame+1)),particles_array);
        else FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/removed_negative_particles.%d",input_water_directory.c_str(),(frame+1)),particles_array);
        for(int i=particles_array.domain.min_corner.x;i<=particles_array.domain.max_corner.x;i++)for(int j=particles_array.domain.min_corner.y;j<=particles_array.domain.max_corner.y;j++)for(int ij=particles_array.domain.min_corner.z;ij<=particles_array.domain.max_corner.z;ij++){
            if(particles_array(i,j,ij)){
                for(int k=0;k<particles_array(i,j,ij)->array_collection->Size();k++){
                    last_frame_removed_particles.array_collection->Append(*particles_array(i,j,ij)->array_collection,k);
                    int id=(*particles_array(i,j,ij)->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID))(k);
                    int entry_number;
                    if(hashtable.Get(id,entry_number)){
                        new_removed_particles_spawn_time.Append(removed_particle_times(entry_number).y);
                        new_removed_particles.array_collection->Append(*particles_array(i,j,ij)->array_collection,k);
                        new_removed_particles.X.Last()-=new_removed_particles.V.Last()*(frame_plus_one_time-removed_particle_times(entry_number).y);}}
                delete particles_array(i,j,ij);}}
        LOG::cout<<"----There are "<<new_removed_particles.array_collection->Size()<<" new escaped particles in frame "<<frame<<std::endl;}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/spray_particles.%d",output_directory.c_str(),frame),spray_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/foam_particles.%d",output_directory.c_str(),frame),foam_particles);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
void Read_Output_Files_Fluids(const int frame)
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::Read_Output_Files_Fluids(frame);
    Load_Simulated_Water(frame);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/spray_particles.%d",output_directory.c_str(),frame),spray_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/foam_particles.%d",output_directory.c_str(),frame),foam_particles);
    printf("read %d spray and %d foam particles\n\n",spray_particles.array_collection->Size(),foam_particles.array_collection->Size());
}
//#####################################################################
};
}
#endif
