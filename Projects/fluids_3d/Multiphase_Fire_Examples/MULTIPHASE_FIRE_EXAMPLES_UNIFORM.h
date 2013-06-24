//#####################################################################
// Copyright 2005, Frank Losasso
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIPHASE_FIRE_EXAMPLES_UNIFORM 
//#####################################################################
#ifndef __MULTIPHASE_FIRE_EXAMPLES_UNIFORM__
#define __MULTIPHASE_FIRE_EXAMPLES_UNIFORM__

#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform_Computations/SMOOTH_UNIFORM.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class MULTIPHASE_FIRE_EXAMPLES_UNIFORM:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::first_frame;using BASE::data_directory;using BASE::Adjust_Phi_With_Source;
    using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;
    using BASE::output_directory;using BASE::restart;using BASE::restart_frame;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;
    
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    int air_region;

    CYLINDER<T> inner_cylinder1,middle_cylinder1,outer_cylinder1;
    CYLINDER<T> inner_cylinder2,middle_cylinder2,outer_cylinder2;
    T source_shutoff_time;
    T sphere_drop_time;

    T reaction_bandwidth;
    bool pseudo_dirichlet;
    
    MULTIPHASE_FIRE_EXAMPLES_UNIFORM(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>(stream_type,0,fluids_parameters.FIRE),
        rigid_body_collection(solid_body_collection.rigid_body_collection),pseudo_dirichlet(false)
    {
    }

    virtual ~MULTIPHASE_FIRE_EXAMPLES_UNIFORM()
    {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add("-pseudo_dirichlet",&pseudo_dirichlet,"pseudo_dirichlet");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    fluids_parameters.Initialize_Number_Of_Regions(Number_Of_Regions(test_number));
    fluids_parameters.write_particles=true;
    PARAMETER_LIST parameters;
    fluids_parameters.use_reacting_flow=true;
    fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
    last_frame=512;frame_rate=96;
    fluids_parameters.temperature_container.Set_Ambient_Temperature(T(283.15));fluids_parameters.temperature_container.Set_Cooling_Constant((T)4000);
    fluids_parameters.density_container.Set_Ambient_Density(0);
    fluids_parameters.temperature_products=3000;fluids_parameters.temperature_fuel=298;
    fluids_parameters.temperature_buoyancy_constant=0;
    fluids_parameters.use_vorticity_confinement_fuel=fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.write_debug_data=fluids_parameters.write_velocity=true;

    if(test_number==1){
        fluids_parameters.incompressible_iterations=100;
        fluids_parameters.densities(1)=(T)1000;fluids_parameters.densities(2)=(T)1;
        fluids_parameters.normal_flame_speeds(1,2)=fluids_parameters.normal_flame_speeds(2,1)=(T).0003;
        fluids_parameters.fuel_region(1)=true;
        fluids_parameters.use_flame_speed_multiplier=true;
        fluids_parameters.confinement_parameters(2)=(T).2;
        air_region=2;}
    else if(test_number==2){
        fluids_parameters.incompressible_iterations=100;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.densities(1)=(T)1000;fluids_parameters.densities(2)=(T)800;fluids_parameters.densities(3)=(T)2000;fluids_parameters.densities(4)=(T)1;
        fluids_parameters.normal_flame_speeds(2,4)=fluids_parameters.normal_flame_speeds(4,2)=(T).0003;
        fluids_parameters.normal_flame_speeds(3,4)=fluids_parameters.normal_flame_speeds(4,3)=(T).0001;
        fluids_parameters.fuel_region(2)=true;fluids_parameters.fuel_region(3)=true;
        fluids_parameters.viscosities(3)=50;
        fluids_parameters.implicit_viscosity=true;
        fluids_parameters.use_flame_speed_multiplier=true;
        fluids_parameters.confinement_parameters(4)=(T).2;
        fluids_parameters.temperature_buoyancy_constant=(T)0.01;
        if(pseudo_dirichlet) fluids_parameters.pseudo_dirichlet_regions(4)=true;
        inner_cylinder1.Set_Endpoints(VECTOR<T,3>(0,(T).55,(T).3),VECTOR<T,3>((T).1,(T).55,(T).3));inner_cylinder1.radius=(T).06;
        middle_cylinder1.Set_Endpoints(VECTOR<T,3>(0,(T).55,(T).3),VECTOR<T,3>((T).1,(T).55,(T).3));middle_cylinder1.radius=(T).09;
        outer_cylinder1.Set_Endpoints(VECTOR<T,3>(0,(T).55,(T).3),VECTOR<T,3>((T).099,(T).55,(T).3));outer_cylinder1.radius=(T).1;
        source_shutoff_time=(T)2.5+(T)1e-5;
        sphere_drop_time=(T)2.9+(T)1e-5;
        air_region=4;}
    else if(test_number==3){
        fluids_parameters.implicit_viscosity_iterations=50;
        fluids_parameters.implicit_viscosity=true;
        frame_rate=48;
        fluids_parameters.densities(1)=(T)2000;
        fluids_parameters.viscosities(1)=(T)200;           
        fluids_parameters.densities(2)=(T)2000;
        fluids_parameters.viscosities(2)=(T)200;           
        fluids_parameters.densities(3)=(T)1000;
        fluids_parameters.densities(4)=(T)50;
        //fluids_parameters.surface_tensions(4,1)=fluids_parameters.surface_tensions(1,4)=(T)50;
        //fluids_parameters.surface_tensions(4,2)=fluids_parameters.surface_tensions(2,4)=(T)50;
        fluids_parameters.surface_tensions(4,3)=fluids_parameters.surface_tensions(3,4)=(T)10;
        fluids_parameters.normal_flame_speeds(1,2)=fluids_parameters.normal_flame_speeds(2,1)=(T).01;
        fluids_parameters.normal_flame_speeds(1,1)=fluids_parameters.normal_flame_speeds(2,2)=(T).01;
        fluids_parameters.fuel_region(1)=true;
        fluids_parameters.fuel_region(2)=true;

        fluids_parameters.use_flame_speed_multiplier=true;
        fluids_parameters.reseeding_frame_rate=5;
        //fluids_parameters.cfl/=2;

        reaction_bandwidth=2;

        inner_cylinder1.Set_Endpoints(VECTOR<T,3>(0,(T).2,(T).5),VECTOR<T,3>((T).1,(T).2,(T).5));inner_cylinder1.radius=(T).1;
        outer_cylinder1.Set_Endpoints(VECTOR<T,3>(0,(T).2,(T).5),VECTOR<T,3>((T).099,(T).2,(T).5));outer_cylinder1.radius=(T).11;
        inner_cylinder2.Set_Endpoints(VECTOR<T,3>((T).9,(T).2,(T).5),VECTOR<T,3>((T)1,(T).2,(T).5));inner_cylinder2.radius=(T).1;
        outer_cylinder2.Set_Endpoints(VECTOR<T,3>((T).899,(T).2,(T).5),VECTOR<T,3>((T)1,(T).2,(T).5));outer_cylinder2.radius=(T).11;}
    else if(test_number==4){
        fluids_parameters.incompressible_iterations=100;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.densities(1)=(T)1000;fluids_parameters.densities(2)=(T)800;fluids_parameters.densities(3)=(T)2000;fluids_parameters.densities(4)=(T)1;
        fluids_parameters.normal_flame_speeds(2,4)=fluids_parameters.normal_flame_speeds(4,2)=(T).0003;
        fluids_parameters.normal_flame_speeds(3,4)=fluids_parameters.normal_flame_speeds(4,3)=(T).0001;
        fluids_parameters.fuel_region(2)=true;fluids_parameters.fuel_region(3)=true;
        fluids_parameters.viscosities(3)=50;
        fluids_parameters.implicit_viscosity=true;
        fluids_parameters.use_flame_speed_multiplier=true;
        fluids_parameters.confinement_parameters(4)=(T).2;
        fluids_parameters.temperature_buoyancy_constant=(T)0.01;
        if(pseudo_dirichlet) fluids_parameters.pseudo_dirichlet_regions(4)=true;
        inner_cylinder1.Set_Endpoints(VECTOR<T,3>(0,(T).55,(T).3),VECTOR<T,3>((T).1,(T).55,(T).3));inner_cylinder1.radius=(T).06;
        middle_cylinder1.Set_Endpoints(VECTOR<T,3>(0,(T).55,(T).3),VECTOR<T,3>((T).1,(T).55,(T).3));middle_cylinder1.radius=(T).09;
        outer_cylinder1.Set_Endpoints(VECTOR<T,3>(0,(T).55,(T).3),VECTOR<T,3>((T).099,(T).55,(T).3));outer_cylinder1.radius=(T).1;
        source_shutoff_time=(T)2.5+(T)1e-5;
        sphere_drop_time=(T)2.9+(T)1e-5;
        air_region=4;}
}
//#####################################################################
// Function Number_Of_Regions
//#####################################################################
static int Number_Of_Regions(int test_number)
{
    if(test_number==1) return 2;
    if(test_number==2) return 4;
    if(test_number==3) return 4;
    LOG::cout<<"Unrecognized example: "<<test_number<<std::endl;
    assert(false);exit(0);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Get_Flame_Speed_Multiplier
//#####################################################################
void Get_Flame_Speed_Multiplier(const T dt,const T time) PHYSBAM_OVERRIDE
{
    T_FACE_ARRAYS_SCALAR& flame_speed_multiplier=fluids_parameters.incompressible->projection.flame_speed_multiplier;
    flame_speed_multiplier.Fill(0);

    if(test_number==1){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
            if((T).45<iterator.Location().x && iterator.Location().x<=(T).55  &&  (T_GRID::dimension==2 || ((T).45<iterator.Location().z && iterator.Location().z<=(T).55)))
                flame_speed_multiplier.Component(iterator.Axis())(iterator.Face_Index())=1;}
    else if(test_number==2){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
            if(time<sphere_drop_time && SPHERE<TV>(VECTOR<T,3>((T).25,(T).6,(T).75),(T).08).Lazy_Inside(iterator.Location()))
                flame_speed_multiplier.Component(iterator.Axis())(iterator.Face_Index())=1;}
    else if(test_number==3){
        ARRAY<ARRAY<T,VECTOR<int,3> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
        T reaction_bandwidth_times_edge_length=reaction_bandwidth*fluids_parameters.grid->dX.Min();

        T_ARRAYS_SCALAR& phi1=phis(1);
        //LEVELSET<TV> levelset1(fluids_parameters.grid,phi1);
        //levelset1.Set_Band_Width(2*reaction_bandwidth_times_edge_length+fluids_parameters.grid.dX.Min());levelset1.Fast_Marching_Method();
        T_ARRAYS_SCALAR& phi2=phis(2);
        //LEVELSET<TV> levelset2(fluids_parameters.grid,phi2);
        //levelset2.Set_Band_Width(2*reaction_bandwidth_times_edge_length+fluids_parameters.grid.dX.Min());levelset2.Fast_Marching_Method();
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            VECTOR<int,3> cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
            T face_phi1=(T).5*(phi1(cell1)+phi1(cell2)),face_phi2=(T).5*(phi2(cell1)+phi2(cell2));
            if(face_phi1<reaction_bandwidth_times_edge_length && face_phi2<reaction_bandwidth_times_edge_length){
                T min_phi=min(face_phi1,face_phi2);
                flame_speed_multiplier.Component(iterator.Axis())(iterator.Face_Index())=clamp(reaction_bandwidth_times_edge_length+min_phi,(T)0,(T)1);}}
        return;}

    T ignition_temperature=1075; // 1100 is .5 as fast as 1000, 1150 is .25 as fast as 1000.
    T_ARRAYS_SCALAR temp_temperature(fluids_parameters.grid->Domain_Indices(3));
    BOUNDARY<TV,T> boundary;boundary.Fill_Ghost_Cells(*fluids_parameters.grid,fluids_parameters.temperature_container.temperature,temp_temperature,dt,time,3);
    SMOOTH::Smooth<T,TV::m>(temp_temperature,5,0);
    if(fluids_parameters.mpi_grid) fluids_parameters.mpi_grid->Exchange_Boundary_Cell_Data(temp_temperature,1);
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        T face_temperature=(T).5*(temp_temperature(iterator.First_Cell_Index())+temp_temperature(iterator.Second_Cell_Index()));
        if(face_temperature>ignition_temperature) flame_speed_multiplier.Component(iterator.Axis())(iterator.Face_Index())=clamp((T).01*(face_temperature-ignition_temperature),(T)0,(T)1);}
}
//#####################################################################
// Function Set_Ghost_Density_And_Temperature_Inside_Flame_Core
//#####################################################################
void Set_Ghost_Density_And_Temperature_Inside_Flame_Core() PHYSBAM_OVERRIDE
{
    if(test_number==3) return;

    T_ARRAYS_SCALAR phi;LEVELSET<TV> levelset(*fluids_parameters.grid,phi);
    T_FACE_ARRAYS_SCALAR& flame_speed_multiplier=fluids_parameters.incompressible->projection.flame_speed_multiplier;
    TEMPERATURE_CONTAINER<T_GRID>& temperature=fluids_parameters.temperature_container;
    DENSITY_CONTAINER<T_GRID>& density=fluids_parameters.density_container;

    fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(fluids_parameters.fuel_region,levelset,false);

    T bandwidth=2*fluids_parameters.grid->dX.Min();
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
        T cell_flame_speed_multiplier=0;
        for(int i=0;i<T_GRID::dimension;i++)
            cell_flame_speed_multiplier+=flame_speed_multiplier.Component(i)(iterator.First_Face_Index(i))+flame_speed_multiplier.Component(i)(iterator.Second_Face_Index(i));
        cell_flame_speed_multiplier/=2*T_GRID::dimension;
        if(-phi(index)<0){
            if(-bandwidth<-phi(index) && fluids_parameters.particle_levelset_evolution_multiple->Levelset(air_region).phi(index)<bandwidth){
                temperature.temperature(index)=LINEAR_INTERPOLATION<T,T>::Linear(temperature.ambient_temperature,fluids_parameters.temperature_products,cell_flame_speed_multiplier);
                density.density(index)=LINEAR_INTERPOLATION<T,T>::Linear(density.ambient_density,fluids_parameters.density,cell_flame_speed_multiplier);}
            else temperature.temperature(index)=temperature.ambient_temperature;}}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    ARRAY<T_ARRAYS_SCALAR>& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;

    if(test_number==1){
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
            phis(1)(iterator.Cell_Index())=iterator.Location().y-(T).2;
        phis(2).Copy(-1,phis(1));}
    if(test_number==2){
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            TV X=iterator.Location();TV_INT cell=iterator.Cell_Index();
            phis(1)(cell)=X.y-(T).2;
            phis(2)(cell)=1;
            phis(3)(cell)=SPHERE<TV>(VECTOR<T,3>((T).25,(T).6,(T).75),(T).08).Signed_Distance(X);
            phis(4)(cell)=-min(phis(1)(cell),phis(2)(cell),phis(3)(cell));}}
    if(test_number==3){
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            VECTOR<T,3> X=iterator.Location();VECTOR<int,3> cell=iterator.Cell_Index();
            phis(1)(cell)=1;
            phis(2)(cell)=1;
            phis(3)(cell)=X.y-(T).35;
            phis(4)(cell)=-min(phis(1)(cell),phis(2)(cell),phis(3)(cell));}}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    if(test_number==2){
        if(time<source_shutoff_time) Adjust_Phi_With_Source(inner_cylinder1,2,MATRIX<T,4>::Identity_Matrix());
        else if(time<source_shutoff_time+((T)1./frame_rate)) Adjust_Phi_With_Source(inner_cylinder1,4,MATRIX<T,4>::Identity_Matrix());
        if(time<sphere_drop_time) Adjust_Phi_With_Source(SPHERE<TV>(VECTOR<T,3>((T).25,(T).6,(T).75),(T).08),3,MATRIX<T,4>::Identity_Matrix());}
    if(test_number==3){
        Adjust_Phi_With_Source(inner_cylinder1,1,MATRIX<T,4>::Identity_Matrix());
        Adjust_Phi_With_Source(inner_cylinder2,2,MATRIX<T,4>::Identity_Matrix());
        ARRAY<ARRAY<T,VECTOR<int,3> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
        T reaction_seed_bandwidth=(T).5*reaction_bandwidth*fluids_parameters.grid->dX.Min();
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            VECTOR<int,3> cell=iterator.Cell_Index();
            T band=max(fabs(phis(1)(cell)),fabs(phis(2)(cell)))-reaction_seed_bandwidth;
            int region1,region2;T min_phi1,min_phi2;
            fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Two_Minimum_Regions(cell,region1,region2,min_phi1,min_phi2);
            if(min_phi1<-3*fluids_parameters.grid->dX.Min()) continue;
            //if(!(region1==1 && region2==2 || region1==2 && region2==1)) continue;
            phis(1)(cell)=max(phis(1)(cell),-band);
            phis(2)(cell)=max(phis(2)(cell),-band);
            phis(3)(cell)=max(phis(3)(cell),-band);
            phis(4)(cell)=min(phis(4)(cell),band);}}
    return false;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
{
    bool first=true;
    if(test_number==2 && time<source_shutoff_time){
        Get_Source_Reseed_Mask(inner_cylinder1,MATRIX<T,4>::Identity_Matrix(),cell_centered_mask,first);first=false;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
   if(test_number==2){
       Get_Source_Velocities(outer_cylinder1,MATRIX<T,4>::Identity_Matrix(),VECTOR<T,3>());
       if(time<source_shutoff_time)
           Get_Source_Velocities(middle_cylinder1,MATRIX<T,4>::Identity_Matrix(),VECTOR<T,3>(2.5,0,0));
       if(time<sphere_drop_time)
           Get_Source_Velocities(SPHERE<TV>(VECTOR<T,3>((T).25,(T).6,(T).75),(T).08),MATRIX<T,4>::Identity_Matrix(),VECTOR<T,3>());}
    if(test_number==3){
        Get_Source_Velocities(outer_cylinder1,MATRIX<T,4>::Identity_Matrix(),VECTOR<T,3>());
        Get_Source_Velocities(outer_cylinder2,MATRIX<T,4>::Identity_Matrix(),VECTOR<T,3>());
        Get_Source_Velocities(inner_cylinder1,MATRIX<T,4>::Identity_Matrix(),VECTOR<T,3>((T).2,0,0));
        Get_Source_Velocities(inner_cylinder2,MATRIX<T,4>::Identity_Matrix(),VECTOR<T,3>(-(T).2,0,0));}
}
//#####################################################################
};    
}
#endif
