//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIPHASE_FIRE_EXAMPLES 
//#####################################################################
#ifndef __MULTIPHASE_FIRE_EXAMPLES__
#define __MULTIPHASE_FIRE_EXAMPLES__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Grid_Tools/Computations/SMOOTH_UNIFORM.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class MULTIPHASE_FIRE_EXAMPLES:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::first_frame;using BASE::data_directory;
    using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;using BASE::Get_Source_Velocities;using BASE::resolution;
    using BASE::output_directory;using BASE::restart;using BASE::restart_frame;using BASE::solid_body_collection;using BASE::test_number;
    using BASE::Adjust_Phi_With_Source;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    bool use_object;

    T reaction_bandwidth;
    
    MULTIPHASE_FIRE_EXAMPLES(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,0,fluids_parameters.FIRE),
         rigid_body_collection(solid_body_collection.rigid_body_collection)
    {
        parse_args.Parse();
        fluids_parameters.Initialize_Number_Of_Regions(Number_Of_Regions(test_number));
        fluids_parameters.write_particles=true;
        PARAMETER_LIST parameters;
        fluids_parameters.use_reacting_flow=true;
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,10*resolution+1),RANGE<TV>(TV(0,0),TV(1,1)));
        output_directory=LOG::sprintf("Multiphase_Fire_Examples/Example_%d__Resolution_%d_%d",test_number,(fluids_parameters.grid->counts.x-1),(fluids_parameters.grid->counts.x-1));
        fluids_parameters.number_particles_per_cell=16;
        last_frame=512;
        fluids_parameters.temperature_container.Set_Ambient_Temperature(T(283.15));fluids_parameters.temperature_container.Set_Cooling_Constant((T)4000);
        fluids_parameters.density_container.Set_Ambient_Density(0);
        fluids_parameters.temperature_products=3000;fluids_parameters.temperature_fuel=298;
        fluids_parameters.temperature_buoyancy_constant=0;
        fluids_parameters.use_vorticity_confinement_fuel=fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.confinement_parameters(2)=(T).2;
        fluids_parameters.write_debug_data=fluids_parameters.write_velocity=true;
        fluids_parameters.object_friction=use_object?(T)1:(T)0;
        fluids_parameters.incompressible_iterations=50;

        if(test_number==1){
            fluids_parameters.densities(1)=(T)1000;
            fluids_parameters.densities(2)=(T)1;
            fluids_parameters.normal_flame_speeds(1,2)=fluids_parameters.normal_flame_speeds(2,1)=(T).0003;
            fluids_parameters.normal_flame_speeds(1,1)=fluids_parameters.normal_flame_speeds(2,2)=(T)0;
            fluids_parameters.fuel_region(1)=true;
            fluids_parameters.fuel_region(2)=false;
            fluids_parameters.use_flame_speed_multiplier=true;}
        if(test_number==2){
            fluids_parameters.implicit_viscosity_iterations=50;
            fluids_parameters.implicit_viscosity=true;

            fluids_parameters.densities(1)=(T)2000;
            fluids_parameters.viscosities(1)=(T)200;           
            fluids_parameters.densities(2)=(T)2000;
            fluids_parameters.viscosities(2)=(T)200;           
            fluids_parameters.densities(3)=(T)1000;          
            fluids_parameters.densities(4)=(T)50;
            //fluids_parameters.surface_tensions(4,1)=fluids_parameters.surface_tensions(1,4)=(T)50;
            //fluids_parameters.surface_tensions(4,2)=fluids_parameters.surface_tensions(2,4)=(T)50;
            fluids_parameters.surface_tensions(4,3)=fluids_parameters.surface_tensions(3,4)=(T)10;
            fluids_parameters.normal_flame_speeds(1,4)=fluids_parameters.normal_flame_speeds(4,1)=(T).01;
            fluids_parameters.normal_flame_speeds(4,2)=fluids_parameters.normal_flame_speeds(2,4)=(T).01;
            fluids_parameters.fuel_region(1)=true;
            fluids_parameters.fuel_region(2)=true;
            fluids_parameters.use_flame_speed_multiplier=true;
            fluids_parameters.reseeding_frame_rate=5;
            //fluids_parameters.cfl/=4;
            
            reaction_bandwidth=2;
        }
    }

    virtual ~MULTIPHASE_FIRE_EXAMPLES()
    {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Number_Of_Regions
//#####################################################################
static int Number_Of_Regions(int test_number)
{
    if(test_number==1) return 2;
    if(test_number==2) return 4;
    LOG::cout<<"Unrecognized example: "<<test_number<<std::endl;
    assert(false);exit(0);
}
//#####################################################################
// Function Get_Flame_Speed_Multiplier
//#####################################################################
void Get_Flame_Speed_Multiplier(const T dt,const T time) override
{
    ARRAY<T,FACE_INDEX<TV::m> >& flame_speed_multiplier=fluids_parameters.incompressible->projection.flame_speed_multiplier;
    flame_speed_multiplier.Fill(0);

    if(test_number==1){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
            if(.45<iterator.Location().x&&iterator.Location().x<=.55) flame_speed_multiplier.Component(iterator.Axis())(iterator.Face_Index())=1;

        T ignition_temperature=1150; // 1100 is .5 as fast as 1000, 1150 is .25 as fast as 1000.
        ARRAY<T,TV_INT> temp_temperature(fluids_parameters.grid->Domain_Indices(3));
        BOUNDARY<TV,T> boundary;boundary.Fill_Ghost_Cells(*fluids_parameters.grid,fluids_parameters.temperature_container.temperature,temp_temperature,dt,time,3);
        SMOOTH::Smooth<T,TV::m>(temp_temperature,5,0);
        if(fluids_parameters.mpi_grid)fluids_parameters.mpi_grid->Exchange_Boundary_Cell_Data(temp_temperature,1);
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            T face_temperature=(T).5*(temp_temperature(iterator.First_Cell_Index())+temp_temperature(iterator.Second_Cell_Index()));
            if(face_temperature>ignition_temperature) flame_speed_multiplier.Component(iterator.Axis())(iterator.Face_Index())=clamp((T).01*(face_temperature-ignition_temperature),(T)0,(T)1);}}
    if(test_number==2){
        ARRAY<ARRAY<T,VECTOR<int,2> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
        T reaction_bandwidth_times_edge_length=reaction_bandwidth*fluids_parameters.grid->dX.Min();

        ARRAY<T,TV_INT>& phi1=phis(1);
        //T_FAST_LEVELSET levelset1(fluids_parameters.grid,phi1);
        //levelset1.Set_Band_Width(2*reaction_bandwidth_times_edge_length+fluids_parameters.grid.dX.Min());levelset1.Fast_Marching_Method();
        ARRAY<T,TV_INT>& phi2=phis(2);
        //T_FAST_LEVELSET levelset2(fluids_parameters.grid,phi2);
        //levelset2.Set_Band_Width(2*reaction_bandwidth_times_edge_length+fluids_parameters.grid.dX.Min());levelset2.Fast_Marching_Method();
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            VECTOR<int,2> cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
            T face_phi1=(T).5*(phi1(cell1)+phi1(cell2)),face_phi2=(T).5*(phi2(cell1)+phi2(cell2));
            if(face_phi1<reaction_bandwidth_times_edge_length&&face_phi2<reaction_bandwidth_times_edge_length){
                T min_phi=min(face_phi1,face_phi2);
                flame_speed_multiplier.Component(iterator.Axis())(iterator.Face_Index())=clamp(reaction_bandwidth_times_edge_length+min_phi,(T)0,(T)1);}}}
/*
        for(FACE_ITERATOR<TV> iterator(fluids_parameters.grid);iterator.Valid();iterator.Next()){
            VECTOR<int,2> cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
            T face_phi1=(T).5*(phis(1)(cell1)+phis(1)(cell2)),face_phi2=(T).5*(phis(2)(cell1)+phis(2)(cell2));
            if(face_phi1<reaction_bandwidth_times_edge_length&&face_phi2<reaction_bandwidth_times_edge_length){
                T min_phi=min(face_phi1,face_phi2);
                flame_speed_multiplier.Component(iterator.Axis())(iterator.Face_Index())=clamp(reaction_bandwidth_times_edge_length+min_phi,(T)0,(T)1);}}}*/
}
//#####################################################################
// Function Set_Ghost_Density_And_Temperature_Inside_Flame_Core
//#####################################################################
void Set_Ghost_Density_And_Temperature_Inside_Flame_Core() override
{
    ARRAY<T,TV_INT> phi;LEVELSET<TV> levelset(*fluids_parameters.grid,phi);
    ARRAY<T,FACE_INDEX<TV::m> >& flame_speed_multiplier=fluids_parameters.incompressible->projection.flame_speed_multiplier;

    if(test_number==1){
        fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(fluids_parameters.fuel_region,levelset,false);
        
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            T cell_flame_speed_multiplier=0;
            for(int i=0;i<2;i++)
                cell_flame_speed_multiplier+=flame_speed_multiplier.Component(i)(iterator.First_Face_Index(i))+flame_speed_multiplier.Component(i)(iterator.Second_Face_Index(i));
            cell_flame_speed_multiplier*=.25;
            if(-phi(iterator.Cell_Index())<0){
                if(-2*fluids_parameters.grid->dX.Min()<-phi(iterator.Cell_Index()))
                    fluids_parameters.temperature_container.temperature(iterator.Cell_Index())=
                        cell_flame_speed_multiplier*fluids_parameters.temperature_container.hot_point+(1-cell_flame_speed_multiplier)*fluids_parameters.temperature_container.ambient_temperature;
                else fluids_parameters.temperature_container.temperature(iterator.Cell_Index())=fluids_parameters.temperature_container.ambient_temperature;}}}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    ARRAY<ARRAY<T,VECTOR<int,2> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    if(test_number==1){
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
            phis(1)(iterator.Cell_Index())=min(iterator.Location().y-(T).2,RANGE<TV>(TV(),TV((T).2,(T).5)).Signed_Distance(iterator.Location()));
        //phis(1)(iterator.Cell_Index())=iterator.Location().y-(T).2;
        phis(2).Copy(-1,phis(1));}
    if(test_number==2){
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            VECTOR<T,2> X=iterator.Location();VECTOR<int,2> cell=iterator.Cell_Index();
            phis(1)(cell)=1;//RANGE<TV>(0,1,.0,.2).Signed_Distance(X);
            phis(2)(cell)=1;//RANGE<TV>(.45,.55,.4,.5).Signed_Distance(X);
            phis(3)(cell)=-min(-X.y+(T).7,phis(1)(cell),phis(2)(cell));
            phis(4)(cell)=-min(phis(1)(cell),phis(2)(cell),phis(3)(cell));}}
//            phis(1)(cell)=fabs(X.y-.2)-.2;
//            phis(2)(cell)=fabs(X.y-.6)-.2;
//            phis(3)(cell)=1;
//            phis(4)(cell)=-min(phis(1)(cell),phis(2)(cell),phis(3)(cell));}}
}

//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) override
{
    if(test_number==2){
        Adjust_Phi_With_Source(RANGE<TV>(TV(0,(T).1),TV((T).1,(T).3)),1,MATRIX<T,3>::Identity_Matrix());
        Adjust_Phi_With_Source(RANGE<TV>(TV((T).9,(T).1),TV(1,(T).3)),2,MATRIX<T,3>::Identity_Matrix());
        ARRAY<ARRAY<T,VECTOR<int,2> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
        T reaction_seed_bandwidth=(T).5*reaction_bandwidth*fluids_parameters.grid->dX.Min();
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            VECTOR<int,2> cell=iterator.Cell_Index();
            T band=max(fabs(phis(1)(cell)),fabs(phis(2)(cell)))-reaction_seed_bandwidth;
            int region1,region2;T min_phi1,min_phi2;
            fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Two_Minimum_Regions(cell,region1,region2,min_phi1,min_phi2);
            if(min_phi1<-3*fluids_parameters.grid->dX.Min())continue;
            //if(!(region1==1&&region2==2||region1==2&&region2==1))continue;
            phis(1)(cell)=max(phis(1)(cell),-band);
            phis(2)(cell)=max(phis(2)(cell),-band);
            phis(3)(cell)=max(phis(3)(cell),-band);
            phis(4)(cell)=min(phis(4)(cell),band);}}
    return false;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) override
{
    if(test_number==2){
        Get_Source_Velocities(RANGE<TV>(TV(0,(T).1),TV((T).1,(T).3)),MATRIX<T,3>::Identity_Matrix(),VECTOR<T,2>((T).2,0));
        Get_Source_Velocities(RANGE<TV>(TV((T).9,(T).1),TV(1,(T).3)),MATRIX<T,3>::Identity_Matrix(),VECTOR<T,2>(-(T).2,0));}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    override
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
};    
}
#endif
