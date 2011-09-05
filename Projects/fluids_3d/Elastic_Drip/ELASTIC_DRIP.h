//#####################################################################
// Copyright 2004-2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTIC_DRIP 
//##################################################################### 
#ifndef __ELASTIC_DRIP__
#define __ELASTIC_DRIP__

#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class ELASTIC_DRIP:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;

    RIGID_BODY_LIST_3D<T> rigid_body_list;
    ARRAY<T,VECTOR<int,3> > phi_object;
    LEVELSET_3D<GRID<TV> > levelset_object;
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > V_object;
    CYLINDER<T> source;
    VECTOR<T,3> source_velocity;
    bool moving_objects;
    int example;

    ELASTIC_DRIP()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.WATER),levelset_object(fluids_parameters.grid,phi_object),
        source(VECTOR<T,3>(.5,.9,.5),VECTOR<T,3>(.5,1,.5),.1),source_velocity(0,-.8,0),moving_objects(false),example(0)
    {
        first_frame=0;last_frame=300;
        frame_rate=24;
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.second_order_cut_cell_method=true;

        restart=false;restart_frame=0;
        fluids_parameters.use_strain=true;fluids_parameters.write_strain=true;
        fluids_parameters.scalar_substeps=1;
        fluids_parameters.object_friction=1;
        fluids_parameters.adhesion_coefficient=0;

        example=4;

        if(example==1){
            output_directory="Elastic_Drip/output1";
            //fluids_parameters.grid.Initialize(60,100,60,.4,1.6,0,2,.4,1.6);
            fluids_parameters.grid.Initialize(30,50,30,.4,1.6,0,2,.4,1.6);
            source.Set_Endpoints(VECTOR<T,3>(1,1.9,1),VECTOR<T,3>(1,2,1));
            source.radius=.1;
            source_velocity=VECTOR<T,3>(0,-.8,0);
            fluids_parameters.viscosity=(T)200;
            fluids_parameters.elastic_modulus=4000;
            fluids_parameters.plasticity_alpha=1;
            fluids_parameters.plasticity_gamma=0;}
        else if(example==2){
            output_directory="Elastic_Drip/output2";
            //fluids_parameters.grid.Initialize(60,100,60,.4,1.6,0,2,.4,1.6);
            fluids_parameters.grid.Initialize(30,50,30,.4,1.6,0,2,.4,1.6);
            source.Set_Endpoints(VECTOR<T,3>(1,1.9,1),VECTOR<T,3>(1,2,1));
            source.radius=.1;
            source_velocity=VECTOR<T,3>(0,-1.6,0);
            fluids_parameters.viscosity=(T)200;
            fluids_parameters.elastic_modulus=4000;
            fluids_parameters.plasticity_alpha=1;
            fluids_parameters.plasticity_gamma=0;}
        else if(example==3){
            output_directory="Elastic_Drip/sphere_output";
            //fluids_parameters.grid.Initialize(60,100,60,.4,1.6,0,2,.4,1.6);
            fluids_parameters.grid.Initialize(30,50,30,.4,1.6,0,2,.4,1.6);
            source.Set_Endpoints(VECTOR<T,3>(1,1.9,1),VECTOR<T,3>(1,2,1));
            source.radius=.1;
            source_velocity=VECTOR<T,3>(0,-1.6,0);
            fluids_parameters.viscosity=(T)200;
            fluids_parameters.elastic_modulus=4000;
            fluids_parameters.plasticity_alpha=1;
            fluids_parameters.plasticity_gamma=0;
            rigid_body_list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere",(T).15,true,true,false);
            rigid_body_particles.Rigid_Body(1).position=VECTOR<T,3>(1,.7,1);}
        else if(example==4){
            output_directory="Elastic_Drip/two_sphere_output";
            //fluids_parameters.grid.Initialize(60,100,60,.4,1.6,0,2,.4,1.6);
            fluids_parameters.grid.Initialize(30,50,30,.4,1.6,0,2,.4,1.6);
            source.Set_Endpoints(VECTOR<T,3>(1,1.9,1),VECTOR<T,3>(1,2,1));
            //source.radius=.1;
            source.radius=.2;
            source_velocity=VECTOR<T,3>(0,-1.6,0);
            fluids_parameters.viscosity=(T)200;
            fluids_parameters.elastic_modulus=4000;
            fluids_parameters.plasticity_alpha=1;
            fluids_parameters.plasticity_gamma=0;
            moving_objects=true;
            rigid_body_list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere",(T).15,true,true,false);
            rigid_body_list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere",(T).15,true,true,false);}
        else{std::cerr<<"Dying inexplicably.\n";exit(1);}
    }

    ~ELASTIC_DRIP()
    {}

///####################################################################
// Function Construct_Levelset_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    if(example==4){
        rigid_body_particles.Rigid_Body(1).position=VECTOR<T,3>(.9,.7,1);
        rigid_body_particles.Rigid_Body(2).position=VECTOR<T,3>(1.1,.7,1);
        rigid_body_particles.Rigid_Body(1).velocity=VECTOR<T,3>(0,-.1,0);
        rigid_body_particles.Rigid_Body(2).velocity=VECTOR<T,3>(0,.1,0);
        for(int r=1;r<=rigid_body_list.rigid_bodies.m;r++)rigid_body_particles.Rigid_Body(r).position+=time*rigid_body_particles.Rigid_Body(r).velocity;}

    static bool initialized=false;
    GRID<TV>& grid=fluids_parameters.grid;
    if(!rigid_body_list.rigid_bodies.m || initialized) return;
    std::cout<<"Constructing levelsets"<<std::endl;
    phi_object.Resize(grid,3,false,false);
    V_object.Resize(grid,3,false,false);
    for(int i=-2;i<=grid.m+3;i++)for(int j=-2;j<=grid.n+3;j++)for(int ij=-2;ij<=grid.mn+3;ij++){
        T min_phi=1;int min_r=1;VECTOR<T,3> X=grid.X(i,j,ij);
        for(int r=1;r<=rigid_body_list.rigid_bodies.m;r++){
            T phi=rigid_body_particles.Rigid_Body(r).Implicit_Surface_Value(X);
            if(min_phi>phi){min_phi=phi;min_r=r;}}
        phi_object(i,j,ij)=min_phi;
        V_object(i,j,ij)=rigid_body_particles.Rigid_Body(min_r).Pointwise_Object_Velocity(X);}
    levelset_object.Compute_Cell_Minimum_And_Maximum(true);
    std::cout<<"done"<<std::endl;
    if(!moving_objects) initialized=true;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    ARRAY<T,VECTOR<int,3> >::copy(1,fluids_parameters.particle_levelset_evolution.phi);
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE
{
    if(!rigid_body_list.rigid_bodies.m) return;
    fluids_parameters.Extrapolate_Phi_Into_Object(phi_object);
}
//#####################################################################
// Function Adjust_Phi_with_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    if(!rigid_body_list.rigid_bodies.m) return;
    fluids_parameters.Adjust_Phi_With_Object(phi_object,V_object,time);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution.phi;
    for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)for(int ij=1;ij<=grid.mn;ij++){
        VECTOR<T,3> X=grid.X(i,j,ij);
        if(source.Lazy_Inside(X))phi(i,j,ij)=min(phi(i,j,ij),source.Signed_Distance(X));}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,3> >(grid);
    for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)for(int ij=1;ij<=grid.mn;ij++)if(source.Lazy_Inside(grid.X(i,j,ij)))(*cell_centered_mask)(i,j,ij)=true;
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR<T,3> >& particles,const int index,VECTOR<T,3>& V,
    const typename PARTICLE_LEVELSET<T,VECTOR<T,3> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(!rigid_body_list.rigid_bodies.m) return true;
    fluids_parameters.Adjust_Particle_For_Object(levelset_object,V_object,particles,index,V,particle_type,dt,time);
    return true;
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR<T,3> >& particles,const typename PARTICLE_LEVELSET<T,VECTOR<T,3> >::PARTICLE_TYPE particle_type,const T time)
{
    T contour_value=particle_type == PARTICLE_LEVELSET<T,VECTOR<T,3> >::NEGATIVE || particle_type == PARTICLE_LEVELSET<T,VECTOR<T,3> >::REMOVED_NEGATIVE ? -fluids_parameters.grid.dx : 0;
    for(int k=particles.array_collection->Size();k>=1;k--)if(levelset_object.Lazy_Inside(particles.X(k),contour_value)) particles.Delete_Particle(k);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV> &grid=fluids_parameters.grid,&u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
    PROJECTION_3D<T>& projection=fluids_parameters.incompressible.projection;
    for(int i=1;i<=u_grid.m;i++)for(int j=1;j<=u_grid.n;j++)for(int ij=1;ij<=u_grid.mn;ij++)if(source.Lazy_Inside(u_grid.X(i,j,ij))){projection.elliptic_solver->psi_N_u(i,j,ij)=true;projection.u(i,j,ij)=source_velocity.x;}
    for(int i=1;i<=v_grid.m;i++)for(int j=1;j<=v_grid.n;j++)for(int ij=1;ij<=v_grid.mn;ij++)if(source.Lazy_Inside(v_grid.X(i,j,ij))){projection.elliptic_solver->psi_N_v(i,j,ij)=true;projection.v(i,j,ij)=source_velocity.y;}
    for(int i=1;i<=w_grid.m;i++)for(int j=1;j<=w_grid.n;j++)for(int ij=1;ij<=w_grid.mn;ij++)if(source.Lazy_Inside(w_grid.X(i,j,ij))){projection.elliptic_solver->psi_N_w(i,j,ij)=true;projection.w(i,j,ij)=source_velocity.z;}
    if(fluids_parameters.use_strain)for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)for(int ij=1;ij<=grid.mn;ij++)if(source.Lazy_Inside(grid.X(i,j,ij)))fluids_parameters.e(i,j,ij)=SYMMETRIC_MATRIX<T,3>();
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(!rigid_body_list.rigid_bodies.m) return;
    fluids_parameters.Extrapolate_Velocity_Into_Object(phi_object,V_object,3,true,time);
}
//#####################################################################
// Function Adjust_Strain
//#####################################################################
void Adjust_Strain(ARRAY<SYMMETRIC_MATRIX<T,3> ,VECTOR<int,3> >& e_ghost,const T time)
{
    if(!rigid_body_list.rigid_bodies.m)return;
    fluids_parameters.Adjust_Strain_For_Object(levelset_object,e_ghost,time);
}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Read_Output_Files_Fluids(frame);
    rigid_body_list.template Read<RW>(output_directory,frame);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Write_Output_Files(frame);
    rigid_body_list.template Write<RW>(output_directory,frame);
}
//#####################################################################
};      
}
#endif    


