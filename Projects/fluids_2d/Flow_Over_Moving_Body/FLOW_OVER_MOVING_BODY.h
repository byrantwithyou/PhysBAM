//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOW_OVER_MOVING_BODY 
//##################################################################### 
#ifndef __FLOW_OVER_MOVING_BODY__
#define __FLOW_OVER_MOVING_BODY__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template<class T,class RW=T>
class FLOW_OVER_MOVING_BODY:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::abort_when_dt_below;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Use_Thin_Shells_Fluid_Coupling_Defaults;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::use_object_pseudo_velocity_for_boundary_conditions;

    bool use_source;
    BOX_2D<T> source;
    MATRIX<T,3> world_to_source;
    VECTOR_2D<T> source_velocity;

    ARRAY<T,VECTOR<int,2> > phi_object;
    bool write_phi_object;
    RANDOM_NUMBERS random;
    T initial_water_level;

    FLOW_OVER_MOVING_BODY()
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.WATER)
    {
        first_frame=0;last_frame=250;
        frame_rate=24;
        restart=false;restart_frame=0;
        fluids_parameters.grid.Initialize(121,121,-3,3,0,6);
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.bias_towards_negative_particles=true;fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;

        initial_water_level=(T)0.021; // some number so the interface won't lie in the center of a cell

        use_source=true;
        source=BOX_2D<T>(-0.6,0.6,-1.2,0);
        MATRIX<T,3> source_to_world=MATRIX<T,3>::Translation_Matrix(VECTOR_2D<T>(-1.5,6));
        world_to_source=source_to_world.Inverse();
        source_velocity=VECTOR_2D<T>(0,-1);

        fluids_parameters.viscosity=0;
        fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=false;fluids_parameters.second_order_cut_cell_method=false;

        fluids_parameters.incompressible_iterations=100;

        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;fluids_parameters.write_debug_data=true;
        fluids_parameters.store_particle_ids=true;
        output_directory="Flow_Over_Moving_Body/output";
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;

        write_phi_object=true;

        solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/circle", (T)1, true, true);
        Update_Rigid_Bodies(0);
    }

    ~FLOW_OVER_MOVING_BODY() 
    {}

//#####################################################################
// Function Update_Rigid_Bodies
//#####################################################################
void Update_Rigid_Bodies(const T time)
{
    INTERPOLATION_CURVE<T,VECTOR_2D<T> > motion_curve;
    motion_curve.Add_Control_Point(0,VECTOR_2D<T>(-1.5,2));
    motion_curve.Add_Control_Point(1,VECTOR_2D<T>(-1.5,2));
    motion_curve.Add_Control_Point(2,VECTOR_2D<T>(-1.5,4.5));
    motion_curve.Add_Control_Point(3,VECTOR_2D<T>(1.5,2));
    for(int i=0;i<solids_parameters.rigid_body_parameters.list.rigid_bodies.m;i++){
        solids_parameters.rigid_body_parameters.list.rigid_bodies(i)->position=motion_curve.Value(time);
        solids_parameters.rigid_body_parameters.list.rigid_bodies(i)->velocity=motion_curve.Derivative(time);
    }
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    GRID<TV>& grid=fluids_parameters.grid;
    // Not so good to set up a heaviside function here because then the interface will
    // be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++)
        fluids_parameters.particle_levelset_evolution.phi(i,j)=grid.y(j)-grid.ymin-initial_water_level;
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    if(solids_parameters.rigid_body_parameters.list.rigid_bodies.m<1) return;
    Update_Rigid_Bodies(time);
    phi_object.Resize(grid,3);
    for(int i=-2;i<=grid.m+3;i++) for(int j=-2;j<=grid.n+3;j++)
        phi_object(i,j)=solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Implicit_Curve_Extended_Value(grid.X(i,j));

    if(write_phi_object){
        ARRAY<T,VECTOR<int,2> > temp_phi(grid);
        ARRAY<T,VECTOR<int,2> >::get(temp_phi,phi_object);
        LEVELSET_2D<T> temp_levelset(grid,temp_phi);
        char filename[256];sprintf(filename,"%s/phi_object.phi",output_directory.c_str());
        std::ofstream output(filename,std::ios::binary);
        temp_levelset.template Write<RW>(output);
    }
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    if(!use_source) return;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){
        VECTOR_2D<T> source_X=world_to_source*grid.X(i,j);
        if(source.Lazy_Inside(source_X)) fluids_parameters.particle_levelset_evolution.phi(i,j)=min(fluids_parameters.particle_levelset_evolution.phi(i,j),source.Signed_Distance(source_X));}
}
//#####################################################################
// Function Adjust_Phi_with_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;INCOMPRESSIBLE_2D<T>& incompressible=fluids_parameters.incompressible;
    if(solids_parameters.rigid_body_parameters.list.rigid_bodies.m<1) return;
    T one_over_two_dx=1/(2*grid.dx),one_over_two_dy=1/(2*grid.dy);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) if(fluids_parameters.particle_levelset_evolution.phi(i,j) < 0 && phi_object(i,j) < 0){
        VECTOR_2D<T> V_fluid=incompressible.V(i,j);
        VECTOR_2D<T> V_object=solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Pointwise_Object_Velocity(grid.X(i,j)); // velocity object should be spatially varying
        VECTOR_2D<T> V_relative=V_fluid-V_object;
        VECTOR_2D<T> normal=VECTOR_2D<T>((phi_object(i+1,j)-phi_object(i-1,j))*one_over_two_dx,
                                         (phi_object(i,j+1)-phi_object(i,j-1))*one_over_two_dy);
        T denominator=normal.Magnitude();if(denominator > 1e-8) normal/=denominator;else normal=VECTOR_2D<T>(1,0);
        T VN=VECTOR_2D<T>::Dot_Product(V_relative,normal),magnitude=V_relative.Magnitude();
        if(VN > .1*magnitude) fluids_parameters.particle_levelset_evolution.phi(i,j)=-phi_object(i,j);
    }
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects(const T time)
{
    if(solids_parameters.rigid_body_parameters.list.rigid_bodies.m<1) return;
    fluids_parameters.Extrapolate_Phi_Into_Object(phi_object);
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    if(!use_source) return;
    if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,2> >(grid);

    T padding=3*grid.max_dx_dy;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) if(!source.Outside(world_to_source*grid.X(i,j),padding)) (*cell_centered_mask)(i,j)=true;
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(VECTOR_2D<T>& X,VECTOR_2D<T>& V,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{          
    if(solids_parameters.rigid_body_parameters.list.rigid_bodies.m<1) return true;
    fluids_parameters.Adjust_Particle_For_Object(*solids_parameters.rigid_body_parameters.list(1),X,V,particle_type,dt,time);
    return true;
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR_2D<T> >& particles,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,
                                     const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    if(solids_parameters.rigid_body_parameters.list.rigid_bodies.m<1) return;
    if(particle_type == PARTICLE_LEVELSET<T,VECTOR_2D<T> >::NEGATIVE || particle_type == PARTICLE_LEVELSET<T,VECTOR_2D<T> >::REMOVED_NEGATIVE){
        for(int k=0;k<particles.Size();k++) if(solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Implicit_Curve_Lazy_Inside_Extended_Levelset(particles.X(k),-grid.dx)) particles.Add_To_Deletion_List(k);}
   else for(int k=0;k<particles.Size();k++) if(solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Implicit_Curve_Lazy_Inside_Extended_Levelset(particles.X(k))) particles.Add_To_Deletion_List(k);
   particles.Delete_Elements_On_Deletion_List(false,true); // already sorted
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    GRID<TV> &u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid;
    PROJECTION_2D<T>& projection=fluids_parameters.incompressible.projection;
    if(!use_source) return;
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) if(source.Lazy_Inside(world_to_source*u_grid.X(i,j))){projection.elliptic_solver->psi_N_u(i,j)=true;projection.u(i,j)=source_velocity.x;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) if(source.Lazy_Inside(world_to_source*v_grid.X(i,j))){projection.elliptic_solver->psi_N_v(i,j)=true;projection.v(i,j)=source_velocity.y;}
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(solids_parameters.rigid_body_parameters.list.rigid_bodies.m<1) return;
    fluids_parameters.Extrapolate_Velocity_Into_Object(phi_object,*solids_parameters.rigid_body_parameters.list.rigid_bodies(1),3,true,time);
}
//#####################################################################
};    
}
#endif  
