//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Level_Sets/REINITIALIZATION.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Forces/FLUID_GRAVITY.h>
#include <Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV_input> PLS_FC_EXAMPLE<TV_input>::
PLS_FC_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),last_frame(100),write_substeps_level(-1),substeps_delay_frame(-1),
    write_output_files(true),viewer_dir("output"),restart(0),number_of_ghost_cells(5),dt(1),time(0),
    time_steps_per_frame(1),use_preconditioner(true),max_iter(100000),solver_tolerance(1e-10),dump_matrix(false),
    sparse_dump_matrix(false),use_advection(true),use_reduced_advection(true),omit_solve(false),number_of_colors(1),
    use_discontinuous_velocity(false),use_p_null_mode(false),use_level_set_method(false),use_pls(false),
    dump_largest_eigenvector(false),save_pressure(false),test_system(false),use_polymer_stress(false),
    num_multigrid_levels(2),use_multigrid(false),grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    collision_bodies_affecting_fluid(*new GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>(grid)),
    particle_levelset_evolution_multiple(*new PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>(grid,collision_bodies_affecting_fluid,number_of_ghost_cells)),
    advection_scalar(*new ADVECTION_HAMILTON_JACOBI_ENO<TV,T>),
    levelset_color(grid,*new ARRAY<T,TV_INT>,*new ARRAY<int,TV_INT>),debug_particles(*new DEBUG_PARTICLES<TV>)
{
    debug_particles.debug_particles.template Add_Array<TV>("V");
    debug_particles.debug_particles.template Add_Array<T>("display_size");
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV_input> PLS_FC_EXAMPLE<TV_input>::
~PLS_FC_EXAMPLE()
{
    delete &debug_particles;
    delete &collision_bodies_affecting_fluid;
    delete &levelset_color.color;
    delete &levelset_color.phi;
    delete &advection_scalar;
    delete &particle_levelset_evolution_multiple;
}
//#####################################################################
// Function Merge_Velocities
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Merge_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& V,const ARRAY<ARRAY<T,FACE_INDEX<TV::m> > > u,const ARRAY<int,FACE_INDEX<TV::m> >& color) const
{
    for(FACE_ITERATOR<TV> it(grid,number_of_ghost_cells);it.Valid();it.Next()){
        int c=color(it.Full_Index());
        if(c<0){
            c=0;
            if(abs(u(c)(it.Full_Index()))>1e10) continue;}
        V(it.Full_Index())=u(c)(it.Full_Index());}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Write_Output_Files()
{
    if(!write_output_files) return;
    ARRAY<T,FACE_INDEX<TV::m> > V(face_velocities(0).domain_indices);
    Merge_Velocities(V,face_velocities,face_color);

    Write_To_File(stream_type,viewer_dir.current_directory+"/mac_velocities",V);
    Write_To_File(stream_type,viewer_dir.output_directory+"/common/grid",grid);
    // particle levelset
    for(int i=0;i<number_of_colors;i++){
        std::string ii=LOG::sprintf("%d",i);
        PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=*particle_levelset_evolution_multiple.particle_levelset_multiple.particle_levelsets(i);
        Write_To_File(stream_type,viewer_dir.current_directory+"/levelset_"+ii,particle_levelset.levelset);
        Write_To_File(stream_type,viewer_dir.current_directory+"/positive_particles_"+ii,particle_levelset.positive_particles);
        Write_To_File(stream_type,viewer_dir.current_directory+"/negative_particles_"+ii,particle_levelset.negative_particles);
        Write_To_File(stream_type,viewer_dir.current_directory+"/removed_positive_particles_"+ii,particle_levelset.removed_positive_particles);
        Write_To_File(stream_type,viewer_dir.current_directory+"/removed_negative_particles_"+ii,particle_levelset.removed_negative_particles);
        Write_To_Text_File(viewer_dir.current_directory+"/last_unique_particle_id_"+ii,particle_levelset.last_unique_particle_id);}
    debug_particles.Write_Debug_Particles(stream_type,viewer_dir);
    if(save_pressure) Write_To_File(stream_type,viewer_dir.current_directory+"/pressure",pressure);
    Write_To_File(stream_type,viewer_dir.current_directory+"/restart_data",
        time,face_color,prev_face_color,face_velocities,prev_face_velocities);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Read_Output_Files()
{
    for(int i=0;i<number_of_colors;i++){
        std::string ii=LOG::sprintf("%d",i);
        PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=*particle_levelset_evolution_multiple.particle_levelset_multiple.particle_levelsets(i);
        Read_From_File(viewer_dir.current_directory+"/levelset_"+ii,particle_levelset.levelset);
        Read_From_File(viewer_dir.current_directory+"/positive_particles_"+ii,particle_levelset.positive_particles);
        Read_From_File(viewer_dir.current_directory+"/negative_particles_"+ii,particle_levelset.negative_particles);
        Read_From_File(viewer_dir.current_directory+"/removed_positive_particles_"+ii,particle_levelset.removed_positive_particles);
        Read_From_File(viewer_dir.current_directory+"/removed_negative_particles_"+ii,particle_levelset.removed_negative_particles);
        Read_From_Text_File(viewer_dir.current_directory+"/last_unique_particle_id_"+ii,particle_levelset.last_unique_particle_id);}
    if(save_pressure) Read_From_File(viewer_dir.current_directory+"/pressure",pressure);
    Read_From_File(viewer_dir.current_directory+"/restart_data",
        time,face_color,prev_face_color,face_velocities,prev_face_velocities);
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    for(int i=0;i<number_of_colors;i++){
        TV& X=particles.X(index);TV X_new=X+dt*V;
        TV min_corner=grid.domain.Minimum_Corner(),max_corner=grid.domain.Maximum_Corner();
        for(int axis=0;axis<TV::m;axis++){
            T shift=((X_new[axis]-min_corner[axis])/(max_corner[axis]-min_corner[axis]));if(shift<0)shift-=(T)1;X_new[axis]-=(max_corner[axis]-min_corner[axis])*(int)shift;X=X_new-dt*V;
        }
    }
}
//#####################################################################
// Function Color_At_Cell
//#####################################################################
template<class TV_input> int PLS_FC_EXAMPLE<TV_input>::
Color_At_Cell(const TV_INT& index) const
{
    for(int i=0;i<bc_phis.m;i++)
        if(bc_phis(i)(index)<=0)
            return ~i;
    return particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.Inside_Region(index);
}
//#####################################################################
// Function Color_At_Cell
//#####################################################################
template<class TV_input> int PLS_FC_EXAMPLE<TV_input>::
Color_At_Cell(const TV_INT& index,T& phi) const
{
    for(int i=0;i<bc_phis.m;i++)
        if(bc_phis(i)(index)<=0){
            phi=-bc_phis(i)(index);
            return ~i;}
    int c=particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.Inside_Region(index,phi);
    phi=-phi;
    return c;
}
//#####################################################################
// Function Enforce_Phi_Boundary_Conditions
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Enforce_Phi_Boundary_Conditions()
{
    ARRAY<ARRAY<T,TV_INT> >& phis=particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.phis;
    for(int i=0;i<phis.m;i++){
        boundary.Fill_Ghost_Cells(grid,phis(i),phis(i),0,0,number_of_ghost_cells);
        if(phis(i).array.Min()>=0) phis(i).array.Fill(number_of_ghost_cells*grid.dX.Max());}
}
//#####################################################################
// Function Rebuild_Levelset_Color
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Rebuild_Levelset_Color()
{
    Make_Levelsets_Consistent();
    for(CELL_ITERATOR<TV> it(grid,number_of_ghost_cells);it.Valid();it.Next())
        levelset_color.color(it.index)=Color_At_Cell(it.index,levelset_color.phi(it.index));
    boundary.Fill_Ghost_Cells(grid,levelset_color.phi,levelset_color.phi,0,0,number_of_ghost_cells);
    boundary_int.Fill_Ghost_Cells(grid,levelset_color.color,levelset_color.color,0,0,number_of_ghost_cells);
}
//#####################################################################
// Function Fill_Levelsets_From_Levelset_Color
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Fill_Levelsets_From_Levelset_Color()
{
    ARRAY<ARRAY<T,TV_INT> >& phis=particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.phis;
    for(CELL_ITERATOR<TV> it(grid,number_of_ghost_cells);it.Valid();it.Next()){
        int c=levelset_color.color(it.index);
        T p=levelset_color.phi(it.index);
        for(int i=0;i<bc_phis.m;i++)
            bc_phis(i)(it.index)=(c==~i)?-p:p;
        for(int i=0;i<phis.m;i++)
            phis(i)(it.index)=(c==i)?-p:p;}

    for(int i=0;i<bc_phis.m;i++){
        LEVELSET<TV> fl(grid,bc_phis(i),number_of_ghost_cells);
        fl.boundary=&boundary;
        Reinitialize(fl,number_of_ghost_cells*2,(T)0,number_of_ghost_cells*grid.dX.Max(),grid.domain.Edge_Lengths().Magnitude(),(T).9,3,5,1);}
    for(int i=0;i<phis.m;i++)
        Reinitialize(*particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.levelsets(i),number_of_ghost_cells*2,(T)0,
            number_of_ghost_cells*grid.dX.Max(),grid.domain.Edge_Lengths().Magnitude(),(T).9,3,5,1);

    for(int i=0;i<bc_phis.m;i++){
        boundary.Fill_Ghost_Cells(grid,bc_phis(i),bc_phis(i),0,0,number_of_ghost_cells);
        if(bc_phis(i).array.Min()>=0) bc_phis(i).array.Fill(number_of_ghost_cells*grid.dX.Max());}
    for(int i=0;i<phis.m;i++){
        boundary.Fill_Ghost_Cells(grid,phis(i),phis(i),0,0,number_of_ghost_cells);
        if(phis(i).array.Min()>=0) phis(i).array.Fill(number_of_ghost_cells*grid.dX.Max());}
}
//#####################################################################
// Function Get_Levelset_Velocity
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::m> >& V_levelset,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Get_Levelset_Velocity
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<TV>& levelset_multiple,ARRAY<T,FACE_INDEX<TV::m> >& V_levelset,const T time) const
{
    ARRAY<T,FACE_INDEX<TV::m> > n(grid,number_of_ghost_cells),m(grid,number_of_ghost_cells);
    Merge_Velocities(n,face_velocities,face_color);
    Merge_Velocities(m,prev_face_velocities,prev_face_color);
    T alpha=(time-this->time)/dt;
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next())
        V_levelset(it.Full_Index())=(1-alpha)*n(it.Full_Index())+alpha*m(it.Full_Index());
    boundary.Apply_Boundary_Condition_Face(grid,V_levelset,time+dt);
}
//#####################################################################
// Function Make_Levelsets_Consistent
//#####################################################################
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Make_Levelsets_Consistent()
{
    ARRAY<ARRAY<T,TV_INT> >& phis=particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.phis;
    for(CELL_ITERATOR<TV> it(grid,number_of_ghost_cells);it.Valid();it.Next()){
        T min1=FLT_MAX,min2=FLT_MAX,bc_min=FLT_MAX;
        for(int i=0;i<bc_phis.m;i++)
            bc_min=min(bc_min,bc_phis(i)(it.index));
        for(int i=0;i<phis.m;i++){
            T p=phis(i)(it.index);
            if(p<min1){min2=min1;min1=p;}
            else if(p<min2) min2=p;}
        T shift=min(bc_min+min1,(T).5*(min2+min1));
        for(int i=0;i<phis.m;i++)
            phis(i)(it.index)-=shift;}
    Enforce_Phi_Boundary_Conditions();
}
//#####################################################################
namespace PhysBAM{
template class PLS_FC_EXAMPLE<VECTOR<float,2> >;
template class PLS_FC_EXAMPLE<VECTOR<float,3> >;
template class PLS_FC_EXAMPLE<VECTOR<double,2> >;
template class PLS_FC_EXAMPLE<VECTOR<double,3> >;
}
