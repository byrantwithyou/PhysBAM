//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Ordinary_Differential_Equations/EULER_STEP_PARTICLES.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include <Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Set_Radius
//#####################################################################
template<class T> inline T VORTEX_PARTICLE_EVOLUTION_3D<T>::
Gaussian_Kernel(const T distance_squared)
{
//    const T normalization=1/(T)pow(2*pi,1.5);
//    return one_over_radius_cubed*normalization*exp(-(T).5*one_over_radius_squared*distance_squared);
    return exp(-distance_squared*4); // returns 1 when distance=0 and almost 0 when distance=1
}
//#####################################################################
// Function Set_Radius
//#####################################################################
/*template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Set_Radius(const T radius_input)
{
    radius=radius_input;radius_squared=sqr(radius);one_over_radius_squared=(T)1/sqr(radius);one_over_radius_cubed=(T)1/cube(radius);scattered_interpolation.Set_Radius_Of_Influence(radius);
}*/
//#####################################################################
// Function Compute_Body_Force
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Compute_Body_Force(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<TV,TV_INT>& force,const T dt,const T time)
{
    ARRAY<T,TV_INT> grid_vorticity_magnitude(grid.Domain_Indices(2),false);
    VORTICITY_UNIFORM<TV>::Vorticity(grid,FACE_LOOKUP_UNIFORM<TV>(face_velocities_ghost),grid_vorticity,grid_vorticity_magnitude);

    if(apply_individual_particle_forces){
        // compute missing vorticity per particle
        VORTICITY_PARTICLES<TV> missing_vorticity_particles;missing_vorticity_particles.Add_Elements(vorticity_particles.Size());
        LINEAR_INTERPOLATION_UNIFORM<TV,TV> vorticity_interpolation;

        for(int p=0;p<vorticity_particles.Size();p++){
            // find missing vorticity
            TV local_grid_vorticity=vorticity_interpolation.Clamped_To_Array(grid,grid_vorticity,vorticity_particles.X(p));
            TV missing_vorticity=vorticity_particles.vorticity(p)-local_grid_vorticity;
            TV sign_check=missing_vorticity*vorticity_particles.vorticity(p);
            for(int a=0;a<3;a++) if(sign_check[a]<0) missing_vorticity[a]=0;
            missing_vorticity_particles.X(p)=vorticity_particles.X(p);
            missing_vorticity_particles.vorticity(p)=missing_vorticity;
            missing_vorticity_particles.radius(p)=vorticity_particles.radius(p);}
        
        if(mpi_grid){
            T max_radius = (T)0;
            if(vorticity_particles.Size()>=1) max_radius=vorticity_particles.radius.Max();
            Exchange_Boundary_Particles_Flat(*mpi_grid,missing_vorticity_particles,max_radius);}

        T small_number=(T)1e-4*grid.dX.Min();
        for(int p=0;p<missing_vorticity_particles.Size();p++){
            T radius=missing_vorticity_particles.radius(p);
            TV_INT box_radius((int)(radius/grid.dX.x)+1,(int)(radius/grid.dX.y)+1,(int)(radius/grid.dX.z)+1);
            TV_INT index=grid.Clamped_Index(missing_vorticity_particles.X(p));
            RANGE<TV_INT> range(clamp_min(index-box_radius,TV_INT()),clamp_max(index+box_radius,grid.counts));
            for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
                TV X_minus_Xp=grid.X(it.index)-missing_vorticity_particles.X(p);
                T distance_squared=X_minus_Xp.Magnitude_Squared();
                if(distance_squared>small_number&&distance_squared<=sqr(radius)){
//                    force(it.index)-=(T)(force_scaling*Gaussian_Kernel(distance_squared)/sqrt(distance_squared))*TV::Cross_Product(X_minus_Xp,missing_vorticity_particles.vorticity(p));}}}
                    T distance=sqrt(distance_squared);
                    force(it.index)-=(T)(force_scaling*Gaussian_Kernel(sqr(distance/radius))/distance)*TV::Cross_Product(X_minus_Xp,missing_vorticity_particles.vorticity(p));}}}}
    else{
        if(mpi_grid) PHYSBAM_NOT_IMPLEMENTED(); // this has not been mpi'd yet
        ARRAY<T,TV_INT> grid_vorticity_particles_magnitude(grid.Domain_Indices(2),false);
    
        // form grid vorticity from vortex particles
        scattered_interpolation.Transfer_To_Grid(vorticity_particles.X,vorticity_particles.vorticity,grid,grid_vorticity_particles);
        if(remove_grid_vorticity_from_particle_vorticity)
            for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next())
                grid_vorticity_particles(it.index)-=grid_vorticity(it.index);
    
        // find vorticity magnitudes
        for(RANGE_ITERATOR<TV::m> it(grid_vorticity.domain);it.Valid();it.Next())
            grid_vorticity_magnitude(it.index)=grid_vorticity(it.index).Magnitude();
        for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next())
            grid_vorticity_particles_magnitude(it.index)=grid_vorticity_particles(it.index).Magnitude();
    
        // compute confinement force
        for(RANGE_ITERATOR<TV::m> it(force.domain);it.Valid();it.Next()){
            TV vortex_normal_vector,particle_vortex_normal_vector;
            for(int d=0;d<TV::m;d++){
                TV_INT a(it.index),b(it.index);
                a(d)++;
                b(d)--;
                vortex_normal_vector(d)=(grid_vorticity_magnitude(a)-grid_vorticity_magnitude(b))*(T).5*grid.one_over_dX(d);
                particle_vortex_normal_vector(d)=(grid_vorticity_particles_magnitude(a)-grid_vorticity_particles_magnitude(b))*(T).5*grid.one_over_dX(d);}
            vortex_normal_vector.Normalize();
            particle_vortex_normal_vector.Normalize();
            force(it.index)=grid_confinement_parameter*TV::Cross_Product(vortex_normal_vector,grid_vorticity(it.index))
                        +particle_confinement_parameter*TV::Cross_Product(particle_vortex_normal_vector,grid_vorticity_particles(it.index));}}
}
//#####################################################################
// Function Compute_Body_Force
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Compute_Body_Force(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time)
{
    ARRAY<TV,TV_INT> cell_force(grid.Domain_Indices(1),false);
    Compute_Body_Force(face_velocities_ghost,cell_force,dt,time);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        force.Component(axis)(iterator.Face_Index())+=(T).5*(cell_force(iterator.First_Cell_Index())[axis]+cell_force(iterator.Second_Cell_Index())[axis]);}
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Euler_Step(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const T dt,const T time)
{
    LOG::Time("Advancing vorticity particles");
    
    ARRAY<TV,TV_INT> two_times_V_ghost(grid.Domain_Indices(2));
    for(CELL_ITERATOR<TV> it(grid,2);it.Valid();it.Next())
        for(int d=0;d<TV::m;d++){
            TV_INT next(it.index);
            next(d)++;
            two_times_V_ghost(it.index)(d)=face_velocities_ghost.Component(d)(it.index)+face_velocities_ghost.Component(d)(next);}

    // vortex stretching/tilting term  - omega dot grad V
    ARRAY<MATRIX<T,3> ,TV_INT> VX(grid.Domain_Indices(1),false);LINEAR_INTERPOLATION_UNIFORM<TV,MATRIX<T,3> > VX_interpolation;
    T one_over_four_dx=1/(4*grid.dX.x),one_over_four_dy=1/(4*grid.dX.y),one_over_four_dz=1/(4*grid.dX.z);
    for(int i=0;i<grid.counts.x+1;i++) for(int j=0;j<grid.counts.y+1;j++) for(int ij=0;ij<grid.counts.z+1;ij++)
        VX(TV_INT(i,j,ij))=MATRIX<T,3>(one_over_four_dx*(two_times_V_ghost(i+1,j,ij)-two_times_V_ghost(i-1,j,ij)),
                                 one_over_four_dy*(two_times_V_ghost(i,j+1,ij)-two_times_V_ghost(i,j-1,ij)),
                                 one_over_four_dz*(two_times_V_ghost(i,j,ij+1)-two_times_V_ghost(i,j,ij-1)));

    if(renormalize_vorticity_after_stretching_tilting) for(int p=0;p<vorticity_particles.Size();p++){
        T old_magnitude=vorticity_particles.vorticity(p).Magnitude();
        vorticity_particles.vorticity(p)+=dt*VX_interpolation.Clamped_To_Array(grid,VX,vorticity_particles.X(p))*vorticity_particles.vorticity(p);
        vorticity_particles.vorticity(p).Normalize();vorticity_particles.vorticity(p)*=old_magnitude;}
    else for(int p=0;p<vorticity_particles.Size();p++){
        vorticity_particles.vorticity(p)+=dt*VX_interpolation.Clamped_To_Array(grid,VX,vorticity_particles.X(p))*vorticity_particles.vorticity(p);}
    // advance vortex particle position
    EULER_STEP_PARTICLES<TV>::Euler_Step_Face(vorticity_particles.X,grid,face_velocities_ghost,dt);
    if(mpi_grid) Exchange_Boundary_Particles_Flat(*mpi_grid,vorticity_particles);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
   LOG::Time("Writing Vortex Specific Data");
   FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/vorticity_particles",output_directory.c_str(),frame),vorticity_particles);
   FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/grid_vorticity",output_directory.c_str(),frame),grid_vorticity);
   FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/grid_vorticity_particles",output_directory.c_str(),frame),grid_vorticity_particles);
   LOG::Stop_Time();
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& input_directory,const int frame)
{
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/vorticity_particles",input_directory.c_str(),frame),vorticity_particles);
}
//#####################################################################
namespace PhysBAM{
template class VORTEX_PARTICLE_EVOLUTION_3D<float>;
template class VORTEX_PARTICLE_EVOLUTION_3D<double>;
}
