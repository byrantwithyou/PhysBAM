#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHEAR_TEST
//#####################################################################
#ifndef __SHEAR_TEST__
#define __SHEAR_TEST__

#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
namespace PhysBAM{

template<class T>
class SHEAR_TEST:public SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_2D<T> >
{
    typedef RLE_GRID_2D<T> T_GRID;typedef VECTOR<T,2> TV;
    typedef SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID> BASE;
public:
    using BASE::grid;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;
    using BASE::fluids_parameters;using BASE::write_frame_title;using BASE::write_substeps;using BASE::ground;using BASE::phi_object;using BASE::V_object;using BASE::levelset_object;
    using BASE::incompressible;using BASE::particle_levelset;

    T water_level;
    T shear_constant;
    ARRAY<VECTOR<int,2> > long_run_range;

    SHEAR_TEST(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_2D<T> >(stream_type)
    {
        grid.Set_Uniform_Grid(GRID<TV>(101,101,0,1,0,1));
        grid.Set_Negative_And_Positive_Bandwidths(5,5);
        grid.Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile();
        first_frame=0;last_frame=2000;
        frame_rate=24*4;
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=true; 
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.incompressible_iterations=50;
        fluids_parameters.write_levelset=fluids_parameters.write_velocity=fluids_parameters.write_particles=fluids_parameters.write_removed_positive_particles=true;
        fluids_parameters.write_removed_negative_particles=fluids_parameters.write_debug_data=true;fluids_parameters.write_ghost_values=true;
        fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.store_particle_ids=true;

        fluids_parameters.gravity=0;
        fluids_parameters.object_friction=1;
        shear_constant=2;
        water_level=0.5103;
        long_run_range.Resize(grid.uniform_grid.m-1);
        RANDOM_NUMBERS<T> random;random.Set_Seed(12345);
        int max_j=(int)(water_level*(grid.uniform_grid.n-1)-grid.negative_bandwidth);
        for(int i=1;i<=grid.uniform_grid.m-1;i++){
            long_run_range(i).x=random.Get_Uniform_Integer(1,max_j/3);
            long_run_range(i).y=random.Get_Uniform_Integer(2*max_j/3,max_j);
            LOG::cout << i << ": " << long_run_range(i).x << ", " << long_run_range(i).y << std::endl;}

        output_directory="Shear_Test/output";
    }

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE {return false;}
    void Get_Source_Reseed_Mask(ARRAY<bool>*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initial_Ground
//#####################################################################
T Initial_Ground(const VECTOR<T,1>& X) PHYSBAM_OVERRIDE
{
    return 0;
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
T Initial_Phi(const VECTOR<T,2>& X) const PHYSBAM_OVERRIDE
{
    return X.y-water_level;
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(typename RLE_GRID_2D<T>::ITERATOR_FACE_X face(grid,0);face;face++){
        if(grid.long_run_faces_horizontal==1) incompressible.V.u(face.Face())=shear_constant*face.X().y;
        else{
            incompressible.V.u(face.Face())=shear_constant*grid.uniform_grid.y(face.j());
            if(face.Long()) incompressible.V.u(face.Face()+1)=shear_constant*grid.uniform_grid.y(face.jmax()-1);}}
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Get_Object_Velocities(dt,time);

    for(typename RLE_GRID_2D<T>::ITERATOR_FACE_X face(grid,BOX<VECTOR<int,1> >(1-grid.number_of_ghost_cells,1));face;face++){
        if(grid.long_run_faces_horizontal==1) incompressible.V.u(face.Face())=shear_constant*face.X().y;
        else{
            incompressible.V.u(face.Face())=shear_constant*grid.uniform_grid.y(face.j());
            if(face.Long()) incompressible.V.u(face.Face()+1)=shear_constant*grid.uniform_grid.y(face.jmax()-1);}}

    for(typename RLE_GRID_2D<T>::ITERATOR_FACE_X face(grid,BOX<VECTOR<int,1> >(grid.uniform_grid.m,grid.uniform_grid.m+grid.number_of_ghost_cells));face;face++){
        if(grid.long_run_faces_horizontal==1) incompressible.V.u(face.Face())=shear_constant*face.X().y;
        else{
            incompressible.V.u(face.Face())=shear_constant*grid.uniform_grid.y(face.j());
            if(face.Long()) incompressible.V.u(face.Face()+1)=shear_constant*grid.uniform_grid.y(face.jmax()-1);}}
}
//#####################################################################
// Function Get_Cell_Should_Be_Long
//#####################################################################
void Get_Cell_Should_Be_Long(ARRAY<bool>& cell_should_be_long,const T time) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Get_Cell_Should_Be_Long(cell_should_be_long,time);
    for(typename RLE_GRID_2D<T>::ITERATOR_CELL cell(grid,0);cell;cell++) if(cell.Center().y<=water_level && cell.Center().y>=0){
        if(cell.j >= long_run_range(cell.i).x && cell.jmax() <= long_run_range(cell.i).y){cell_should_be_long(cell.Cell())=true;if(cell.Long()) cell_should_be_long(cell.Cell()+1)=true;}
        else if(cell.jmax() < long_run_range(cell.i).x || cell.j > long_run_range(cell.i).y){cell_should_be_long(cell.Cell())=false;if(cell.Long()) cell_should_be_long(cell.Cell()+1)=false;}}
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const VECTOR<T,2>& X) const PHYSBAM_OVERRIDE
{
    return X.y;
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Construct_Levelsets_For_Objects(time);
    for(typename RLE_GRID_2D<T>::ITERATOR_CELL cell(grid,3);cell;cell++){int c=cell.Cell();
        phi_object(c)=min(phi_object(c),grid.uniform_grid.y_minus_half(cell.jmax()));
        if(cell.Long()) phi_object(c+1)=phi_object(c);}
    levelset_object.Compute_Normals();
    levelset_object.Compute_Block_Minimum_And_Maximum();
}
//#####################################################################
};    
}
#endif
#endif
