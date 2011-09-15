#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FALLING_DROP
//#####################################################################
#ifndef __FALLING_DROP__
#define __FALLING_DROP__

#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
namespace PhysBAM{

template<class T_input>
class FALLING_DROP:public SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_3D<T_input> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef RLE_GRID_3D<T> T_GRID;
    typedef SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_3D<T> > BASE;
public:
    using BASE::grid;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::fluids_parameters;
    using BASE::write_frame_title;using BASE::write_substeps;using BASE::ground;using BASE::phi_object;
    using BASE::V_object;using BASE::levelset_object;using BASE::incompressible;
    typedef typename RLE_GRID_3D<T>::ITERATOR_CELL ITERATOR_CELL; typedef typename RLE_GRID_3D<T>::ITERATOR_FACE_X ITERATOR_FACE_X;

    FALLING_DROP()
    {
        grid.Set_Uniform_Grid(GRID<TV>(101,101,101,0,1,0,1,0,1));
        grid.Set_Negative_And_Positive_Bandwidths(5,3);
        //grid.Set_Negative_And_Positive_Bandwidths(5,3);
        //grid.Set_Negative_And_Positive_Bandwidths(100,3);
        grid.Set_Linear_Pressure_And_Constant_Horizontal_Velocity_Profile();
        first_frame=0;last_frame=100;
        frame_rate=24*4;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.incompressible_iterations=50;
        write_frame_title=true;
        fluids_parameters.write_levelset=fluids_parameters.write_velocity=fluids_parameters.write_particles=fluids_parameters.write_removed_positive_particles=true;
        fluids_parameters.write_removed_negative_particles=fluids_parameters.write_debug_data=true;fluids_parameters.write_ghost_values=true;
        fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        //fluids_parameters.scalar_substeps=2.9;
        fluids_parameters.scalar_substeps=(T).9;
        fluids_parameters.cfl=(T).9;

        output_directory="Falling_Drop/output";
    }

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(ARRAY<bool>*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initial_Ground
//#####################################################################
T Initial_Ground(const VECTOR<T,2>& X) PHYSBAM_OVERRIDE
{
    return 0;
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
T Initial_Phi(const VECTOR<T,3>& X) const PHYSBAM_OVERRIDE
{
    static SPHERE<T> sphere((VECTOR<T,3>((T).5,(T).7,(T).5)),(T).2);
    return min(sphere.Signed_Distance(X),X.y-(T).21);
    //return circle.Signed_Distance(X);
    //return BOX_3D<T>(0.3,0.7,0,0.4).Signed_Distance(X);
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const VECTOR<T,3>& X) const PHYSBAM_OVERRIDE
{
    return X.y;
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Construct_Levelsets_For_Objects(time);
    for(ITERATOR_CELL cell(grid,3);cell;cell++){int c=cell.Cell();
        phi_object(c)=min(phi_object(c),grid.uniform_grid.y_minus_half(cell.jmax()));
        if(cell.Long()) phi_object(c+1)=phi_object(c);}
    levelset_object.Compute_Normals();
    levelset_object.Compute_Block_Minimum_And_Maximum();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    //for(ITERATOR_FACE_X face(grid,0);face;face++)
    //  incompressible.V.u(face.Face())=1;
}
//#####################################################################
};    
}
#endif
#endif
